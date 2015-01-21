#include "Equations_conf.h"  // <--- Equations configure
#ifdef TBK_EQUATIONS

// ============================================
#if (TBK_EQUATIONS == 0) // one equation model
// ============================================

// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverTBK.h"       // Navier-Stokes class header file


// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include -----------------------------------
// #include "MGMesh.h"          // Mesh class
#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsSystem.h"  // Equation map class
#include "EquationSystemsExtendedM.h"  // Equation map class

// standard lib -----------------------------------------------
#include <string.h>          // string library

// local alg lib -----------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ==============================================================
#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "InterfaceFunctionM.h"
#endif


// ==========================================================
/// This routine constructs the FSI class:
MGSolTBK::MGSolTBK(
  MGEquationsSystem& mg_equations_map_in,
  int             nvars_in[],         // KLQ number of variables
  std::string eqname_in,
  std::string varname_in
):  MGSolDA(mg_equations_map_in,nvars_in,eqname_in,varname_in),
  /// A) reading class parameters
  // mesh params ------------
  _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes (top level)
  // phys params ------------
  _dt(_mgutils.get_par("dt")),       // parameter  dt
  _uref(mg_equations_map_in.get_par("Uref")),    // parameter  u reference
  _lref(mg_equations_map_in.get_par("Lref")),    // parameter  l reference
  _rhof(mg_equations_map_in.get_par("rho0")),    // parameter density
  _muf(mg_equations_map_in.get_par("mu0")) {     // parameter viscosity
  // ================================================

  /// B) setting class variables
  _dir=0;
  
#if (TBK_EQUATIONS==0)
  _var_names[0]="kt";    _refvalue[0]=_uref; // kappa only equation
#else
#if (TBK_EQUATIONS%2==0)
  if(!varname_in.compare("kt")) _dir=0;  // kappa
  if(!varname_in.compare("et")) _dir=1;  // omega
  _var_names[0]=varname_in;   _refvalue[0]=_uref;
#else
  _var_names[0]="kt";    _refvalue[0]=_uref; // kappa
  _var_names[1]="et";    _refvalue[1]=_uref; // omega
#endif
#endif
  /// C ) setting solver type
  for(int l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVER_TBK);
  /// D) setting nondimensional parameters
  _IRe=_muf/(_rhof*_lref*_uref);          // Reynolds number

  return;
}





// ====================================================================
/// This function assembles the matrix and the rhs:
//  ======================================================================
  void  MGSolTBK::GenMatRhs(const double time, const int
                            Level,const  int mode) {
    // ===============================================
    /// a) Set up
    // geometry -----
    const int    ndim = DIMENSION;                                           //dimension
    double      xx_qnds[DIMENSION*NDOF_FEM], xxb_qnds[DIMENSION*NDOF_FEMB]; // element node coords
    int         el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];                     // element connectivity
    int         el_neigh[NDOF_FEM];                                         // element connectivity
    const int   offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
    const int   el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
    int         sur_toply[NDOF_FEMB];                                       // boundary topology

    // gauss integration ------------------
    const int    el_ngauss = _fe[2]->_NoGauss1[ndim-1];   //elem gauss points
    const int    elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];   //elem gauss points
    double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];            // Jac, Jac*w Jacobean
    double dphijdx_g[3][DIMENSION];
    double dphiidx_g[3][DIMENSION]; // global derivatives at g point

    //number of constant[0]-linear[1]-quadratic[2] element dof
    int el_ndof[3];  el_ndof[0]=1;  int elb_ndof[3]; elb_ndof[0]=1;
    int el_mat_nrows =0;
    for(int ideg=1; ideg<3; ideg++) {
      el_ndof[ideg]=_fe[ideg]->_NoShape[ndim-1];
      elb_ndof[ideg]=_fe[ideg]->_NoShape[ ndim-2];
      el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
    };
    int el_mat_ncols = el_mat_nrows;                     //square matrix
    std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector

    // fields -> Navier-Stokes
    const int ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NS_F]];
    const int k_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[K_F]];  // kappa-omega equations
    double    normal[DIMENSION]; double val_tbg[10];
    double    vel_g[DIMENSION]; double vel_dxg[DIMENSION*DIMENSION]; double mod2_vel=1.e-20;
    double Pe_h[DIMENSION]; double h_eff[DIMENSION]; double f_upwind[DIMENSION];
    for(int idim=0; idim< DIMENSION; idim++) {
      vel_g[idim]=(1-idim)*0.+idim*.7;
      mod2_vel += vel_g[idim]*vel_g[idim];
      for(int jdim=0; jdim< DIMENSION; jdim++) vel_dxg[idim+jdim*DIMENSION]=0.;
    }
    double IRe_eff[2];
    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) -----------
    A[Level]->zero();    if(mode ==1) b[Level]->zero();
    DenseMatrixM KeM;    DenseVectorM FeM;
    KeM.resize(el_mat_nrows,el_mat_ncols);    FeM.resize(el_mat_nrows);

    int ndof_lev=0;
        for (int pr=0;pr <_mgmesh._iproc; pr++) {
            int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
            ndof_lev +=delta;
        }
    
    /// b) Element  Loop over the volume (n_elem)
    //element type set up
    const int  nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1];
    const int  nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];
    for(int iel=0; iel < (nel_e - nel_b); iel++) {

      // set to zero matrix and rhs
      KeM.zero();         FeM.zero();

      // geometry element quantities ---------------------------
      // Element connectivity  and coordinates (xx_qnds)
      _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);
      _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);

      // set element-nodes variables  bc (bc_q_dofs)
      get_el_dof_bc(Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

      // element nodes coordinates
      for (int idim=0; idim<DIMENSION; idim++) {
        for (int d=0; d< NDOF_FEM; d++)  _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];
        // element grid distance
        h_eff[idim]=fabs(xx_qnds[idim*NDOF_FEM+2+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+0]);
        double h_1=fabs(xx_qnds[idim*NDOF_FEM+3+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+1]);
        if (h_eff[idim] <  h_1) h_eff[idim]=h_1; // Max dx diagonal term
      }
      

      // external fields
      for(int eq=0; eq<_data_eq[2].n_eqs; eq++) {
        _data_eq[2].mg_eqs[eq]->
        get_el_sol(0,_data_eq[2].indx_ub[eq+1]-_data_eq[2].indx_ub[eq],
                   el_ndof[2],el_conn,offset,_data_eq[2].indx_ub[eq],
                   _data_eq[2].ub);
      }

#ifdef DIST_FIX
      _y_dist=DIST_FIX;
#else
      _y_dist=_mgmesh._dist[ iel+nel_b];
#endif


      // element fields ----------------------------------
      double phase=1.;    double rho =  phase;
     
      // *****************************************************************
      /// c) gaussian integration loop (n_gauss)
      for(int qp=0; qp<  el_ngauss; qp++) {

        // shape functions at gaussian points -----------------------------------
        for(int ideg=2; ideg<3; ideg++) {  // linear-quadratic  [1]=linear [2]=quad
          det[ideg]      = _fe[ideg]->Jac(qp,xx_qnds,InvJac[ideg]);     // Jacobian
          JxW_g[ideg] =det[ideg]*_fe[ideg]->_weight1[ndim-1][qp];       // weight
          _fe[ideg]->get_phi_gl_g(ndim,qp,_phi_g[ideg]);               // shape funct
          _fe[ideg]->get_dphi_gl_g(ndim,qp,InvJac[ideg],_dphi_g[ideg]); // global coord deriv
        }

        // fields ------------------------------------------
        interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],
                      _phi_g[2],el_ndof[2],_ub_g[2]);                          // quadratic
#ifdef AXISYM
        JxW_g[2]  *=_ub_g[2][0];    
#endif

        // velocity ---------------------------------------------------
        _sP=1.e-20;          // stress tensor 
	mod2_vel =1.e-20;    // velocity field
#ifdef NS_EQUATIONS   
        // velocity field
        for(int idim=0; idim< DIMENSION; idim++) {
          vel_g[idim]=_ub_g[2][ns_idx+idim];
          mod2_vel += vel_g[idim]*vel_g[idim];
        }
        // derivatives
        interp_el_gdx(_data_eq[2].ub,ns_idx,DIMENSION,_dphi_g[2],el_ndof[2],vel_dxg);
#endif
        // stress tensor and velocity square modulus
        for(int idim=0; idim< DIMENSION; idim++) {
          for(int jdim=0; jdim< DIMENSION; jdim++) {  // turbulence source --
            const double sss=vel_dxg[idim+jdim*DIMENSION]+vel_dxg[jdim+idim*DIMENSION];
            _sP += sss*sss;
          } // upwind
          Pe_h[idim]=0.5*mod2_vel*h_eff[idim]/_IRe;
          f_upwind[idim]=UP_WIND_K*rho*0.5*(1./tanh(Pe_h[idim])-1./Pe_h[idim])*h_eff[idim]/sqrt(mod2_vel);
        }


        // eff visc turb at g pt
        _mu_turb=0.;  _kappa_g[0]= _ub_g[2][k_idx]; _kappa_g[1]=  _ub_g[2][k_idx+1];
        interp_el_gdx(_data_eq[2].ub,k_idx,2,_dphi_g[2],el_ndof[2],_ub_dxg);
        f_mu(val_tbg);  
        IRe_eff[0] = _IRe*(1.+_mu_turb);//  visc_eff  at g pt ;
//         IRe_eff[1] = _IRe*(1.+_mu_turb);//  visc_eff  at g pt ;
	
        /// d) Assemblying K-E equation
        for(int i=0; i< el_ndof[2]; i++) {  // +++++++++++
          // set up row
          for(int idim=0; idim<ndim; idim++) {
            dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
          }

          const double phii_g=_phi_g[2][i];
          // loop k,e
          for(int  ivar=0; ivar< _nvars[2]; ivar++)    {
            

            int  indx=i+ivar*el_ndof[2];//ivar=0;
            double dtxJxW_g=JxW_g[2]*_bc_vol[i+ivar*NDOF_FEM];

            // Assemblying rhs ----------------------------
            if(mode == 1)                {
              FeM(indx)  +=  dtxJxW_g*rho*(
                                _kappa_g[ivar+_dir]*phii_g/_dt    // _ub_g[2][k_idx+ivar+_dir]*phii_g/_dt   // time
                                +.5*_c_source[ivar+_dir]*_sP*phii_g  // source term
//                                 +_c_diss[ivar+_dir]*_ub_g[2][k_idx+ivar+_dir]*phii_g
                             );
            }

            // Assemblying Matrix quad ---------------------------------
            // ---------------------------------------------------------
            for(int j=0; j<el_ndof[2]; j++) {
              const double phij_g= _phi_g[2][j];
              // set up
              double Lap_g=0.,Adv_g=0.;

              for(int kdim=0; kdim<ndim; kdim++) {
                dphijdx_g[2][kdim] =_dphi_g[2][j+kdim*el_ndof[2]];
                Adv_g +=vel_g[kdim]*dphijdx_g[2][kdim];
                Lap_g += (IRe_eff[ivar+_dir]+f_upwind[kdim]*vel_g[kdim]*vel_g[kdim]
                        )*dphijdx_g[2][kdim]*dphiidx_g[2][kdim];
              }
              // diagonal blockce emus_op s [1-5-9]
              KeM(indx,j+ivar*el_ndof[2]) += dtxJxW_g*rho*(
                                               phij_g*phii_g/_dt                                         // time
                                               + Adv_g*phii_g  // advection
                                               + Lap_g     // viscous Laplacian
                                               + _c_diss[ivar+_dir]* phij_g*phii_g
                                             );
            } // end A element matrix quad -quad (end loop on j)----------------------
          } // end loop ivar
        } // end loop i

      } // end of the quadrature point qp-loop


      // ====================== end volume (element) =======================================

      // ======================================================================
      // =========================  boundary  =================================
      // ======================================================================
      for(int  iside=0; iside< el_sides; iside++) {
        if(el_neigh[iside] == -1) {

          // setup boundary element -> connectivity+coordinates
          for(int   lbnode=0; lbnode<NDOF_FEMB; lbnode++) {
            int lnode=_mgmesh._GeomEl._surf_top[lbnode+NDOF_FEMB*iside];// local nodes
            sur_toply[lbnode]=lnode;          // lbnode -> lnode
            elb_conn[lbnode]=el_conn[lnode];  // connctivity el_conn->elb_conn

            for(int   idim=0; idim<ndim; idim++) {  // coordinates
              double xyz=xx_qnds[idim*NDOF_FEM+lnode];
              _data_eq[2].ub[idim*NDOF_FEM+lbnode]=xxb_qnds[idim*NDOF_FEMB+lbnode]=xyz;
            }
          }
          // normal
          _fe[2]->normal_g(xxb_qnds,normal);
        
          // loop over kappa omega
          for(int ivar=0; ivar< _nvars[2]; ivar++)    {
//             IRe_eff = _IRe*(1.+_mu_turb*((ivar+_dir)* SIGMA_E+(1-ivar+_dir)*SIGMA_K));//  visc_eff  at g pt
            // Dirichlet boundary conditions  ***********************
            if(_bc_vol[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] ==0) {

              double det= _fe[2]->JacSur(elb_ngauss-1,xxb_qnds,InvJac[2]); // jacobian
              double Ipenalty=det*_dt;                           // Dirichlet bc flag

              for(int i=0; i< elb_ndof[2]; i++) {  // +++++++++++++++++++++++

                const int  indx_sol=sur_toply[i]+(ivar+_dir)*el_ndof[2];//indx sol;
                const int  indx_row=sur_toply[i]+ivar*el_ndof[2];       //indx row;
                double val=_data_eq[2].ub[k_idx*NDOF_FEM+indx_sol];     // old value

                // flag setup (\int_[S bc_normal*(u \dot n)+bc_tg*(u-(u.n)n)+bc_rhs*val] dS)
                int bc_s=(int)_bc_bd[indx_row];//  int bc_v=(int)_bc_vol[indx_sol];
                int bc_rhs=((bc_s&2)>>1);  // (1??) nonhomogeneous
//              int bc_normal=(bc_s%2);    // (??1) normal

                // Assemblying rhs // non-homogeneous ----------------------------
                if(mode == 1)  {
                  double kw_i=_data_eq[2].ub[k_idx*NDOF_FEM+indx_sol];
		  if(kw_i <1.e-20) kw_i =1.e-20;
                  double k_i=fabs(_data_eq[2].ub[k_idx*NDOF_FEM+sur_toply[i]]);
		  if(k_i <1.e-20) k_i =1.e-20;
		  double val=0.400772603053*k_i*sqrt(k_i)/(_y_dist);
// 		  double val=2.*_muf/(_rhof*BETASTAR*_y_dist*_y_dist);//sqrt(k_i)/(BETASTAR*_y_dist);
// 		  if(val >20000) val =20000;
// 		   val =1.e-0;
		  val *=CMU*k_i;
// 		  printf("val %lf k_i %lf \n",val,k_i);
                  FeM(indx_row) += Ipenalty*( // Dirichlet -> non-homogeneous flag (1??)
                                     bc_rhs*bc_rhs*kw_i// _u_old[indx_sol] (100) single comp bc
                                     +(1.-bc_rhs)*((1-ivar-_dir)*0.
                                                   +(ivar+_dir)*val
// 						   0.400772603053*k_i*sqrt(k_i)/(sqrt(_IRe)*_y_dist)
//            +(ivar+_dir)*2.*_muf/(BETASTAR*_y_dist*_y_dist)
                                                  )// w=1/delta// w=1/delta
                                   );
                }
                // Assemblying Matrix ---------------------------------
                KeM(indx_row,indx_row) += Ipenalty;  // (?1) Normal Dirichlet bc
              }// i loop
            } // end if  Dirichlet  ***********************************_c_diss[ivar+_dir]

            //  Neumann boundary conditions  ***********************************
            else if(_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] !=0) {

              const double alpha_t=1.e+12;//(1-ivar-_dir)/_y_dist-(ivar+_dir)/_y_dist;
              for(int qp=0; qp<  elb_ngauss; qp++) {  //gaussian integration loop (n_gauss)

                // quad/linear  [2]=quad [1]=linear------------------------------------
                det[2]     = _fe[2]->JacSur(qp,xxb_qnds,InvJac[2]);   // local coord _phi_g and jac
                JxW_g[2]=det[2]*_fe[2]->_weight1[ndim-2][qp]; // weight
                _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);   // global coord _phi_g
#ifdef AXISYM   // axisymmetric  (index ->0)
		interp_el_sol(_data_eq[2].ub,0,DIMENSION,_phi_g[2],elb_ndof[2],_ub_g[2]);
                JxW_g[2]  *=_ub_g[2][0];
#endif
                // ***********************************
                for(int i=0; i< elb_ndof[2]; i++) {  // Assemblying NS equation
                  // set up row i
                  const double phii_g=_phi_g[2][i];
                  const int indx_sol=sur_toply[i]+(ivar+_dir)*el_ndof[2];// volume dof index
                  const int indx_row=sur_toply[i]+ivar*el_ndof[2];// volume dof index
                  // boundary flag
                  int bc_s=(int)_bc_bd[indx_row];  int bc_v=(int)_bc_vol[indx_row];
                  double dtJxW_g=JxW_g[2]*((bc_v==0)?0:1);
                  int bc_rhs   =((bc_s&2)>>1); // (1??) -> nonhomogeneous
                  int bc_normal=(bc_s%2);      // (??1) -> normal

                  // Assemblying rhs ----------------------------
                  if(mode == 1)   {
// 		     double k_i=fabs(_data_eq[2].ub[k_idx*NDOF_FEM+sur_toply[i]]);_c_diss[ivar+_dir]
// 		  if(k_i <1.e-20) k_i =1.e-20;
//                       FeM(indx_sol)  += bc_rhs*alpha_t*dtJxW_g*_phi_g[2][i]*10000.; //0.400772603053/* *k_i*sqrt(k_i)*//(_y_dist);
                  }

                  // Assemblying Matrix ---------------------------------
                  for(int j=0; j<elb_ndof[2]; j++) {  // quad -quad
                    const double phij_g= _phi_g[2][j];
//                     KeM(indx_row,sur_toply[j]+(ivar+_dir)*el_ndof[2]) += dtJxW_g*alpha_t*phij_g*phii_g*bc_normal; // (??1) -> normal
                    KeM(indx_row,sur_toply[j]+ivar*el_ndof[2]) += 6.*IRe_eff[ivar+_dir]*dtJxW_g*phij_g*phii_g*bc_normal/_y_dist;// (??1) -> normal
                  }
                }// i-loop

              }// end gaussian  integration



            }  //  end if Neumann  ***********************************

          } // end if ivar   +++++++++++++++++++++++++++++++


        } //end if side
      } //end for iside

      // ==========================================================================
      // ====================== end boundary =====================================
      // ==========================================================================

      /// e) Add them to the global matrix ------------------------
      A[Level]->add_matrix(KeM,el_dof_indices);
      if(mode == 1) b[Level]->add_vector(FeM,el_dof_indices);
    } // end of element loop
    // clean ---------------------------------------------------
    el_dof_indices.clear();
     A[Level]->close(); if(mode == 1) b[Level]->close();
    // ----------------------------------------------------------
#ifdef PRINT_INFO
    std::cout << " GenMatRhs(NS): Assembled  " << Level << " (Level) "  << std::endl;
#endif

    return;
  }


// =========================================================================================
/// This function controls the assembly and the solution of the NS_equation system:
  void MGSolTBK::MGTimeStep(const double time, const int
                            /*iter*/) { // ------------------------------------

    std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;
    /// B) Assemblying of the rhs (top level in surface and volume with MGSolTBK::GenRhs,GenRhsB),
#if PRINT_TIME==1
    std::clock_t start_time=std::clock();
#endif
    
    GenMatRhs(time,_NoLevels-1,1);
    /// C) Assemblying of the  MGmatrices (MGSolTBK::GenMatrix),
    for(int Level = 0 ; Level < _NoLevels-1; Level++) GenMatRhs(time,Level,0); // matrix

#if PRINT_TIME==1
    std::clock_t end_time=std::clock();
    std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif

// //     / E) Solution of the linear MGsystem (MGSolTBK::MGSolve).
     MGSolve(1.e-6,40);

#if PRINT_TIME==11
    end_time=std::clock();
    std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif
    /// F) Update of the old solution at the top Level  (MGSolTBK::OldSol_update),
//     x_old[_NoLevels-1]->localize(*x_oold[_NoLevels-1]);
    x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);

    return;
  }

// =====================================================
/// This function computes the turbulent viscosity for
// =====================================================
  void MGSolTBK:: f_mu(double val1[]) {
    /// A) standard kappa-epsilon model
   double sP=sqrt(_sP);    double f_b=1.;        // Pv
   if (_kappa_g[0]< 1.e-20) _kappa_g[0]= 1.e-20; // kappa
   if (_kappa_g[1]< 1.e-20) _kappa_g[1]= 1.e-20; // epsilon

   // turbulent time=1/omega
  double tau_k=1.;  // turbulent time
  
  double l_mix=0.07*0.03025;
  
  
  #if (TBK_EQUATIONS==0)   //kappa  (Prandtl length)--------------------------
   tau_k=l_mix/sqrt(_kappa_g[0]);
  #endif  // end kappa -------------------------------------------------------
  #if (TBK_EQUATIONS/2==1)  // kappa-epsilon (TBK_EQUATIONS=2 or 3) ---------- 
  tau_k= CMU*_kappa_g[0]/_kappa_g[1]; // tau=1/omega=CMU*kappa/epsilon
  #endif   // kappa-epsilon --------------------------------------------------
  #if (TBK_EQUATIONS/2==2)       // kappa-omega ------------------------------
  tau_k=1/_kappa_g[1]; // tau=1/omega
  #endif
//   if (tau_k> MAX_TAU) tau_k= MAX_TAU;
  
  // turbulent viscosity=  k/tau
  double Re_t= _kappa_g[0]*tau_k/_IRe;       //  viscosity ratio vt/v
  if (Re_t > MU_TOP)  Re_t =MU_TOP;          //  limits
  if (Re_t < MU_LOW)  Re_t =MU_LOW;          //  limits  
  _mu_turb=Re_t;

  // Boundary corrections
  double R_t=Re_t/CMU;                                             // turbulent Reynolds number
  double y_star =_y_dist*sqrt(_kappa_g[0]/sqrt(R_t))/_IRe; 
  
    
  
 #ifdef LOWRE     /// B) Low-Reynolds  correction model
  double f_mu =(1.-exp(-1.*y_star/(14.)))*(1.-exp(-1.*y_star/14.));
  _mu_turb *=f_mu*(1.+5./pow(R_t,0.75)*exp(-1.*R_t*R_t/40000.));    // correction
 #endif
//     const double alpha_star=1.*(0.024+ 0.166666666667*Re_t)/(1.+0.166666666667*Re_t);
//     _mu_turb *=alpha_star;
//     // Turbulent viscosity
//     double sP=sqrt(_sP);    double f_b=1.;
//     
//     if(_kappa_g[0]< 1.e-20) _kappa_g[0]= 1.e-20;           // kappa
// //     if(_kappa_g[1]< 1.e-20) _kappa_g[1]= 1.e-10;           // epsilon
//     double omega=_kappa_g[1]/(CMU*_kappa_g[0]);
//     if(omega < 1.) omega=1.;
//     double Re_t=_kappa_g[0]/omega;// CMU*_kappa_g[0]*_kappa_g[0]/_kappa_g[1];  // turbulent viscosity   
  
    // coefficients
//    double f_e=1.;
//       double f_e=(1.-exp(-1.*y_star/(3.1)))*(1.-exp(-1.*y_star/(3.1)))*(1.-.3*exp(-1.*R_t*R_t/42.25));
    _c_source[0]=  _mu_turb*_IRe;//*exp(-_kappa_g[0]);
//     _c_source[1] =CMU*_c_source[0]/tau_k; //*exp(-_kappa_g[1]);  // standard kappa-omega model
    _c_diss[0]=   0.08*sqrt(_kappa_g[0])/l_mix;      //_kappa_g[1]/_kappa_g[0];//*exp(_kappa_g[0]);
//     _c_diss[1]=f_e*CMU*CMU/(tau_k*tau_k);//*exp(_kappa_g[1]);
    
  
    
 
    return;
  }

#endif // ----------  end TBK_EQUATIONS == 0 
#endif // ----------  end TBK_EQUATIONS  