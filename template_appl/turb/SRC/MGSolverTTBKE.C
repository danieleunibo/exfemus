
#include "Equations_conf.h"  // <--- Equations configure
 #ifdef TTBK_EQUATIONS

// ============================================
#if (TTBK_EQUATIONS/2==1) // 3D-2D kappa equation
// ============================================
// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverTTBK.h"       // Navier-Stokes class header file


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
MGSolTTBK::MGSolTTBK(
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
  _muf(mg_equations_map_in.get_par("mu0")),
  _Tref(mg_equations_map_in.get_par("Tref")),
  _cp0(mg_equations_map_in.get_par("cp0")),      // parameter  Cp reference
  _kappa0(mg_equations_map_in.get_par("kappa0")) // parameter  conductivity reference  
  {     // parameter viscosity
  // ================================================

  /// B) setting class variables
  _dir=0;
#if (TTBK_EQUATIONS%2==0)
  if(!varname_in.compare("kh")) _dir=0;  // kappa
  if(!varname_in.compare("eh")) _dir=1;  // epsilon
  _var_names[0]=varname_in;   _refvalue[0]=_Tref*_Tref;
#else
  _var_names[0]="kh";    _refvalue[0]=_Tref*_Tref; // kappa
  _var_names[1]="eh";    _refvalue[1]=_Tref*_Tref*_uref*_lref; // epsilon
#endif
  /// C ) setting solver type
  for(int l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVERTBK);
  /// D) setting nondimensional parameters
  _IRe=_muf/(_rhof*_lref*_uref);          // Reynolds number
  _alpha=_kappa0/(_rhof*_cp0);            // thermal diffusitvity
  _IPrdl=_rhof*_alpha/_muf;               // Prandtl number
  _IPrdl_turb=1./PRT;              // Turbulent Prandtl number
  _alpha_turb=0.;                         // Turbulent thermal diffusitvity
  
  return;
}


  
  



// ====================================================================
/// This function assembles the matrix and the rhs:
//  ======================================================================
  void  MGSolTTBK::GenMatRhs(const double time, const int
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
    const int kT_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[KTT_F]];  // energy turb equations   
    double    normal[DIMENSION]; double val_tbg[10];
    
    double T_dxg[DIMENSION];
    double Pe_h[DIMENSION]; double h_eff[DIMENSION]; double f_upwind[DIMENSION];
    double mod2_vel=1.e-20;  // velocity modulus
    double    vel_g[DIMENSION]; double vel_dxg[DIMENSION*DIMENSION];
    for(int idim=0; idim< DIMENSION; idim++) {
      vel_g[idim]=(1-idim)*0.+idim*.7;
      mod2_vel += vel_g[idim]*vel_g[idim];
      for(int jdim=0; jdim< DIMENSION; jdim++) vel_dxg[idim+jdim*DIMENSION]=0.;
    }

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
//      _dt=0.00001;

      // element fields ----------------------------------
      double phase=1.;    
double  alpha_eff = _IPrdl*_IRe;
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
        JxW_g[2]  *=_ub_g[2][0];   /*   JxW_g[1]  *=_ub_g[2][0];*/
#endif

      // Velocity field -> [NS_F] -> (quad, _indx_eqs[NS_F]) ------------------------------
      mod2_vel=1.e-20;  // velocity modulus
      _sP=1.e-20;       // turbulence production        
#ifdef NS_EQUATIONS   
      const int ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NS_F]];  // Navier-Stokes
   
      // velocity field derivatives (->vel_dxg[][])
      interp_el_gdx(_data_eq[2].ub,ns_idx,DIMENSION,_dphi_g[2],el_ndof[2],vel_dxg);
      
      // turbulence production term _sP
      for (int idim=0; idim< DIMENSION; idim++) {
	vel_g[idim] =_ub_g[2][ns_idx+idim]; // velocity field
        mod2_vel +=vel_g[idim]*vel_g[idim]; // velocity modulus
        for (int jdim=0; jdim< DIMENSION; jdim++) {
          _sP += (vel_dxg[jdim+idim*DIMENSION]+vel_dxg[jdim+idim*DIMENSION])*
                 (vel_dxg[idim+jdim*DIMENSION]+vel_dxg[jdim+idim*DIMENSION]);
        }
      }	
      mod2_vel =sqrt(mod2_vel); // velocity modulus
     
#ifdef TBK_EQUATIONS    // Turbulent viscosity  -> [K_F] -> (quad, _indx_eqs[K_F])  
      const int k_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[K_F]];  // kappa-equations
      _kappa_g[0]= _ub_g[2][k_idx];// val_tbg[0]=
      _kappa_g[1]= _ub_g[2][k_idx+1]; //  val_tbg[1]
#endif  // TBK_EQUATIONS        

#endif  // NS_EQUATIONS  ---------------------------------------------------------------

      // Temperature field -> [T_F] -> (quad, _indx_eqs[T_F]) --------------------------
      _sT=1.e-20;  // turbulence production  
#ifdef T_EQUATIONS 	
	 int nT_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F]];
	 interp_el_gdx(_data_eq[2].ub,nT_idx,1,_dphi_g[2],el_ndof[2],T_dxg);
        // stress tensor and velocity square modulus
	 T_dxg[1] += //36000*2.*0.03025/(10340.0*145.75*1.2*(0.03025*0.03025-0.004*0.004));
	               360./(10340.*145.75*V_MID*0.03025);
	 
        for(int idim=0; idim< DIMENSION; idim++) {
          _sT += T_dxg[idim]* T_dxg[idim];
          Pe_h[idim]=0.5*mod2_vel*h_eff[idim]/(_IRe*_IPrdl);
          f_upwind[idim]=UP_WIND_TTK*0.5*(1./tanh(Pe_h[idim])-1./Pe_h[idim])*h_eff[idim]/mod2_vel;
        }
#endif  // ------------------------------------------------------------------------------
        
	_kappaT_g[0]= _ub_g[2][kT_idx]; _kappaT_g[1]=  _ub_g[2][kT_idx+1];
        f_mu(val_tbg);  
	
// 	double IRe_eff = _IRe*(1.+_nut_ratio);//  visc_eff  at g pt
      


        /// d) Assemblying Kt-Et equation
        for(int i=0; i< el_ndof[2]; i++) {  // +++++++++++
          // set up row
          for(int idim=0; idim<ndim; idim++) {
            dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
          }

          const double phii_g=_phi_g[2][i];
          // loop u,v,w
          for(int  ivar=0; ivar< _nvars[2]; ivar++)    {
            alpha_eff = _IPrdl*_IRe+_nut_ratio*_IRe*_IPrdl_turb *SIGMA_KH;
	    if(ivar==0)   alpha_eff = _IPrdl*_IRe+_nut_ratio*_IRe*_IPrdl_turb*SIGMA_EH;
            int  indx=i+ivar*el_ndof[2];//ivar=0;
            double dtxJxW_g=JxW_g[2]*_bc_vol[i+ivar*NDOF_FEM];

            // Assemblying rhs ----------------------------
            if(mode == 1)                {
              FeM(indx)  +=  dtxJxW_g*(
                              1* _kappaT_g[ivar+_dir]*phii_g/_dt    // _ub_g[2][kT_idx+ivar+_dir]*phii_g/_dt   // time
                                  +_c_source[ivar+_dir]*_sT*phii_g
                                  + 0.5*(ivar+_dir)*_turb_source[ivar+_dir]*_sP*phii_g  // mechanical source for epsilon

//                                 +_c_diss[ivar+_dir]*_ub_g[2][kT_idx+ivar+_dir]*phii_g
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
                Lap_g += (alpha_eff+f_upwind[kdim]*vel_g[kdim]*vel_g[kdim]
                        )*dphijdx_g[2][kdim]*dphiidx_g[2][kdim];
              }
              // diagonal blockce emus_op s [1-5-9]
              KeM(indx,j+ivar*el_ndof[2]) += dtxJxW_g*(
                                             1*   phij_g*phii_g/_dt                                         // time
                                                + Adv_g*phii_g  // advection
                                                + Lap_g     // viscous Laplacian
                                                + _turb_source[0]*(ivar+_dir)*phij_g* phii_g 
//                                              + (IRe_eff+f_upwind*vel_g[ivar+_dir]*vel_g[ivar+_dir])
//                                               *(dphijdx_g[2][ivar+_dir]*dphiidx_g[2][ivar+_dir])        // viscous tensor
                                             );
               int idimp1=(ivar+1+_dir)%2; // block +1 [2-6-7] ---------------
 #if (TTBK_EQUATIONS%2==0)    // splitting
               FeM(indx) += -_data_eq[2].ub[kT_idx*NDOF_FEM+j+idimp1*el_ndof[2]]*//_u_old[j+idimp1*el_ndof[2]]*
 #else       // no splitting
               KeM(indx,j+idimp1*el_ndof[2]) +=
 #endif
//                             dtxJxW_g* //rho*(alpha_eff+f_upwind*vel_g[_dir]*vel_g[ivar+idimp1])*
//                             dphijdx_g[2][_dir]*dphiidx_g[2][idimp1];
                          _c_diss[ivar+_dir]*dtxJxW_g*phij_g*phii_g;			    

              // ----------------------------------------------------------------------(ret*0.001841/(10340*0.03025))
            } // end A element matrix quad -quad (end loop on j)----------------------
            // end B^T element matrix ------------------------------------
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
          
            // Dirichlet boundary conditions  ***********************
            if(_bc_vol[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] ==0) {

              double det= _fe[2]->JacSur(elb_ngauss-1,xxb_qnds); // jacobian
              double Ipenalty=det/(_dt);                           // Dirichlet bc flag
         
              for(int i=0; i< elb_ndof[2]; i++) {  // +++++++++++++++++++++++

                const int  indx_sol=sur_toply[i]+(ivar+_dir)*el_ndof[2];//indx sol;
                const int  indx_row=sur_toply[i]+ivar*el_ndof[2];       //indx row;
                double val=_data_eq[2].ub[kT_idx*NDOF_FEM+indx_sol];     // old value
                // flag setup (\int_[S bc_normal*(u \dot n)+bc_tg*(u-(u.n)n)+bc_rhs*val] dS)
                int bc_s=(int)_bc_bd[indx_row];//  int bc_v=(int)_bc_vol[indx_sol];
                int bc_rhs=((bc_s&2)>>1);  // (1??) nonhomogeneous
//              int bc_normal=(bc_s%2);    // (??1) normal

                // Assemblying rhs // non-homogeneous ----------------------------
                if(mode == 1)  {
                  double kw_i=_data_eq[2].ub[kT_idx*NDOF_FEM+indx_sol];
		  if(kw_i <1.e-20) kw_i =1.e-20;
		  double eps_i=fabs(_data_eq[2].ub[(kT_idx+1)*NDOF_FEM+sur_toply[i]]);
		  if(eps_i <1.e-20) eps_i =1.e-20;
                  double k_i=fabs(_data_eq[2].ub[kT_idx*NDOF_FEM+sur_toply[i]]);
		  if(k_i <1.e-20) k_i =1.e-20;
// 		    if(k_i >0.1) k_i =0.1;
// 		  double val=0;//.400772603053*k_i*sqrt(k_i)/(_y_dist);
		   val=2.*_muf/(_rhof*BETASTAR*_y_dist*_y_dist);//sqrt(k_i)/(BETASTAR*_y_dist)
                     if (val >20000)  val =20000.;
		  val *= BETASTAR*k_i*_IPrdl;

// 		printf(" I am here %g %g %g  \n  ", val,_IPrdl, k_i);
//            +(ivar+_dir)*2.*_muf/(BETASTAR*_y_dist*_y_dist)		  

// 		  double limit= eps_i/(_kappa_g[1]*_IPrdl);
// 		  if(k_i < limit) k_i = limit; 
		  
//  		  if(val >2000) val =2000;
                  FeM(indx_row) += Ipenalty*( // Dirichlet -> non-homogeneous flag (1??)
                                     bc_rhs*kw_i     // _u_old[indx_sol] (100) single comp bc
                                     +(1.-bc_rhs)*((1-ivar-_dir)*kw_i
                                                   +(ivar+_dir)*val
                                                  )// w=1/delta// w=1/delta
                                   );
                }
                // Assemblying Matrix ---------------------------------
                KeM(indx_row,indx_row) += Ipenalty;  // (?1) Normal Dirichlet bc
              }// i loop
            } // end if  Dirichlet  ***********************************

            //  Neumann boundary conditions  ***********************************
            else if(_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] !=0) {

              const double alpha_t=1.;//(1-ivar-_dir)/_y_dist-(ivar+_dir)/_y_dist;
              for(int qp=0; qp<  elb_ngauss; qp++) {  //gaussian integration loop (n_gauss)

                // quad/linear  [2]=quad [1]=linear------------------------------------
                det[2]     = _fe[2]->JacSur(qp,xxb_qnds);   // local coord _nut_ratio_phi_g and jac
                JxW_g[2]=det[2]*_fe[2]->_weight1[ndim-2][qp]; // weight
                _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);   // global coord _phi_g
		
#ifdef AXISYM   // axisymmetric  (index ->0)
	        interp_el_sol(_data_eq[2].ub,0,DIMENSION,_phi_g[2],elb_ndof[2],_ub_g[2]);
                JxW_g[2]  *=_ub_g[2][0];
#endif
                // ***********************************
                for(int i=0; i< elb_ndof[2]; i++) { 
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
// 		     double k_i=fabs(_data_eq[2].ub[kT_idx*NDOF_FEM+sur_toply[i]]);
// 		  if(k_i <1.e-20) k_i =1.e-20;
//                       FeM(indx_sol)  += bc_rhs*alpha_t*dtJxW_g*_phi_g[2][i]*10000.; //0.400772603053/* *k_i*sqrt(k_i)*//(_y_dist);
                  }

                  // Assemblying Matrix ---------------------------------
                  for(int j=0; j<elb_ndof[2]; j++) {  // quad -quad
                    const double phij_g= _phi_g[2][j];
                    KeM(indx_row,sur_toply[j]+(ivar+_dir)*el_ndof[2]) +=
                  (1-ivar+_dir)*  2.*_IPrdl*_IRe *dtJxW_g*phij_g*phii_g*bc_normal/_y_dist+
		  (ivar+_dir)*  2.*_IPrdl*_IRe *dtJxW_g*phij_g*phii_g*bc_normal/_y_dist;
                    
//                     dtJxW_g*alpha_t*phij_g*phii_g*bc_normal; // (??1) -> normal
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
    A[Level]->close();
    if (mode == 1) b[Level]->close();
    // ----------------------------------------------------------
#ifdef PRINT_INFO
    std::cout << " GenMatRhs(NS): Assembled  " << Level << " (Level) "  << std::endl;
#endif

    return;
  }


// =========================================================================================
/// This function controls the assembly and the solution of the NS_equation system:
  void MGSolTTBK::MGTimeStep(const double time, const int
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

    /// E) Solution of the linear MGsystem (MGSolTTBK::MGSolve).
       MGSolve(1.e-6,40);
//    std::cout << std::endl<< std::endl << " !!!!NOT SOLVING TK, PAY ATTENTION!!!!" << std::endl<< std::endl<< std::endl; 

#if PRINT_TIME==11
    end_time=std::clock();
    std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif
    /// F) Update of the old solution at the top Level  (MGSolTTBK::OldSol_update),
      x_old[_NoLevels-1]->localize(*x_oold[_NoLevels-1]);
      x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);
      
    return;
  }

// =====================================================
/// This function computes the turbulent viscosity for
// =====================================================
  void MGSolTTBK:: f_mu(double val1[]) {
    
 /// A) Momentum  Turbulent viscosity
_nut_ratio=0.;
#ifdef TBK_EQUATIONS
 if (_kappa_g[0]< 1.e-20) _kappa_g[0]= 1.e-20; // kappa
 if (_kappa_g[1]< 1.e-20) _kappa_g[1]= 1.e-20; // kappa
  double tau_k=1.;  // turbulent time

  #if (TBK_EQUATIONS==0)   //kappa  (Prandtl length)--------------------------
   tau_k=_y_dist/sqrt(_kappa_g[0]);
  #endif  // end kappa -------------------------------------------------------

  #if (TBK_EQUATIONS/2==1)  // kappa-epsilon (TBK_EQUATIONS=2 or 3) ---------- 
  tau_k= CMU*_kappa_g[0]/_kappa_g[1]; // tau=1/omega=CMU*kappa/epsilon
  #endif   // kappa-epsilon --------------------------------------------------

  #if (TBK_EQUATIONS/2==2)       // kappa-omega ------------------------------
  tau_k=1/_kappa_g[1]; // tau=1/omega
  #endif
  if (tau_k> MAX_TAU) tau_k= MAX_TAU;
  
  
  // turbulent viscosity
  double Re_t= _kappa_g[0]*tau_k/_IRe;
  if (Re_t > MU_TOP)  Re_t =MU_TOP;
  if (Re_t < MU_LOW)  Re_t =MU_LOW;
  _nut_ratio=Re_t;
  
      // Boundary corrections
   double R_t=Re_t/CMU;                                             // turbulent Reynolds number
   double R_eps = _y_dist*sqrt(_kappa_g[0]/sqrt(R_t))/_IRe;
   double y_plus = _y_dist*0.547722557505*sqrt(_kappa_g[0])/_IRe;
//    double R_eps =_y_dist/sqrt((_muf/_rhof)*sqrt((_muf/_rhof)/_kappa_g[1]));    // *sqrt(_kappa_g[0]/sqrt(R_t))/_IRe; 
   double f_mu =(1.-exp(-1.*R_eps/(14.)))*(1.-exp(-1.*R_eps/14.));
  _nut_ratio *=f_mu*(1.+5./pow(R_t,0.75)*exp(-1.*R_t*R_t/40000.));    // correction
  
//  // #ifdef LOWRE     /// B) Low-Reynolds  correction model
//     const double alpha_star=1.*(0.024+ 0.166666666667*Re_t)/(1.+0.166666666667*Re_t);
//     _nut_ratio *=alpha_star;
// // #endif  
  
  
 
  
#endif  // -----------------------------------------------------------------    
    /// A) Energy  Turbulent viscosity
    // Energy  Turbulent viscosity
	if(_kappaT_g[1]< 1.e-20) _kappaT_g[1]= 1.e-5;           // kappa
	if(_kappaT_g[0]< 1.e-20) _kappaT_g[0]= 1.e-5;           // kappa

    double tauT_k= CMU*_kappaT_g[0]/_kappaT_g[1];
    if (tauT_k> MAX_TAUT) tauT_k= MAX_TAUT; // tau=1/omega=CMU*kappa/epsilon
    double rT=tauT_k/tau_k; 
    
    
      double f_alpha =(1.-exp(-1.*R_eps/(19.*sqrt(_IPrdl)) ))*(1.-exp(-1.*R_eps/14.));
      double f_d=exp(-1.*R_t*R_t/(200.*200.));
      double f_asym =(1.-exp(-1.*R_eps/(19.*sqrt(_IPrdl)) ))*(1.-exp(-1.*R_eps/14.));

       double a_wall = sqrt(2.*rT)*1.3*_IPrdl/pow(R_t,0.75)*f_d*f_alpha;
       double a_inter =2.*rT/(.3 +rT)*f_alpha*exp(-1.*R_t*R_t/(500.*500.));
       double asymp = 0.9*f_asym;

      double b_nu = f_mu*(1.+5./pow(R_t,0.75)*exp(-1.*R_t*R_t/40000.));
      
      _IPrdl_turb=1.11111*(a_inter+a_wall+asymp)/b_nu;
//       if(_IPrdl_turb>1.) _IPrdl_turb=1.;
//       mu=10340*kt*kt/(et*0.001841)
//       reps=(0.03025-coordsX)*sqrt(kt/sqrt(mu/0.09))/(0.001841/10340)
//       f_alpha= (1-exp(-reps/(14/sqrt(0.023))))*(1-exp(-reps/14))
//       fd=exp(-(mu*mu/40000))
//       Prt=0.9/(falpha*(1/sqrt(R/0.023)+0.7*2*R/(0.5+R)+sqrt(2*R)*3*fd/(0.023*mu^0.75)))

         
//        
// 	 _IPrdl_turb=1./1.5;
//  _IPrdl_turb= 1./( 0.9+0.7*_IPrdl/_nut_ratio); 
     // coefficients

//       double cd2var =1.;//(1.9*(1.-0.3*exp(-R_t*R_t/(6.5*6.5)))-1.)*(1.-exp(-R_eps/5.7))*(1.-exp(-R_eps/5.7));//0.025*sqrt(R_t)*(1+3.*sqrt(1/(_IPrdl*rT))));
//       double cd1var =1.;//(1.-exp(-R_eps))*(1.-exp(-R_eps));//0.025*sqrt(R_t)*(1+3.*sqrt(1/(_IPrdl*rT))));
//       double cp1var = 1.;
//       double cp2var = 1.;
    _c_source[0]=_nut_ratio*_IRe*_IPrdl_turb;
    _turb_source[0] =CD2*CMU/tau_k;// mechanical dissipation for epsilon
    _c_source[1] =CP1*CMU*_c_source[0]/tauT_k;//*_kappaT_g[1]/_kappaT_g[0];
    _turb_source[1] = CP2*_kappaT_g[1]*tau_k;
    _c_diss[0]=1.;//_kappaT_g[1]/_kappaT_g[0];
    _c_diss[1]=CD1*CMU*CMU/(tauT_k*tauT_k);//_kappaT_g[1]/_kappaT_g[0];
//  		+ cd2var*CD2*CMU/tau_k; // mechanical dissipation for epsilon
                              
    

    return;
  }


// end ---   f_mu -------------------------
#endif

#endif // ----------  end TBK_EQUATIONS  
