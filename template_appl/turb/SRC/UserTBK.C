

#include "Equations_conf.h"  // <--- Equations configure
#ifdef TBK_EQUATIONS
// ============================================
#if (TBK_EQUATIONS/2==0) // 3D-2D kappa equation
// ============================================
#include <sstream>
#include "MGGeomEl.h"
// configuration files -----------
#include "Domain_conf.h"
#include "MGSTconf.h"
#include "MGFEconf.h"
#include "MGSTBKconf.h"
#include "MGSTTBKconf.h"   // turbulent energy
#include "Printinfo_conf.h"
#include "MGEquationsMap.h"

// class local configuration -------
#include "MGSolverTBK.h"

// local include -------------------
#include "MGMesh.h"
#include "MGSystem.h"
#include "numeric_vectorM.h"
#include "dense_vectorM.h"
#include "sparse_matrixM.h"
#include "dense_matrixM.h"
#include "linear_solverM.h"

#include "parallelM.h"
// ===============================
// Initial and boundary conditions
// ===============================



// ==========================================================
/// This routine constructs the FSI class:
MGSolTBK::MGSolTBK(
  MGEquationsMap& mg_equations_map_in,
  std::string eqname_in,
  std::string varname_in
):  MGSolDA(mg_equations_map_in,1,0,eqname_in,varname_in),
  /// A) reading class parameters
  // mesh params ------------
  _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes (top level)
  // phys params ------------
  _dt(_mgutils.get_par("dt")),       // parameter  dt
  _uref(_mgphys.get_par("Uref")),    // parameter  u reference
  _lref(_mgphys.get_par("Lref")),    // parameter  l reference
  _rhof(_mgphys.get_par("rho0")),    // parameter density
  _muf(_mgphys.get_par("mu0")) {     // parameter viscosity
  // ================================================

  /// B) setting class variables
#if (TBK_EQUATIONS==0)  
  _dir=0;  _var_names[0]="kt";    _refvalue[0]=_uref*_uref; // kappa
#endif
#if (TBK_EQUATIONS==1)  
  _dir=0;  _var_names[0]="vt";    _refvalue[0]=_uref*_lref; // kappa
#endif

  /// C ) setting solver type
  for(int l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVERTBK);
  /// D) setting nondimensional parameters
  _IRe=_muf/(_rhof*_lref*_uref);          // Reynolds number

  return;
}





// =========================================================
/// This function generates the initial conditions for the NS system:
void MGSolTBK::ic_read(
  double xp[],
  double u_value[]
) {// =======================================
//xp[] is the NON-DIMENSIONAL node coordinate
  double ILref = 1./_lref;

#if DIMENSION==2 // --------------  2D //--------------------------
  // xp[]=(xp,yp) u_value[]=(u,v,p)

#if (TBK_EQUATIONS%2==0)                // splitting 
  if(_dir==0)  { //  kappa
    u_value[0] = .1;
  }
  if(_dir==1)   { // omega
    u_value[0] = .1;
//     if(xp[1]<0.01) {
//        if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) { u_value[0] =4.; }
//      }
  }
#else                          // --------- coupled ---------------------------
  //  kappa -> 0  omega -> 1
  u_value[0] = 1.;
  u_value[1] = 1.;
//   if(xp[1]<0.01) {
//     u_value[1] =4.;
//     if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) {  u_value[0] =.1; }
//   }
//    if (xp[1]<0.001 && xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL)
//     u_value[1] =6852.79932546*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref;
//      if (xp[1]<0.001 && xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL)
//   if(xp[1]<0.01)
//     if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) u_value[1] =0.777+ 0*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref;
#endif


#endif  // //-----------------------//------------------------------------

// ===================================================================
#if DIMENSION==3  // --------------  3D //--------------------------

#if (TBK_EQUATIONS%2==0)    // ------- splitting ------------------- 
  if(_dir==0)  {//  kappa +++++++++++++++++++++++++
    u_value[0] = 0.;
  }
  if(_dir==1)   { // omega ++++++++++++++++++++++++
    u_value[0] = 0.;
    if(xp[1]<0.01) {
      if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) { u_value[0] =0.777+ 0*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref; }
    }
  }
#else               // --------- coupled ---------------------------
  // xp[]=(xp,yp,zp)   kappa -> 0  omega -> 1
  u_value[0]= 0.;
  u_value[1] =0.;
  if(xp[1]<0.01) {
    if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) { u_value[0] =0.777+ 0*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref; }
  }
#endif

#endif // //-------------------------------------------
  return;
}

// ========================================
/// This function  defines the boundary conditions for the NS system:
void MGSolTBK::bc_read(
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], // normal
  int bc_flag[]         // boundary condition flag
) {// ===================================
  //     0 ->  single component
  //     +4 ->  nonhomogeneous
  //     +2 -> tg
  //     +1 ->  normal
// ===================================
  double ILref = 1./_lref;

#if DIMENSION==2  // ----- 2D boundary conditions ----------
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann

#if (TBK_EQUATIONS%2==0)        // splitting 
  if(_dir==0) { //  kappa
    if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  // top 
    if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  //  left
    if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  //  right
    if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=2; bc_Neum[0]=0; }  // bottom
  }
  if(_dir==1) {  // omega
    if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  // top
    if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=1; bc_Neum[0]=1; }  //  left
    if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=1; bc_Neum[0]=1; }  //  right
    if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=2; bc_Neum[0]=0; }  // bottom
  }

#else           // --------- coupled ---------------------------
  //  kappa -> 0  omega -> 1
  if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=1; bc_Neum[1]=1;} // top
  if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=1; bc_Neum[1]=0;} //  left
  if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=1; bc_Neum[1]=0;} //  right
  if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=2; bc_flag[1]=2; bc_Neum[0]=0; bc_Neum[1]=0;} //  bottom
#endif // //---------------------------------------------

#endif // //---------------------------------------------

#if DIMENSION==3 // -------- 3D boundary conditions ---------
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
  //    boundary conditions box
#if (TBK_EQUATIONS%2==0)          // splitting 
  if(_dir==0) {  //  kappa
    if(xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  //  left
    if(xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  //  right
    if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  // top
    if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  // bottom
    if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  //  left
    if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  //  right
  }
  if(_dir==1) { // omega
    if(xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  //  left
    if(xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  //  right
    if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  // top
    if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=2; bc_Neum[0]=0; }  // bottom
    if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  //  left
    if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  //  right
#else            // --------- coupled ---------------------------
  if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=2; bc_Neum[0]=0; bc_Neum[1]=0;} //  INLET
  if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=0; bc_Neum[1]=1;} //  OUTLET
  if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=0; bc_Neum[1]=0;}
  if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=0; bc_Neum[1]=0;}
  if(xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=1; bc_Neum[1]=1;} //current
  if(xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=1; bc_Neum[1]=1;} //current
#endif
#endif  // //-----------------------//---------------------------

    return;
  }
  #if (TBK_EQUATIONS==0) // 3D-2D kappa equation
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
    const int k_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[K_F]];  // kappa-omega equations
    double    normal[DIMENSION]; double val_tbg[10];
    double    vel_g[DIMENSION]; double vel_dxg[DIMENSION*DIMENSION]; double mod2_vel=1.e-20;
    double Pe_h[DIMENSION]; double h_eff[DIMENSION]; double f_upwind[DIMENSION];
    for(int idim=0; idim< DIMENSION; idim++) {
      vel_g[idim]=(1-idim)*0.+idim*.7;
      mod2_vel += vel_g[idim]*vel_g[idim];
      for(int jdim=0; jdim< DIMENSION; jdim++) vel_dxg[idim+jdim*DIMENSION]=0.;
    }

    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) -----------
    A[Level]->zero();    if(mode ==1) b[Level]->zero();
    DenseMatrixM KeM;    DenseVectorM FeM;
    KeM.resize(el_mat_nrows,el_mat_ncols);    FeM.resize(el_mat_nrows);

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

      // local system fields
      get_el_dof_bc(Level,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

      // element nodes coordinates
      for(int idim=0; idim<ndim; idim++) {
        for(int d=0; d< NDOF_FEM; d++) {
          _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];
        } // diagonal element dimesion
        h_eff[idim]=fabs(xx_qnds[idim*NDOF_FEM+2+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+0]);
      }
      // external fields
      for(int eq=0; eq<_data_eq[2].n_eqs; eq++) {
        _data_eq[2].mg_eqs[eq]->
        get_el_sol(0,_data_eq[2].indx_ub[eq+1]-_data_eq[2].indx_ub[eq],
                   el_ndof[2],el_conn,offset,_data_eq[2].indx_ub[eq],
                   _data_eq[2].ub);
      }

      _y_dist=_mgmesh._dist[iel+nel_b];



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
        JxW_g[2]  *=_ub_g[0];      JxW_g[1]  *=_ub_g[0];
#endif
//  h_eff[0]=h_eff[1]=sqrt(JxW_g[2]);

        // velocity ---------------------------------------------------
        _sP=1.e-20;
#ifdef NS_EQUATIONS   
        int ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NSX_F]];
        // velocity field
        mod2_vel =1.e-20;
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



        // upwind
        /*    double  mod2_grad[2];    mod2_grad[0]=1.e-20;   mod2_grad[1] =1.e-20;
        double Pe_h[DIMENSION];    for(int idim=0; idim< DIMENSION; idim++) { Pe_h[idim]=0.5*fabs(vel_g[idim]+1.e-10)*h_eff[idim]/_IRe;
        fact[idim]=0.5*(1./tanh(Pe_h[idim])-1./Pe_h[idim])*h_eff[idim]/(fabs(vel_g[idim]+1.e-10));
        }*/
//    double Pe_h=0.5*mod2_vel*h_eff/_IRe;
//         const double f_upwind=UP_WIND_K*0.5*(1./tanh(Pe_h)-1./Pe_h)/sqrt(mod2_vel);

        // eff visc turb at g pt
        _mu_turb=0.;  _kappa_g[0]= _ub_g[2][k_idx]; 
        interp_el_gdx(_data_eq[2].ub,k_idx,1,_dphi_g[2],el_ndof[2],_ub_dxg);
//     for(int idim=0; idim< DIMENSION; idim++) {
//           mod2_grad[0] += _ub_dxg[idim]*_ub_dxg[idim];
//    mod2_grad[1] += _ub_dxg[idim+DIMENSION]*_ub_dxg[idim+DIMENSION];
//     }
//    mod2_grad[0] = sqrt( mod2_grad[0]);    mod2_grad[1] = sqrt( mod2_grad[1]);
//    double fupwind2[2];
//    Pe_h=0.5*mod2_grad[0]*h_eff/_IRe; fupwind2[0]=0.5*(1./tanh(Pe_h)-1./Pe_h)*UP_WIND_K/ mod2_grad[0];
//    Pe_h=0.5*mod2_grad[1]*h_eff/_IRe; fupwind2[1] =0.5*(1./tanh(Pe_h)-1./Pe_h)*UP_WIND_K/ mod2_grad[1];
//
//

        f_mu(val_tbg);  double IRe_eff = _IRe*(1.+_mu_turb);//  visc_eff  at g pt


        /// d) Assemblying K-E equation
        for(int i=0; i< el_ndof[2]; i++) {  // +++++++++++
          // set up row
          for(int idim=0; idim<ndim; idim++) {
            dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
          }

          const double phii_g=_phi_g[2][i];
          // loop u,v,w
          for(int  ivar=0; ivar< _nvars[2]; ivar++)    {

//    double Grad2_gi= 0;
//           for(int idim=0; idim<ndim; idim++) {
//      Grad2_gi +=IRe_eff*_ub_dxg[idim+(ivar+_dir)*DIMENSION]*_ub_dxg[idim+(ivar+_dir)*DIMENSION];
//           }
            int  indx=i+ivar*el_ndof[2];//ivar=0;
            double dtxJxW_g=JxW_g[2]*_bc_vol[i+ivar*NDOF_FEM];

            // Assemblying rhs ----------------------------
            if(mode == 1)                {
              FeM(indx)  +=  dtxJxW_g*rho*(
                                _kappa_g[0]*phii_g/_dt    // _ub_g[2][k_idx+ivar+_dir]*phii_g/_dt   // time
                                +.5*_c_source[0]*_sP*phii_g  // x-gravity
//                                 +_c_diss[ivar+_dir]*_ub_g[2][k_idx+ivar+_dir]*phii_g
                             );
            }

            // Assemblying Matrix quad ---------------------------------
            // ---------------------------------------------------------
            for(int j=0; j<el_ndof[2]; j++) {
              const double phij_g= _phi_g[2][j];
              // set up
              double Lap_g=0.,Adv_g=0.,Div_g= 0,Grad2_g= 0;
#ifdef  AXISYM
              Lap_g =(IRe_eff+fact[ivar+_dir]*UP_WIND_K*vel_g[ivar+_dir]*vel_g[ivar+_dir])*
                     (1-ivar)*phii_g*phij_g/(_ub_g[2][0]*_ub_g[2][0]);// axysimmetry
#endif
              for(int kdim=0; kdim<ndim; kdim++) {
                dphijdx_g[2][kdim] =_dphi_g[2][j+kdim*el_ndof[2]];
                Adv_g +=vel_g[kdim]*dphijdx_g[2][kdim];
                Lap_g += (IRe_eff+f_upwind[kdim]*vel_g[kdim]*vel_g[kdim]
                        )*dphijdx_g[2][kdim]*dphiidx_g[2][kdim];
              }
              // diagonal blockce emus_op s [1-5-9]
              KeM(indx,j+ivar*el_ndof[2]) += dtxJxW_g*rho*(
                                               phij_g*phii_g/_dt                                         // time
                                               + Adv_g*phii_g  // advection
                                               + Lap_g     // viscous Laplacian
                                               + _c_diss[0]* phij_g*phii_g
                                             );

              // ----------------------------------------------------------------------
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
              double Ipenalty=det/_dt;                           // Dirichlet bc flag

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
                  double k_i=fabs(_data_eq[2].ub[k_idx*NDOF_FEM+sur_toply[i]]);
                  FeM(indx_row) += Ipenalty*( // Dirichlet -> non-homogeneous flag (1??)
                                     bc_rhs*bc_rhs*kw_i// _u_old[indx_sol] (100) single comp bc
                                     +(1.-bc_rhs)*0.
                                   );
                }
                // Assemblying Matrix ---------------------------------
                KeM(indx_row,indx_row) += Ipenalty;  // (?1) Normal Dirichlet bc
              }// i loop
            } // end if  Dirichlet  ***********************************

            //  Neumann boundary conditions  ***********************************
            else if(_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] !=0) {

              const double alpha_t=(1-ivar-_dir)/_y_dist-(ivar+_dir)/_y_dist;
              for(int qp=0; qp<  elb_ngauss; qp++) {  //gaussian integration loop (n_gauss)

                // quad/linear  [2]=quad [1]=linear------------------------------------
                det[2]     = _fe[2]->JacSur(qp,xxb_qnds);   // local coord _phi_g and jac
                JxW_g[2]=det[2]*_fe[2]->_weight1[ndim-2][qp]; // weight
                _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);   // global coord _phi_g
#ifdef AXISYM   // axisymmetric  (index ->0)
                JxW_g[2]  *=_ub_g[0];
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
//                     FeM(indx_sol)  += -bc_rhs*dtJxW_g*_phi_g[0][i];
                  }

                  // Assemblying Matrix ---------------------------------
                  for(int j=0; j<elb_ndof[2]; j++) {  // quad -quad
                    const double phij_g= _phi_g[0][j];
                    KeM(indx_row,sur_toply[j]+ivar*el_ndof[2]) += dtJxW_g*alpha_t*phij_g*phii_g*bc_normal; // (??1) -> normal
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
    // ----------------------------------------------------------
#ifdef PRINT_INFO
    std::cout << " GenMatRhs(K): Assembled  " << Level << " (Level) "  << std::endl;
#endif

    return;
  }
#endif
#if (TBK_EQUATIONS==1) // 3D-2D kappa equation
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
    const int k_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[K_F]];  // kappa-omega equations
    double    normal[DIMENSION]; double val_tbg[10];
    double    vel_g[DIMENSION]; double vel_dxg[DIMENSION*DIMENSION]; double mod2_vel=1.e-20;
    double Pe_h[DIMENSION]; double h_eff[DIMENSION]; double f_upwind[DIMENSION];
    for(int idim=0; idim< DIMENSION; idim++) {
      vel_g[idim]=(1-idim)*0.+idim*.7;
      mod2_vel += vel_g[idim]*vel_g[idim];
      for(int jdim=0; jdim< DIMENSION; jdim++) vel_dxg[idim+jdim*DIMENSION]=0.;
    }

    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) -----------
    A[Level]->zero();    if(mode ==1) b[Level]->zero();
    DenseMatrixM KeM;    DenseVectorM FeM;
    KeM.resize(el_mat_nrows,el_mat_ncols);    FeM.resize(el_mat_nrows);

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

      // local system fields
      get_el_dof_bc(Level,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

      // element nodes coordinates
      for(int idim=0; idim<ndim; idim++) {
        for(int d=0; d< NDOF_FEM; d++) {
          _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];
        } // diagonal element dimesion
        h_eff[idim]=fabs(xx_qnds[idim*NDOF_FEM+2+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+0]);
      }
      // external fields
      for(int eq=0; eq<_data_eq[2].n_eqs; eq++) {
        _data_eq[2].mg_eqs[eq]->
        get_el_sol(0,_data_eq[2].indx_ub[eq+1]-_data_eq[2].indx_ub[eq],
                   el_ndof[2],el_conn,offset,_data_eq[2].indx_ub[eq],
                   _data_eq[2].ub);
      }

      _y_dist=_mgmesh._dist[iel+nel_b];



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
        JxW_g[2]  *=_ub_g[0];      JxW_g[1]  *=_ub_g[0];
#endif
//  h_eff[0]=h_eff[1]=sqrt(JxW_g[2]);

        // velocity ---------------------------------------------------
        _sP=1.e-20;
#ifdef NS_EQUATIONS   
        int ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NSX_F]];
        // velocity field
        mod2_vel =1.e-20;
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



        // upwind
        /*    double  mod2_grad[2];    mod2_grad[0]=1.e-20;   mod2_grad[1] =1.e-20;
        double Pe_h[DIMENSION];    for(int idim=0; idim< DIMENSION; idim++) { Pe_h[idim]=0.5*fabs(vel_g[idim]+1.e-10)*h_eff[idim]/_IRe;
        fact[idim]=0.5*(1./tanh(Pe_h[idim])-1./Pe_h[idim])*h_eff[idim]/(fabs(vel_g[idim]+1.e-10));
        }*/
//    double Pe_h=0.5*mod2_vel*h_eff/_IRe;
//         const double f_upwind=UP_WIND_K*0.5*(1./tanh(Pe_h)-1./Pe_h)/sqrt(mod2_vel);

        // eff visc turb at g pt
        _mu_turb=0.;  _kappa_g[0]= _ub_g[2][k_idx]; _kappa_g[1]=  _ub_g[2][k_idx+1];
        interp_el_gdx(_data_eq[2].ub,k_idx,2,_dphi_g[2],el_ndof[2],_ub_dxg);
//     for(int idim=0; idim< DIMENSION; idim++) {
//           mod2_grad[0] += _ub_dxg[idim]*_ub_dxg[idim];
//    mod2_grad[1] += _ub_dxg[idim+DIMENSION]*_ub_dxg[idim+DIMENSION];
//     }
//    mod2_grad[0] = sqrt( mod2_grad[0]);    mod2_grad[1] = sqrt( mod2_grad[1]);
//    double fupwind2[2];
//    Pe_h=0.5*mod2_grad[0]*h_eff/_IRe; fupwind2[0]=0.5*(1./tanh(Pe_h)-1./Pe_h)*UP_WIND_K/ mod2_grad[0];
//    Pe_h=0.5*mod2_grad[1]*h_eff/_IRe; fupwind2[1] =0.5*(1./tanh(Pe_h)-1./Pe_h)*UP_WIND_K/ mod2_grad[1];
//
//

        f_mu(val_tbg);  double IRe_eff = _IRe*(1.+_mu_turb);//  visc_eff  at g pt


        /// d) Assemblying K-E equation
        for(int i=0; i< el_ndof[2]; i++) {  // +++++++++++
          // set up row
          for(int idim=0; idim<ndim; idim++) {
            dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
          }

          const double phii_g=_phi_g[2][i];
          // loop u,v,w
          for(int  ivar=0; ivar< _nvars[2]; ivar++)    {

//    double Grad2_gi= 0;
//           for(int idim=0; idim<ndim; idim++) {
//      Grad2_gi +=IRe_eff*_ub_dxg[idim+(ivar+_dir)*DIMENSION]*_ub_dxg[idim+(ivar+_dir)*DIMENSION];
//           }
            int  indx=i+ivar*el_ndof[2];//ivar=0;
            double dtxJxW_g=JxW_g[2]*_bc_vol[i+ivar*NDOF_FEM];

            // Assemblying rhs ----------------------------
            if(mode == 1)                {
              FeM(indx)  +=  dtxJxW_g*rho*(
                                _kappa_g[ivar+_dir]*phii_g/_dt    // _ub_g[2][k_idx+ivar+_dir]*phii_g/_dt   // time
                                +.5*_c_source[ivar+_dir]*_sP*phii_g  // x-gravity
//                                 +_c_diss[ivar+_dir]*_ub_g[2][k_idx+ivar+_dir]*phii_g
                             );
            }

            // Assemblying Matrix quad ---------------------------------
            // ---------------------------------------------------------
            for(int j=0; j<el_ndof[2]; j++) {
              const double phij_g= _phi_g[2][j];
              // set up
              double Lap_g=0.,Adv_g=0.,Div_g= 0,Grad2_g= 0;
#ifdef  AXISYM
              Lap_g =(IRe_eff+fact[ivar+_dir]*UP_WIND_K*vel_g[ivar+_dir]*vel_g[ivar+_dir])*
                     (1-ivar)*phii_g*phij_g/(_ub_g[2][0]*_ub_g[2][0]);// axysimmetry
#endif
              for(int kdim=0; kdim<ndim; kdim++) {
                dphijdx_g[2][kdim] =_dphi_g[2][j+kdim*el_ndof[2]];
                Adv_g +=vel_g[kdim]*dphijdx_g[2][kdim];
//      Grad2_g += IRe_eff*_ub_dxg[kdim+(ivar+_dir)*DIMENSION]*dphijdx_g[2][kdim];
                Lap_g += (IRe_eff+f_upwind[kdim]*vel_g[kdim]*vel_g[kdim]
                        )*dphijdx_g[2][kdim]*dphiidx_g[2][kdim];
                Div_g +=_ub_dxg[kdim+kdim*ndim];
              }
              // diagonal blockce emus_op s [1-5-9]
              KeM(indx,j+ivar*el_ndof[2]) += dtxJxW_g*rho*(
                                               phij_g*phii_g/_dt                                         // time
                                               + Adv_g*phii_g  // advection
                                               + Lap_g     // viscous Laplacian
                                               + _c_diss[ivar+_dir]* phij_g*phii_g
//                                               + (IRe_eff+f_upwind*vel_g[ivar+_dir]*vel_g[ivar+_dir])
//                                               *(dphijdx_g[2][ivar+_dir]*dphiidx_g[2][ivar+_dir])        // viscous tensor
                                             );
//               int idimp1=(ivar+1+_dir)%ndim; // block +1 [2-6-7] ---------------
// #if (TBK_EQUATIONS%2==0)    // splitting
//               FeM(indx) += -_data_eq[2].ub[k_idx*NDOF_FEM+j+idimp1*el_ndof[2]]*//_u_old[j+idimp1*el_ndof[2]]*
// #else       // no splitting
//               KeM(indx,j+idimp1*el_ndof[2]) +=
// #endif
//                            dtxJxW_g* rho*(IRe_eff+f_upwind*vel_g[_dir]*vel_g[ivar+idimp1])*
//                            dphijdx_g[2][_dir]*dphiidx_g[2][idimp1];

              // ----------------------------------------------------------------------
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
              double Ipenalty=det/_dt;                           // Dirichlet bc flag

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
                  double k_i=fabs(_data_eq[2].ub[k_idx*NDOF_FEM+sur_toply[i]]);
                  FeM(indx_row) += Ipenalty*( // Dirichlet -> non-homogeneous flag (1??)
                                     bc_rhs*bc_rhs*kw_i// _u_old[indx_sol] (100) single comp bc
                                     +(1.-bc_rhs)*((1-ivar-_dir)*0.
                                                   +(ivar+_dir)*0.400772603053*k_i*sqrt(k_i)/(_y_dist)
// 						   0.400772603053*k_i*sqrt(k_i)/(sqrt(_IRe)*_y_dist)
//            +(ivar+_dir)*2.*_muf/(BETASTAR*_y_dist*_y_dist)
                                                  )// w=1/delta// w=1/delta
                                   );
                }
                // Assemblying Matrix ---------------------------------
                KeM(indx_row,indx_row) += Ipenalty;  // (?1) Normal Dirichlet bc
              }// i loop
            } // end if  Dirichlet  ***********************************

            //  Neumann boundary conditions  ***********************************
            else if(_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] !=0) {

              const double alpha_t=(1-ivar-_dir)/_y_dist-(ivar+_dir)/_y_dist;
              for(int qp=0; qp<  elb_ngauss; qp++) {  //gaussian integration loop (n_gauss)

                // quad/linear  [2]=quad [1]=linear------------------------------------
                det[2]     = _fe[2]->JacSur(qp,xxb_qnds);   // local coord _phi_g and jac
                JxW_g[2]=det[2]*_fe[2]->_weight1[ndim-2][qp]; // weight
                _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);   // global coord _phi_g
#ifdef AXISYM   // axisymmetric  (index ->0)
                JxW_g[2]  *=_ub_g[0];
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
//                     FeM(indx_sol)  += -bc_rhs*dtJxW_g*_phi_g[0][i];
                  }

                  // Assemblying Matrix ---------------------------------
                  for(int j=0; j<elb_ndof[2]; j++) {  // quad -quad
                    const double phij_g= _phi_g[0][j];
                    KeM(indx_row,sur_toply[j]+(ivar+_dir)*el_ndof[2]) += dtJxW_g*alpha_t*phij_g*phii_g*bc_normal; // (??1) -> normal
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
    // ----------------------------------------------------------
#ifdef PRINT_INFO
    std::cout << " GenMatRhs(NS): Assembled  " << Level << " (Level) "  << std::endl;
#endif

    return;
  }
#endif

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
    A[_NoLevels-1]->close();
    /// C) Assemblying of the  MGmatrices (MGSolTBK::GenMatrix),
    for(int
        Level = 0 ; Level < _NoLevels-1; Level++) {
      GenMatRhs(time,Level,0); // matrix
      A[Level]->close();
    }
#if PRINT_TIME==1
    std::clock_t end_time=std::clock();
    std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif

    /// E) Solution of the linear MGsystem (MGSolTBK::MGSolve).
    MGSolve(1.e-6,40);

#if PRINT_TIME==1
    end_time=std::clock();
    std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif
    /// F) Update of the old solution at the top Level  (MGSolTBK::OldSol_update),
    x_old[_NoLevels-1]->localize(*x_oold[_NoLevels-1]);
    x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);

    return;
  }

// =====================================================
/// This function computes the turbulent viscosity for
// =====================================================
  void MGSolTBK:: f_mu(double val1[]) {
    /// A) standard kappa-omega model
    // Turbulent viscosity
    double sP=sqrt(_sP);    double f_b=1.;
    if(_kappa_g[0]< 1.e-20) _kappa_g[0]= 1.e-10;   // kappa

#if (TBK_EQUATIONS==0)  // Mixing length ----------------------
    double Re_t= sqrt(_kappa_g[0])*_y_dist;    
    
    //  limits
    if( Re_t > MU_TOP*_IRe)  Re_t =MU_TOP*_IRe;
    if( Re_t < MU_LOW*_IRe)  Re_t =MU_LOW*_IRe;
     
    // coefficients
    _c_source[0]= Re_t;//*exp(-_kappa_g[0]);
    _c_diss[0]=CMU*sqrt(_kappa_g[0])/_y_dist;//*exp(_kappa_g[0]);
 #endif // ------------------------------------------------------
 
 #if (TBK_EQUATIONS==1)  // Spalart-Allarmas --------------------
    double Re_t=_kappa_g[0];    
    
    //  limits
    if( Re_t > MU_TOP*_IRe)  Re_t =MU_TOP*_IRe;
    if( Re_t < MU_LOW*_IRe)  Re_t =MU_LOW*_IRe;
     
    // coefficients
    _c_source[0]= Re_t;
    _c_diss[0]=CMU*sqrt(_kappa_g[0])/_y_dist;
 #endif // ------------------------------------------------------
 
    _mu_turb=Re_t/_IRe;
    
    
    
// #ifdef LOWRE     /// B) Low-Reynolds  correction model
//     double nu_eff=Re_t;
//     // diss coeff val[5]  (kappa eq) ---------------
//     _c_diss=BETA_W;//val[5]=
//     // eff visc val[6](k) val[7](w) (all eq) -------------------
//     _mu_turb=Re_t*0.5*(0.024+ 0.166666666667* Re_t)/(1.+0.166666666667* Re_t);//=val[6]val[6]
//     // source coeff val[4] (kappa eq) --------------
//     _c_source=ALFA_W*(0.111111111111+ 0.338983050847* Re_t)/(1.+0.338983050847* Re_t);//val[4]
//     // mixed
//     _c_F1=1.;//val[7]=1.;
// #endif
// #ifdef SST  /// C) SST model 
//     double nu_eff=Re_t;
// //      // F1 ----------------
//     double Phi1=500.*_muf/(_rhof*_y_dist*_y_dist*_omega_g);
//     double Phi2= Phi1;
//     double Phib=sqrt(_kappa_g)/(0.09*_omega_g*_y_dist);
//     if(Phi1 < Phib) Phi1=Phib;
//     double dplus=2*_rhof*_chi_k/1.168;  if(dplus < 1.e-10) dplus= 1.e-10;
//     double Phi1b=4.*_rhof*_kappa_g/(1.168*_y_dist*_y_dist*dplus);
//     if(Phi1 > Phi1b) Phi1 = Phi1b; Phi1 *=Phi1;
// 
//     // F1 and F2 coeff
//     double F1=tanh(Phi1*Phi1);       _c_F1= F1;                      //  F1val[7]
//     if(Phi2 <2.*Phib) Phi2=2.*Phib; const  double F2=tanh(Phi2*Phi2);// F2
// 
//     // eff visc coeff  fact_mu (all eq) ------------------------------------
//     nu_eff *=(0.024+ 0.166666666667* Re_t)/(1.+0.166666666667* Re_t);
//     double fact_mu=1;                            // turbulence factor  f_mu
//     double fact_mub=0.31*_omega_g/(sqrt(0.5*_sP)* F2+1.e-20);
//     if(fact_mu > fact_mub) fact_mu = fact_mub;
//     nu_eff *=fact_mu;                             // effective turbulent viscosity
//     _mu_turb=nu_eff*(0.5*F1+(1-F1)/1.168);        // turbulent viscosity =val[6]
//     // equation terms -------------------------------------------------
//     _c_source=(0.111111111111+ 0.338983050847* Re_t)/(1.+0.338983050847* Re_t)*
//               (F1*0.555555555556 +(1.-F1)*0.44);  // source coeff val[5]  (kappa eq)val[4]=
//     _c_diss=F1*0.075+(1.-F1)*0.0828;              // diss coeff val[6] (kappa eq)
// 
// #endif  // ====================================================  
// 
// 
// 
// 
// 
// 
// 
// 
// #ifdef LOWRE  /// B)  Low-Reynolds model    
//     // diss coeff val[5]  (kappa eq) ---------------
//     _chi_k *=1./(_kappa_g[1]*_kappa_g[1]);
//     if(_chi_k >0.) f_b=(1.+680.*_chi_k*_chi_k)*(1.+80.*_chi_k*_chi_k);
//     double  Re_t4=Re_t*0.125; Re_t4 *= Re_t4; Re_t4 *= Re_t4;
//     f_b=1.; _c_diss[0]=BETASTAR*f_b*(0.266666666667+Re_t4)/(1.+Re_t4);
// 
//     // eff visc val[6](k) val[7](w) (all eq) -------------------
//     double nu_t =Re_t*(0.024+ 0.166666666667* Re_t)/(1.+0.166666666667* Re_t);
//     _mu_turb=nu_t*0.5; // =val[6]
// 
//     // source coeff val[4] (kappa eq) --------------
//     _c_source[0]=nu_t*_IRe;// *(0.111111111111+ 0.338983050847* Re_t)/(1.+0.338983050847* Re_t);
// 
// #endif
// #ifdef SST  /// B)  SST-kappa-Omega model  
//     // diss coeff val[6] (kappa eq)
//     double  Re_t4=Re_t*0.125; Re_t4 *= Re_t4; Re_t4 *= Re_t4;
//     _c_diss[0]=BETASTAR*(0.266666666667+Re_t4)/(1.+Re_t4);
// 
//     // F1 ----------------
//     double Phi1=500.*_muf/(_rhof*_y_dist*_y_dist*_omega_g);
//     double Phi2=Phi1;
//     const double Phib=sqrt(_kappa_g)/(0.09*_omega_g*_y_dist);
//     if(Phi1 < Phib) Phi1=Phib;
//     double dplus=2*_rhof*_chi_k/1.168;  if(dplus < 1.e-10) dplus= 1.e-10;
//     const double Phi1b=4.*_rhof*_kappa_g/(1.168*_y_dist*_y_dist*dplus);
//     if(Phi1 > Phi1b)Phi1 = Phi1b; Phi1 *=Phi1;
//     const  double F1=tanh(Phi1*Phi1);
// 
//     if(Phi2 <2.*Phib) Phi2=2.*Phib; const  double F2=tanh(Phi2*Phi2);
//     // eff visc coeff  fact_mu (all eq
//     double nu_t  =Re_t*(0.024+ 0.166666666667* Re_t)/(1.+0.166666666667* Re_t);
//     double   fact_mu = 1.;  double fact_mub=0.31*_omega_g/(sqrt(0.5*sP)* F2+1.e-20);
//     if(fact_mub <1) fact_mu = fact_mub;
// //        if(fact_mu > 10*fact_mub) { fact_mu *=0.1;}
// // //        std::cout << fact_mu << "  bound corr " << fact_mub << "  " << std::endl;
// //       else      fact_mu = fact_mub;
// //     }
//     nu_t  *=fact_mu;
//     _mu_turb=nu_t*(F1/1.176+(1.-F1));// =val[6]
// 
//     // source coeff val[5]  (kappa eq)
//     _c_source[0]= nu_t*_IRe*(0.111111111111+ 0.338983050847* Re_t)/(1.+0.338983050847* Re_t);
// 
// #endif

    return;
  }









// end ---   f_mu -------------------------
#endif // ----------  end TBK_EQUATIONS  

 #endif

