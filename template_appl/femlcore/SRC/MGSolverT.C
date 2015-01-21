

#include "Equations_conf.h"
#ifdef  T_IS_PERSONAL

// ============================================
#ifdef T_EQUATIONS // 3D-2D Energy equation
// ============================================
#include <sstream>

// class local configuration -------
#include "MGSclass_conf.h"   // File conf class T
#include "MGSolverT.h"       // header file class T
#include "ReactData.h"       // Reactor data

#include "MGFE_conf.h" 
// #include "datagen.h"

#include "MGGeomEl.h"
#include "Printinfo_conf.h"
#include "MGEquationsMap.h"




// local include -------------------
#include "MGMesh.h"
#include "MGSystem.h"
#include "numeric_vectorM.h"
#include "dense_vectorM.h"
#include "sparse_matrixM.h"
#include "dense_matrixM.h"
#include "linear_solverM.h"

#include "parallelM.h"


// // ======================================================
// //  This function returns the planar fuel zone
// int zone(
//   double xm,
//   double ym
// ) {// ==================================================
//
//     int i_ass = (int)((ym+HR)/(2*HR));
//     int j_ass = (int)((xm+HR*((i_ass%2)))/(2*HR));
//     int zone_f= mat_zone[i_ass][j_ass];
//     if (zone_f == DUMMY) cerr<< "warn dummy zone \n" << std::endl;
//     return zone_f;
// }// ==================================================
//
// // ======================================================
// //    This function returns the planar fuel power peak
// double asby_power(
//   double xm,
//   double ym
// ) {// ==================================================
//     double peak_f;
//     int i_ass = (int)((ym+HR)/(2*HR));
//     int j_ass = (int)((xm+HR*(i_ass%2))/(2*HR));
//     peak_f =mat_pf[i_ass][j_ass];
//     return peak_f;
// }// ==================================================
//
// // ======================================================
//  //    This function returns the axial  fuel power peak
// // --------------------------------
// double axial_pf(
//   double x,
//   double y,
//   double z
// ) {// ==================================================
//
//     if (z >HIN && z< HOUT) {
//         int izone= zone(x,y);
//         int k=(int)(NZC*(z-HIN)/(HOUT-HIN));
//         if (izone > 2) return 1.;
//         else return axpf[k][izone];
//     }
//     else return 0.;
// }// ===================================================
//
// // ======================================================
// //    This function returns the presssure  loss
// //    over a plane
// double lossf(
//   double x,
//   double y,
//   double z
// ) {// ===================================================
//
//     // int izone= zone(x,y);
//     int k = 0;
//     if (z < HIN) k =(int)((2*z)/HIN);
//     else if (z < HOUT) k = 2+(int)(NZC*(z-HIN)/(HOUT-HIN));
//     else k = 12+(int)(2*(z-HOUT)/(2.-HOUT));
//
//     return axlpf[k][0];
// }// // ===================================================
//
// // ======================================================
// //    This function returns the  Axial  presssure  loss
// // ----------------------------------
// double axial_lf(
//   double z
// ) {// ===================================================
//     int k = 0;
//     if (z < HIN) k =(int)((2*z)/HIN);
//     else if (z < HOUT) k = 2+(int)(NZC*(z-HIN)/(HOUT-HIN));
//     else k = 12+(int)(2*(z-HOUT)/(2.-HOUT));
//
//     return axllf[k];
//
// } // ====================================================



// ======================================================
/// This function constructs the 3d-2D MGSolT class
MGSolT::MGSolT(MGEquationsMap& mg_equations_map_in,
               int nvars_in[],
               std::string eqname_in,
               std::string varname_in):
  MGSolDA(mg_equations_map_in,nvars_in,eqname_in,varname_in),
  /// A) reading parameters
  // mesh params ------------
  _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes
  // phys params ------------
  _dt(_mgutils.get_par("dt")),       // parameter  dt
  _uref(_mgphys.get_par("Uref")),    // parameter  vel reference
  _lref(_mgphys.get_par("Lref")),    // parameter  length reference
  _Tref(_mgphys.get_par("Tref")),    // parameter  temperature reference
  _rhof(_mgphys.get_par("rho0")),    // parameter  density reference
  _muf(_mgphys.get_par("mu0")),       // parameter  viscosity reference
  _cp0(_mgphys.get_par("cp0")),      // parameter  Cp reference
  _kappa0(_mgphys.get_par("kappa0"))// parameter  conductivity reference
//   _h_conv(_mgphys.get_par("hconv")),  // parameter  heat transfer convection coefficient
//   _T_inf(_mgphys.get_par("T_inf"))
{
  //  =================================================
  /// B) setting class ariables
  _var_names[0]=varname_in;
  _refvalue[0]=_Tref;

  /// C ) setting solver type
  for(int l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVERT);

  /// D) setting nondimensional parameters
  _alpha=_kappa0/(_rhof*_cp0);
  _IPrdl=_rhof*_alpha/_muf;
  _IRe=_muf/(_rhof*_uref*_lref);
  _IPrdl_turb=PRT;
  _alpha_turb=0.;
//   _Nusselt=(_h_conv*_lref)/_kappa0;
  _qheat=_mgphys.get_par("qheat")*_lref/(_rhof*_cp0*_Tref*_uref);
  _qs=_mgphys.get_par("qs")/(_rhof*_cp0*_Tref*_uref);

  return;
}



// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
// =================================================
// void MGSolT::ic_read(double xp[],double u_value[]) {
// // =================================================
// #if DIMENSION==2
// // xp[]=(xp,yp) u_value[]=(u,v,p)
// //   u_value[0] =573.15;
// //   if (xp[1]<0.1)
//   u_value[0] = 573.15;//*(1.-xp[1]);
// #else
// // =================================================
//   // xp[]=(xp,ypzp) u_value[]=(u,v,w,p)
//   u_value[0] = 0.*xp[0]* (1-xp[0]) *xp[1]* (2-xp[1]);
// //   if (xp[2]<0.1)
//   u_value[0] = 573.15;
// #endif
// }

// ========================================================
/// This function  defines the boundary conditions for the energy equation:
// void MGSolT::bc_read(double xp[],int bc_Neum[],int bc_flag[]) {
//   // =================================================
//   const double Lref = _mgphys.get_par("Lref");
//   double ILref = 1./Lref;
// #if DIMENSION==2
// // xp[]=(xp,yp) bc_flag[T]=0-Dirichlet 1-Neumann
//   //    boundary conditions box
//
//   if(xp[1]< LYB*ILref+BDRY_TOLL) {bc_Neum[0]=0; bc_flag[0]=2;}   // bottom
//   if(xp[0]< LXB*ILref+BDRY_TOLL) {bc_Neum[0]=1; bc_flag[0]=2;}   // left
//   if(xp[0]> LXE*ILref-BDRY_TOLL) {bc_Neum[0]=1; bc_flag[0]=0;}   // right
//   if(xp[1]> LYE*ILref-BDRY_TOLL) {bc_Neum[0]=1; bc_flag[0]=0;}   // top
// #endif
//
// #if DIMENSION==3
// // =================================================
//   // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
//   //    boundary conditions box
//
//   if(xp[1]> LYE*ILref-BDRY_TOLL)  {bc_flag[0]=0; bc_Neum[0]=1;}// top
//   if(xp[1]< LYB*ILref+BDRY_TOLL)  {bc_flag[0]=2; bc_Neum[0]=1;}// bottom
//   if(xp[2]< LZB*ILref+BDRY_TOLL)  {bc_flag[0]=0; bc_Neum[0]=1;}
//   if(xp[2]> LZE*ILref-BDRY_TOLL)  {bc_flag[0]=0; bc_Neum[0]=1;}
//   if(xp[0]< LXB*ILref+BDRY_TOLL)  {bc_flag[0]=0; bc_Neum[0]=1;}// left
//   if(xp[0]> LXE*ILref-BDRY_TOLL)  {bc_flag[0]=0; bc_Neum[0]=1;}// right
// #endif
//
//   return;
// } // end boundary conditions ==========================


//  ====================================================
/// This function assembles the matrix and the rhs:
//  ====================================================
void  MGSolT::GenMatRhs(
  const double /*time*/, // time  <-
  const int Level,  // Level <-
  const  int mode    // mode  <- (1=rhs) (0=only matrix)
) {  // ===============================================

  // -------------------------------------------------------------------------------------
  /// Set up
  // -------------------------------------------------------------------------------------

  // geometry ----------------------------------------------------------------------------
  const int  ndim = DIMENSION;                       //dimension
  const int  offset = _mgmesh._NoNodes[_NoLevels-1]; // mesh nodes
  const int  el_sides= _mgmesh._GeomEl._n_sides[0];  // element sides
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB]; // element connectivity
  int        el_neigh[NDOF_FEM];                     // bd element connectivity
  int        sur_toply[NDOF_FEMB];                   // boundary topology
  double     xx_qnds[DIMENSION*NDOF_FEM];            // element node coords
  double     xxb_qnds[DIMENSION*NDOF_FEMB];          // bd element node coords
  double     xm[DIMENSION];                          // mid element point

  // gauss integration  -----------------------------------------------------------------
  const int el_ngauss = _fe[2]->_NoGauss1[ndim-1];        // elem gauss points
  const int elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];  // bd elem gauss points
  double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];  // Jac, Jac*w Jacobean
  double dphijdx_g[3][DIMENSION];   // global shape-derivatives at g point
  double dphiidx_g[3][DIMENSION];   // global shape test-derivatives at g point

  // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------
  int el_ndof[3];    el_ndof[0]=1;     // number of local dofs on volume element
  int elb_ndof[3];  elb_ndof[0]=1;     // number of local dofs on bondarry element
  int el_mat_nrows =0;                 // number of global dofs (mat row)
  for(int ideg=2; ideg<3; ideg++) {
    el_ndof[ideg]=_fe[ideg]->_NoShape[ndim-1];
    elb_ndof[ideg]=_fe[ideg]->_NoShape[ ndim-2];
    el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
  };
  int el_mat_ncols = el_mat_nrows;                    // square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector

  // coupling  fields -------------------------------------------------------------------
  // Temperature
  int T_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F]];   // T index
  double Pe_h[3],f_upwind[3],h_eff[3]; // local Peclet,upwind,h_eff
  // Navier-Stokes
  const int ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NS_F]];// NS index
  double vel_g[DIMENSION];                                        // velocity
  for(int idim=0; idim<DIMENSION; idim++) vel_g[idim] =0.;        // velocity not coupled
  // Turbulence upwind
  double  val_tbg[3];

  // From primary
  double T_in=608.15; double p_in=0.; double u_in=0.;double v_in=0.;double w_in=0.;
// #ifdef COUPLED_MESH
//   FILE *fp; fp = fopen("cmsh2.txt", "r");       // open the file
//   fscanf(fp,"%lf,%lf,%lf", &T_in, &v_in,&p_in); // read from file
//   fclose(fp);                                   // close the file
// #endif
 #ifdef COUPLED_MESH
  u_in=_mgmesh._mesh2_data[5+0]; v_in=_mgmesh._mesh2_data[5+1];
  w_in=_mgmesh._mesh2_data[5+2]; p_in=_mgmesh._mesh2_data[5+DIMENSION];
  T_in=_mgmesh._mesh2_data[5+DIMENSION+1];
  #endif
//   vel_g[2]=-sqrt(u2_in*u2_in+v2_in*v2_in+w2_in*w2_in); //set velocity along z-axis


  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) ----------------
  A[Level]->zero();    if(mode ==1) b[Level]->zero();   // global matrix+rhs
  DenseMatrixM KeM;    DenseVectorM FeM;                // local  matrix+rhs


  KeM.resize(el_mat_nrows,el_mat_ncols);
  FeM.resize(el_mat_nrows);        // resize  local  matrix+rhs

  // ------------------------------------------------------------------------------------
  /// b) Element  Loop over the volume (n_elem)
  // ------------------------------------------------------------------------------------

  const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
  const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
  for(int iel=0; iel < (nel_e - nel_b); iel++) {

    // set to zero matrix and rhs and center
    KeM.zero();        FeM.zero();

    // geometry and element  fields -------------------------------------------------------
    // Element Connectivity (el_conn)  and coordinates (xx_qnds)
    _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);
    _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);

    // set element-nodes variables  bc (bc_q_dofs)
    get_el_dof_bc(Level,iel,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

    // grid and gaussian points -----------------------------------------------------------
    for(int idim=0; idim<DIMENSION; idim++) {
      double sum=0.;                                             // element middle point
      for(int d=0; d< NDOF_FEM; d++) {
        _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];// element nodes xxg (DIM)
        sum += xx_qnds[idim*NDOF_FEM+d];
      }
      xm[idim]=sum/NDOF_FEM;                                     // element middle point
      // element length (h_eff[idim])
      h_eff[idim]=fabs(xx_qnds[idim*NDOF_FEM+2+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+0]);
      double h_1=fabs(xx_qnds[idim*NDOF_FEM+3+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+1]);
      if(h_eff[idim] <  h_1) h_eff[idim]=h_1;  // Max dx diagonal term
      f_upwind[idim]=0.;                       // upwind

    }
    // element field values ---------------------------------------------------------------
    for(int deg=0; deg<3; deg++) {  // loop over polinomial degree
      for(int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
        _data_eq[deg].mg_eqs[eq]->get_el_sol(
          0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],
          el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub);
      }
    }

    //  external node (gaussian point ) quantities ----------------------------------------
    // Non-dimensional Physical properties
    double phase=1.; double rhocp=1.;   double    rho =1.; //
    // Heat source
    double Qad =_qheat*_mgphys._sys_data[0][iel+nel_b];    // volume  heat source
    double qn =0;                                         // surface  heat source

// Old stuff !!!!!!!!!!!!!!!
//     double peak_f = asby_power(xm[0],xm[1])*axial_pf(xm[0],xm[1],xm[2]);
//     double loss_f =lossf(xm[0],xm[1],xm[2])/axial_lf(xm[2]);
//
//             if (xm[2]>HIN &&  xm[2]<HOUT) {
// //                 std::cout << iel+nel_b << " " << xm[0] << " " << xm[1] << " " << xm[2] << " " << peak_f  << " " << loss_f << endl;
// //                 aux_data_pf[Level2][iel]=peak_f ;
// //                 aux_data_pl[Level2][iel]=loss_f;
// 		 _mgphys._sys_data[0][iel+nel_b]=peak_f;
// 		  _mgphys._sys_data[1][iel+nel_b]=loss_f;
//             }
//             else {
//
// 		_mgphys._sys_data[0][iel+nel_b]=peak_f=0.;
// 		_mgphys._sys_data[1][iel+nel_b]=loss_f=0.;
//             }
// if (_ub_g[2][2]> 0.5 &&_ub_g[2][2]< 1.5) Qad=_qheat;//*exp(-1*time*30);

    //  external cell (element ) properties -------------------------------------
#ifdef TBK_EQUATIONS     // wall distance (for turbulence  field)  
    _y_dist=_mgmesh._dist[ iel+nel_b];
#endif // ------------------------------------------------------------  





    // ======================================================================
    // Volume =============================================================
    // ======================================================================

    // ---------------------------------------------
    /// c) gaussian integration loop (n_gauss)
    // --------------------------------------------
    for(int qp=0; qp< el_ngauss; qp++) {

      // shape functions at gaussian points -----------------------------------------
      for(int ideg=2; ideg<3; ideg++) {  // linear-quadratic  [1]=linear [2]=quad
        det[ideg]      = _fe[ideg]->Jac(qp,xx_qnds,InvJac[ideg]);     // Jacobian
        JxW_g[ideg] =det[ideg]*_fe[ideg]->_weight1[ndim-1][qp];       // weight
        _fe[ideg]->get_phi_gl_g(ndim,qp,_phi_g[ideg]);               // shape funct
        _fe[ideg]->get_dphi_gl_g(ndim,qp,InvJac[ideg],_dphi_g[ideg]);// global deriv
      }

      //  fields -----------------------------------------------------------
      interp_el_sol(_data_eq[2].ub,0,_data_eq[2].
                    indx_ub[_data_eq[2].n_eqs],_phi_g[2],el_ndof[2],_ub_g[2]); // quadratic

      double alpha_eff =_IPrdl*_IRe;
      double dens=_rhof;
#ifdef AXISYM   // axisymmetric (index -> 0)
      JxW_g[2]  *=_ub_g[2][0];
#endif

      // Temperature field  (T_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F]])
      double Temp_g=_ub_g[2][T_idx]; // local temp
#ifdef TEMPERATURE_DEPENDENCE_T
      dens = densityT(Temp_g);       // density (densityNS(x) in MGSclass_conf.h)
      rho =  dens/_rhof;              // normlized  local density
      rhocp = rho;               // normlized local viscosity
      if (Temp_g > T_FREEZE_SOLID && Temp_g < T_FREEZE_SOLID+DELTA_T_FREEZE) rhocp=100.;
	
      // printf(" Temp %f; rho %f %f ;  visc %f %f",Temp_g, rho, dens, mu, visc);
      //       alpha_eff *=kappaT(uold_g[0]);// nondimensional viscosity (mu/muref)
#endif

      // Velocity field (ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NS_F]])
      double mod2_vel=1.e-20;    // velocity modulus
      _sP=1.e-20;                // turbulence production
#ifdef NS_EQUATIONS
      double vel_dxg[DIMENSION*DIMENSION];
      interp_el_gdx(_data_eq[2].ub,ns_idx,DIMENSION,_dphi_g[2],el_ndof[2],vel_dxg);
      // turbulence production term _sP
      for(int idim=0; idim< DIMENSION; idim++) {
        vel_g[idim] =_ub_g[2][ns_idx+idim];// velocity vector (vel_g[])
        mod2_vel +=vel_g[idim]*vel_g[idim];// velocity modulus (mod2_vel)
        for(int jdim=0; jdim< DIMENSION; jdim++) {
          _sP += (vel_dxg[jdim+idim*DIMENSION]+vel_dxg[jdim+idim*DIMENSION])*
                 (vel_dxg[idim+jdim*DIMENSION]+vel_dxg[jdim+idim*DIMENSION]);
        }
        // upwind
        Pe_h[idim]=0.5*mod2_vel*h_eff[idim]/_IRe;     //  local Peclet number
        f_upwind[idim]=UP_WIND_T*rhocp*0.5*(1./tanh(Pe_h[idim])-1./Pe_h[idim])*
                       h_eff[idim]/sqrt(mod2_vel);    // f_upwind[] upwind
      }
      mod2_vel =sqrt(mod2_vel);  // velocity modulus
#ifdef TBK_EQUATIONS    // Turbulent viscosity  (k_idx=quad.indx_ub[_data_eq[2].tab_eqs[K_F]])  
      const int k_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[K_F]];  // kappa-equations
      _kappa_g[0]= _ub_g[2][k_idx];// val_tbg[0]=
      _kappa_g[1]= _ub_g[2][k_idx+1]; //  val_tbg[1]
      f_mu(val_tbg);
      alpha_eff = _IPrdl*_IRe*(1.+ _mu_turb*_IPrdl_turb/_IPrdl);
#endif  // TBK_EQUATIONS  
#endif  // NS_EQUATIONS       

      // ----------------------------------------------------------------------------
      /// d) Local (element) assemblying energy equation
      // ----------------------------------------------------------------------------

      for(int i=0; i<el_ndof[2]; i++)     {
        // set up row i
        const double phii_g=_phi_g[2][i];
        for(int idim=0; idim< ndim; idim++)  dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
        double dtxJxW_g=JxW_g[2]*_bc_vol[i];

        // Rhs Assemblying  ----------------------------------------------------------------------
        if(mode == 1) {         // rhs
          FeM(i) += dtxJxW_g*(
                      rhocp*Temp_g*phii_g/_dt   // term        (T_t)
                      +Qad*phii_g               // heat source (Q )
// 	              +Qad*(0.2146*pow(_dt,-0.3239))*phii_g    // decay heat source (Q )
                    );
        } // --------------------------------------------------------------------------------------
        // Matrix Assemblying ---------------------------------------------------------------------
        for(int j=0; j<el_ndof[2]; j++) {
          double phij_g= _phi_g[2][j];                            // shape function
          double Adv=0.;          double Lap=0.;
          for(int idim=0; idim< ndim; idim++) {
            dphijdx_g[2][idim]=_dphi_g[2][j+idim*el_ndof[2]];
            Adv +=vel_g[idim]*dphijdx_g[2][idim];                     // advection
            Lap +=(alpha_eff+f_upwind[idim]*vel_g[idim]*vel_g[idim])* //(alpha_eff+f*v^2)*
                  dphijdx_g[2][idim]*dphiidx_g[2][idim];              // Laplacian
          }
//                     printf(" mm %lf %lf  %lf \n", alpha_eff, f_upwind[1]*vel_g[1]*vel_g[1], f_upwind[0]*vel_g[0]*vel_g[0]);
          // energy-equation
          KeM(i,j) +=dtxJxW_g*(
                       rhocp*phii_g*phij_g/_dt // time term        (T_t)
                       + rhocp*Adv*phii_g      // advection term   (u . DT)
                       + Lap                   // diffusion term   (D[aDT])
                     );
        }
      } // ----------------------------------------------------------------------------------------
    } // end of the quadrature point qp-loop ***********************


    // ======================================================================
    // ====================== boundary ======================================
    // ======================================================================

    for(int iside=0; iside< el_sides; iside++)  {
      if(el_neigh[iside] == -1) {

        double alpha_eff = _IPrdl*_IRe;

        for(int idof=0; idof<NDOF_FEMB; idof++) {
          sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
          int idofb=sur_toply[idof];
          elb_conn[idof]=el_conn[idofb];     // connectivity vector
          for(int idim=0; idim<DIMENSION; idim++) {
            xxb_qnds[idim*NDOF_FEMB+idof]=xx_qnds[idim*NDOF_FEM+idofb];   // coordinates
            _data_eq[2].ub[idim*NDOF_FEM+idof]=xxb_qnds[idim*NDOF_FEMB+idof];
          }
        }
        // ======================================================================
        // Dirichlet boundary conditions
        if(_bc_vol[sur_toply[NDOF_FEMB-1]] ==0) {

          //  [NDOF_FEMB-1] is the midpoint of a quadaratic FEM element (HEX9 (2D) or HEX27 (3D))
          int bc_s=(int)_bc_bd[sur_toply[NDOF_FEMB-1]];     // b cond
          double det= _fe[2]->JacSur(elb_ngauss-1,xxb_qnds);// jacobian
          double Ipenalty=det/_dt;                          // Dirichlet bc flag
          // local boundary loop   ---------------------------------------
          for(int lb_node=0; lb_node< elb_ndof[2]; lb_node++) {
            int lv_node= sur_toply[lb_node]; // local vol index

            // flag setup: ---------------------------------------------------------------------
            // surface flag ->  \int (bc_var(?1)*T+bc_val(1?)*val)ds
            int  bc_val = (int)((bc_s&2)>>1);  // (1?)non homogeneous
            int  bc_var = (int)(bc_s%2);    // (?1) variable

            // Assemblying  Matrix & rhs --------------------------------------------------------
            if(mode == 1) {
              double value=0;

//               value= ((xxb_qnds[2*NDOF_FEMB+8]<0.01) ? T_in:_data_eq[2].ub[T_idx*NDOF_FEM+lv_node]);//get the values only at the inlet
              
	      value= ((xxb_qnds[2*NDOF_FEMB+sur_toply[8]]<1.3) ? 
// 	      ((xm[0]<xm[1])? 598.7:673.15):
	      ((xm[0]<xm[1])? 599.498:608.15):
	      _data_eq[2].ub[T_idx*NDOF_FEM+lv_node]);//get the values only at the inlet
value=_data_eq[2].ub[T_idx*NDOF_FEM+lv_node];//standard told

              FeM(lv_node) += bc_val*Ipenalty*value;
            }
            KeM(lv_node,lv_node) += Ipenalty;  //  Dirichlet bc
            // ---------------------------------------------------------------------------------
          }// lb_node -end  local boundary loop -------------------------
        } // end if Dirichlet  boundary conditions

        // ======================================================================
        // Neumann boundary conditions
        else if(_bc_bd[sur_toply[NDOF_FEMB-1]] !=0) {


          // gaussian integration loop (n_gauss)
          // -----------------------------------------------
          for(int qp=0; qp<  elb_ngauss; qp++) {

            // quad/linear  [2]=quad [1]=linear------------------------------------
            det[2]     = _fe[2]->JacSur(qp,xxb_qnds);    // local coord _phi_g and jac
            JxW_g[2]=det[2]*_fe[2]->_weight1[ndim-2][qp];// weight
            _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);   // global coord _phi_g
#ifdef AXISYM   // axisymmetric  
            JxW_g  *=_ub_g[0];
#endif
            // ***********************************
            // local side loop (over the node face)
            for(int lsi_node=0; lsi_node< elb_ndof[2]; lsi_node++) {

              // set up row i
              const double phii_g=_phi_g[2][lsi_node]; // boundary test function
              int lei_node= sur_toply[lsi_node]; // local element index

              // flag setup: ---------------------------------------------------------------------
              // volume flag ->   \int_vol T_eqs dv *bc_v
              int bc_v=(int)_bc_vol[lei_node];
              double dtxJxW_g=JxW_g[2]*bc_v;      // Neumann bc flag and Jac
              // surface flag ->  \int (bc_var(?1)*T+bc_val(1?)*val)ds
              int bc_s=(int)_bc_bd[lei_node];
              int  bc_val = (int)((bc_s&2)>>1);  // (1?) non-homogeneous
              int  bc_var = (int)(bc_s%2);       // (?1) tau=A*variable
              // ---------------------------------------------------------------------------------

              // Assemblying rhs ----------------------------------------------------------------
//               if (mode == 1){

// 		if(xxb_qnds[2*NDOF_FEM+NDOF_FEM-1]<0.001) qn=1.;
//                 FeM(lei_node) += bc_val*dtxJxW_g*phii_g*(// non-homogeneous
//                                      0.
//                                    bc_var*qn // flusso
//                                    +(1 -bc_var)*qn                   // Neumann heat flux bc
//                                  );
// 	      } // -------------------------------------------------------------------------------
              // Assemblying Matrix --------------------------------------------------------------
//               for (int lsj_node=0; lsj_node< elb_ndof[2];  lsj_node++) {
//                 KeM(lei_node,sur_toply[lsj_node]) += dtxJxW_g*bc_var*alpha_eff*phii_g*_phi_g[2][lsj_node]; // Robin bc  (k*dt/dn = h*(-T))
//               }// end j  -------------------------------------------------------------------------

            }// i   +++++++++++++++++++++++++++++++
          } // end of the quadrature point qp-loop **********************


        } // Neumann non homog

      } //end if side
    } // ======================  end for boundary =======================================

    /// e) Global assemblying energy equation
    A[Level]->add_matrix(KeM,el_dof_indices);                  // global matrix
    if(mode == 1)   b[Level]->add_vector(FeM,el_dof_indices);  // global rhs

  } // end of element loop
  // clean
  el_dof_indices.clear();
  A[Level]->close();
  if(mode == 1)  b[Level]->close();
#ifdef PRINT_INFO
  std::cout<< " Matrix Assembled(T)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif
//   A[Level]->close();
//   double nnnorm=A[Level]->linfty_norm();
//  std::cout<< "norm " << nnnorm << "\n";


  return;
}

// =========================================================================================
/// This function compute average surface  values on  over mesh[0] and
///  transfer data to mesh[1]
void MGSolT::print_ext_data(
  double mesh2mesh_data[]
) {
// ========================================================================================

  int flag_probe=_mgphys._flag_probe;
  const int Level=_NoLevels-1;
  const int offset= _mgmesh._NoNodes[_NoLevels-1];
  // Compute temperature on surface defined in _mgphys._node_glob_probe
  double sum_u=0.;
  for(int inode=0; inode< flag_probe; inode++) {
    sum_u += (*x_old[Level])(_node_dof[Level][_mgphys._node_glob_probe[inode]])*_refvalue[0];
  }
  sum_u /=flag_probe;
  // Writing average temperature from surface over mesh[0] to mesh[1]
  _mgmesh._mesh2_data[DIMENSION+1]=sum_u;  // writing on mesh[0]
  mesh2mesh_data[DIMENSION+1]=sum_u;       // writing on mesh[1]


#ifdef PRINT_INFO
  std::cout << " probe T computed  " << Level << " (Level) with "  << flag_probe  << " (nodes) "  << std::endl;
  
  std::cout << " 0n T probe_count " << flag_probe << std::endl;
  std::cout << " Temperature = " <<  mesh2mesh_data[DIMENSION+1]  << std::endl;
  
#endif
  return;
}


#include <iomanip>
// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolT::MGTimeStep(
  const double time,  // time
  const int /*iter*/  // Number of max inter
) {
// =========================================================================================

/// A) Set up the time step
  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;

  /// B) Assemblying of the Matrix-Rhs
#if PRINT_TIME==1
  std::clock_t start_time=std::clock();
#endif
  GenMatRhs(time,_NoLevels-1,1);                                                // matrix and rhs
  for(int Level = 0 ; Level < _NoLevels-1; Level++)   GenMatRhs(time,Level,0);  // matrix
#if PRINT_TIME==1
  std::clock_t end_time=std::clock();
  std::cout << "  Assembly time -----> ="<< std::setprecision(7) << double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif

  /// C) Solution of the linear MGsystem (MGSolT::MGSolve).
  MGSolve(1.e-6,40);
#if PRINT_TIME==1
  end_time=std::clock();
  std::cout << " Assembly+solution time -----> ="<< std::setprecision(7) << double(end_time- start_time) / CLOCKS_PER_SEC
            << "s "<< std::endl;
#endif

  /// D) Update of the old solution at the top Level  (MGSolT::OldSol_update),
  x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);
//    compute_probeT(_mgphys._flag_probe,_mgphys._node_glob_probe);
  return;
}// =======================================================================================





// =====================================================
#ifdef TBK_EQUATIONS
// =====================================================
void  MGSolT::f_mu(double val[]) {

  // Turbulent viscosity
  if(_kappa_g[0]< 1.e-20) _kappa_g[0]= 1.e-20;  // kappa
  if(_kappaT_g[0]< 1.e-20) _kappaT_g[0]= 1.e-10;           // kappa theta
  if(_kappaT_g[1]< 1.e-20) _kappaT_g[1]= 1.e-10;           // epsilon theta
  double tau_k=1.;  // turbulent time

#if (TBK_EQUATIONS==0)   //kappa  (Prandtl length)--------------------------
  tau_k=_y_dist/sqrt(_kappa_g[0]);
#endif  // end kappa -------------------------------------------------------

#if (TBK_EQUATIONS/2==1)  // kappa-epsilon (TBK_EQUATIONS=2 or 3) ---------- 
  tau_k= CMU*_kappa_g[0]/_kappa_g[1]; // tau=1/omega=CMU*kappa/epsilon
#endif   // kappa-epsilon --------------------------------------------------

#if (TBK_EQUATIONS/2==2)       // kappa-omega ------------------------------
  if(_kappa_g[1]> 1.)   tau_k=1/_kappa_g[1];  // tau=1/omega

#endif  // -----------------------------------------------------------------
  // turbulent viscosity
  double Re_t= _kappa_g[0]*tau_k/_IRe;
  // limits
  if(Re_t > MU_TOP)  Re_t =MU_TOP;
  if(Re_t < MU_LOW)  Re_t =MU_LOW;
  _mu_turb=Re_t;

#ifdef TTBK_EQUATIONS
  double rT=sqrt(sqrt(_IPrdl))*_kappa_g[1]*_kappaT_g[0]/_kappaT_g[1]*_kappa_g[0];
  _IPrdl_turb=20.*rT/(PRT*(19*rT+0.4));

#endif

//   // Turbulent viscosity
//   if(_kappa_g< 1.e-20) _kappa_g= 1.e-20;
//   if(_omega_g< 0.1) _omega_g=0.1;
// //     double y_dist=val[2];double sP=val[3];
//   double Re_t= _kappa_g*_rhof/(_omega_g*_muf+1.e-20);
//
//   _mu_turb= Re_t;// val[4]=val[4]=
// #ifdef LOWRE
//   // eff visc val[6](k) val[7](w) (all eq) -------------------
//   _mu_turb= Re_t*(0.024+ 0.166666666667* Re_t)/(1.+0.166666666667* Re_t);// val[4]=
// #endif
// #ifdef SST
//   // Phi2 - F2 ---------------
//   double Phi2=500.*_muf/(_rhof*_y_dist*_y_dist*_omega_g);
//   const double Phib=2.*sqrt(_kappa_g)/(0.09*_omega_g*_y_dist);
//   if(Phi2 < Phib) Phi2=Phib; const  double F2=tanh(Phi2*Phi2);
//
//   // eff visc coeff  fact_mu (all eq)
//   double fact_mu=0.31*_omega_g/(sqrt(0.5*_sP)* F2+1.e-20);
//   if(fact_mu > 1) fact_mu=1.
//                             _mu_turb= Re_t*fact_mu*(0.024+ 0.166666666667* Re_t)/(1.+0.166666666667* Re_t);// val[4]
//
// #endif  // end SST ---------------------
  return;
}
// end ---   f_mu -------------------------
#endif // ----------  end TBK_EQUATIONS  

//set solver piecewice for temperature
MGSolTC::MGSolTC(MGEquationsMap& mg_equations_map_in,
                 int nvars_in[],
                 std::string eqname_in,
                 std::string varname_in):
  MGSolDA(mg_equations_map_in,nvars_in,eqname_in,varname_in),
  /// A) reading parameters
  // mesh params ------------
  _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes
  // phys params ------------
  _dt(_mgutils.get_par("dt")),      // parameter  dt
  _uref(_mgphys.get_par("Uref")),   // parameter  vel reference
  _lref(_mgphys.get_par("Lref")),   // parameter  length reference
  _Tref(_mgphys.get_par("Tref")),   // parameter  temperature reference
  _rhof(_mgphys.get_par("rho0")),   // parameter  density reference
  _muf(_mgphys.get_par("mu0")),     // parameter  viscosity reference
  _cp0(_mgphys.get_par("cp0")),     // parameter  Cp reference
  _kappa0(_mgphys.get_par("kappa0"))// parameter  conductivity reference
//   _h_conv(_mgphys.get_par("hconv")),  // parameter  heat transfer convection coefficient
//   _T_inf(_mgphys.get_par("T_inf"))
{
  //  =================================================
  /// B) setting class ariables
  _var_names[0]=varname_in;
  _refvalue[0]=_Tref;

  /// C ) setting solver type
  for(int l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVERT);

  /// D) setting nondimensional parameters
  _alpha=_kappa0/(_rhof*_cp0);
  _IPrdl=_rhof*_alpha/_muf;
  _IRe=_muf/(_rhof*_uref*_lref);
  _IPrdl_turb=PRT;
  _alpha_turb=0.;
//   _Nusselt=(_h_conv*_lref)/_kappa0;
  _qheat=_mgphys.get_par("qheat_sg")*_lref/(_rhof*_cp0*_Tref*_uref);
  _qs=_mgphys.get_par("qs")/(_rhof*_cp0*_Tref*_uref);

  return;
}








//  ====================================================
/// This function assembles the matrix and the rhs:
//  ====================================================
void  MGSolTC::GenMatRhs(
  const double /*time*/, // time  <-
  const int Level,  // Level <-
  const  int mode    // mode  <- (1=rhs) (0=only matrix)
) {  // ===============================================

  /// Set up
  // geometry ---------------------------------------------------------------------------------------
  const int  ndim = DIMENSION;                                           //dimension
  const int  offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
  const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element sides
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];                     // element connectivity
  int        el_neigh[NDOF_FEM];                                         // bd element connectivity
  int        sur_toply[NDOF_FEMB];                                       // boundary topology
  double     xx_qnds[DIMENSION*NDOF_FEM], xxb_qnds[DIMENSION*NDOF_FEMB]; // element node coords
  double xm[DIMENSION];
  // gauss integration  -----------------------------------------------------------------------------
  const int el_ngauss = _fe[2]->_NoGauss1[ndim-1];                   // elem gauss points
  const int elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];             // bd elem gauss points
  double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];             // Jac, Jac*w Jacobean
  double dphijdx_g[3][DIMENSION];
  double dphiidx_g[3][DIMENSION];   // global derivatives at g point

  // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
  int el_ndof[3];
  el_ndof[0]=1;
  int elb_ndof[3];
  elb_ndof[0]=1;   // number of el dofs
  int el_mat_nrows =0;                                               // number of mat rows (dofs)
  for(int ideg=2; ideg<3; ideg++) {
    el_ndof[ideg]=_fe[ideg]->_NoShape[ndim-1];
    elb_ndof[ideg]=_fe[ideg]->_NoShape[ ndim-2];
    el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
  };
  int el_mat_ncols = el_mat_nrows;                    // square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector

  // coupling  fields -------------------------------------------------------------------------------
  int T_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F]];    // Temperature
  double  val_tbg[3];                                        // turbulence
  double vel_g[DIMENSION];                                   // velocity
  for(int idim=0; idim<DIMENSION; idim++) vel_g[idim] =0.;   // velocity not coupled
  double Pe_h[3],f_upwind[3],h_eff[3];                       // local Peclet, upwind, h_eff

  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) ---------------------------
  A[Level]->zero();  if(mode ==1) b[Level]->zero();   // global matrix+rhs
  DenseMatrixM KeM;  DenseVectorM FeM;                // local  matrix+rhs
  KeM.resize(el_mat_nrows,el_mat_ncols);
  FeM.resize(el_mat_nrows);        // resize  local  matrix+rhs

  // reading from mesh 1 (reactor)
//    double t_out=780;// reading  //_data_eq[2].ub[T_idx*NDOF_FEM+lv_node];
//    double tmp2=673.15;
// #ifdef COUPLED_MESH
//             FILE *fp; fp = fopen("cmsh1.txt", "r");// open the file
//             fscanf(fp,"%lf,%lf", &t_out, &tmp2);   // read from file
//             fclose(fp);                            // close the file
// #endif

  // From reactor
  double u2_in=-1.6;double v2_in=-1.6;
  double w2_in=-1.6; double T2_in=608.15; double p2_in=0.;
// #ifdef COUPLED_MESH
//   FILE *fp;   fp = fopen("cmsh1.txt", "r"); // open the file
//   fscanf(fp,"%lf,%lf,%lf", &T2_in, &v2_in, &p2_in);     // read frm file
//   fclose(fp);                               // close the file
// #endif
  u2_in=_mgmesh._mesh2_data[0]; v2_in=_mgmesh._mesh2_data[1];
  w2_in=_mgmesh._mesh2_data[2]; p2_in=_mgmesh._mesh2_data[DIMENSION];
  T2_in=_mgmesh._mesh2_data[DIMENSION+1];
  vel_g[2]=-sqrt(u2_in*u2_in+v2_in*v2_in+w2_in*w2_in); //set velocity along z-axis
  



  /// b) Element  Loop over the volume (n_elem)
  const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
  const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
  for(int iel=0; iel < (nel_e - nel_b); iel++) {

    // set to zero matrix and rhs and center
    KeM.zero();    FeM.zero();

    // geometry and element  fields ------------------------------------
    // Element Connectivity (el_conn)  and coordinates (xx_qnds)
    _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);
    _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);

    // set element-nodes variables  bc (bc_q_dofs)
    get_el_dof_bc(Level,iel,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

    // grid and gaussian points ---------------------------------------------------------------
    for(int idim=0; idim<DIMENSION; idim++) {
      double sum=0.;
      for(int d=0; d< NDOF_FEM; d++) {
        _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];// element nodes xxg (DIM)
        sum += xx_qnds[idim*NDOF_FEM+d];
      }
      xm[idim]=sum/NDOF_FEM;
      // element grid distance
      h_eff[idim]=fabs(xx_qnds[idim*NDOF_FEM+2+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+0]);
      double h_1=fabs(xx_qnds[idim*NDOF_FEM+3+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+1]);
      if(h_eff[idim] <  h_1) h_eff[idim]=h_1;  // Max dx diagonal term
      f_upwind[idim]=0.;                       // upwind
    }
    // element field values --------------------------------------------------------
    for(int deg=0; deg<3; deg++) {  // loop over polinomial degree
      for(int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
        _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],
                                             el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub);
      }
    }

    //  external node quantities -------------------------------------
    //  external cell properties -------------------------------------

//     double peak_f = asby_power(xm[0],xm[1])*axial_pf(xm[0],xm[1],xm[2]);
//     double loss_f =lossf(xm[0],xm[1],xm[2])/axial_lf(xm[2]);
//
//             if (xm[2]>HIN &&  xm[2]<HOUT) {
// //                 std::cout << iel+nel_b << " " << xm[0] << " " << xm[1] << " " << xm[2] << " " << peak_f  << " " << loss_f << endl;
// //                 aux_data_pf[Level2][iel]=peak_f ;
// //                 aux_data_pl[Level2][iel]=loss_f;
// 		 _mgphys._sys_data[0][iel+nel_b]=peak_f;
// 		  _mgphys._sys_data[1][iel+nel_b]=loss_f;
//             }
//             else {
//
// 		_mgphys._sys_data[0][iel+nel_b]=peak_f=0.;
// 		_mgphys._sys_data[1][iel+nel_b]=loss_f=0.;
//             }



    double phase=1.;
//     double Qad =_qheat*_mgphys._sys_data[0][iel+nel_b];  //solo con flusso medio
    double Qad =0.;//-_data_eq[2].ub[t_idx*NDOF_FEM+t1_ind];
    double rho=1.;  double rhocp=1.;
    double qn =0;
    // ======================================================================
    // Volume =============================================================
    // ======================================================================

    // ---------------------------------------------
    /// c) gaussian integration loop (n_gauss)
    // --------------------------------------------
    for(int qp=0; qp< el_ngauss; qp++) {

      // shape functions at gaussian points -----------------------------------
      for(int ideg=2; ideg<3; ideg++) {  // linear-quadratic  [1]=linear [2]=quad
        det[ideg]      = _fe[ideg]->Jac(qp,xx_qnds,InvJac[ideg]);     // Jacobian
        JxW_g[ideg] =det[ideg]*_fe[ideg]->_weight1[ndim-1][qp];       // weight
        _fe[ideg]->get_phi_gl_g(ndim,qp,_phi_g[ideg]);               // shape funct
        _fe[ideg]->get_dphi_gl_g(ndim,qp,InvJac[ideg],_dphi_g[ideg]); // global coord deriv
      }

      double alpha_eff =_IPrdl*_IRe;

      //  fields -----------------------------------------------------------
      interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],
                    _phi_g[2],el_ndof[2],_ub_g[2]); // quadratic

      // upwind
      double mod2_vel=fabs(vel_g[2]); 
        double Pe_h=0.5*mod2_vel*h_eff[2]/_IRe;     //  local Peclet number
        f_upwind[2]=UP_WIND_T*rhocp*0.5*(1./tanh(Pe_h)-1./Pe_h)*h_eff[2]/sqrt(mod2_vel);    // f_upwind[] upwind
          
      
      // Temperature field (T_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F])---
#ifndef TEMPERATURE_DEPENDENCE
      double Temp_g=_ub_g[2][T_idx]; // local temp
      rho=densityT(Temp_g)/_rhof; //  density (densityT see MGconfclass.h)
      rhocp =rho;
//       alpha_eff *=kappaT(uold_g[0]);// nondimensional viscosity (mu/muref)
#endif
      // Velocity field -> [NS_F] -> (quad, _indx_eqs[NS_F]) -----------------------
//       double mod2_vel=1.e-20;

      /// d) Local (element) assemblying energy equation
      // *********************** *******************************
      for(int i=0; i<el_ndof[2]; i++)     {
        // set up row i
        const double phii_g=_phi_g[2][i];
        for(int idim=0; idim< ndim; idim++)  dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
        double dtxJxW_g=JxW_g[2]*_bc_vol[i];
 
        // Rhs Assemblying  ----------------------------------------------------------------------
        if(mode == 1) {         // rhs
          // vapour generator
          if(xm[2]< H_GV+LONG_GV && xm[2]> H_GV) { // (see Reactor.h)
// 		      Qad=-8.*_qheat*((_ub_g[2][T_idx]-673.15)>0 ? ((_ub_g[2][T_idx]-675.15)/145):0);
            Qad=-1.*Q_GV*(((_ub_g[2][T_idx]-608.15) >0.) ? _ub_g[2][T_idx]-608.15:0.);   //*exp(-1*time*30);
          }
          FeM(i) += dtxJxW_g*(
                      rhocp*_ub_g[2][T_idx]*phii_g/_dt   // time
                      +Qad*phii_g                       // heat source
                    );

        } // --------------------------------------------------------------------------------------
        // Matrix Assemblying ---------------------------------------------------------------------
    
    
       if( _ub_g[2][2] < H_GV){
	 
           KeM(i,i)=1;
           FeM(i)=608.15;
       }
       
       else{
	  for(int j=0; j<el_ndof[2]; j++) {
	    double phij_g= _phi_g[2][j];
	    double Adv=0.;
	    double Lap=0.;
  //           for (int idim=0; idim< ndim; idim++) {
  //             dphijdx_g[2][idim]=_dphi_g[2][j+idim*el_ndof[2]];
  //             Adv +=vel_g[idim]*dphijdx_g[2][idim];                     // advection
  //             Lap +=(alpha_eff+f_upwind[idim]*vel_g[idim]*vel_g[idim])*//(alpha_eff+f_upwind[idim]*vel_g[idim]*vel_g[idim])* // upwind
  //              dphijdx_g[2][idim]*dphiidx_g[2][idim];              // Laplacian
  //           }
	    for(int idim=0; idim< ndim; idim++) {
	      dphijdx_g[2][idim]=_dphi_g[2][j+idim*el_ndof[2]];
	      Adv +=vel_g[idim]*dphijdx_g[2][idim];                     // advection
	      Lap +=(alpha_eff+f_upwind[idim]*vel_g[idim]*vel_g[idim])*//(alpha_eff+f_upwind[idim]*vel_g[idim]*vel_g[idim])* // upwind
		    dphijdx_g[2][idim]*dphiidx_g[2][idim];              // Laplacian
	    }
	    // energy-equation
	    KeM(i,j) +=dtxJxW_g*(
			rhocp*phii_g*phij_g/_dt // time term
			+rhocp*Adv*phii_g      //adv
			+ Lap                  //diff
		      );
	  }
        }
      } // ----------------------------------------------------------------------------------------
    } // end of the quadrature point qp-loop ***********************



    // ======================================================================
    // ====================== boundary ======================================
    // ======================================================================

    for(int iside=0; iside< el_sides; iside++)  {
      if(el_neigh[iside] == -1) {

        double alpha_eff = _IPrdl*_IRe;

        for(int idof=0; idof<NDOF_FEMB; idof++) {
          sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
          int idofb=sur_toply[idof];
          elb_conn[idof]=el_conn[idofb];     // connectivity vector
          for(int idim=0; idim<DIMENSION; idim++) {
            xxb_qnds[idim*NDOF_FEMB+idof]=xx_qnds[idim*NDOF_FEM+idofb];   // coordinates
            _data_eq[2].ub[idim*NDOF_FEM+idof]=xxb_qnds[idim*NDOF_FEMB+idof];
          }
        }
        // ======================================================================
        // Dirichlet boundary conditions
        if(_bc_vol[sur_toply[NDOF_FEMB-1]] ==0) {

          //  [NDOF_FEMB-1] is the midpoint of a quadaratic FEM element (HEX9 (2D) or HEX27 (3D))
          int bc_s=(int)_bc_bd[sur_toply[NDOF_FEMB-1]];     // b cond
          double det= _fe[2]->JacSur(elb_ngauss-1,xxb_qnds);// jacobian
          double Ipenalty=det/_dt;                          // Dirichlet bc flag
          // local boundary loop   ---------------------------------------
          for(int lb_node=0; lb_node< elb_ndof[2]; lb_node++) {
            int lv_node= sur_toply[lb_node]; // local vol index

            // flag setup: ---------------------------------------------------------------------
            // surface flag ->  \int (bc_var(?1)*T+bc_val(1?)*val)ds
            int  bc_val = (int)((bc_s&2)>>1);  // (1?)non homogeneous
            int  bc_var = (int)(bc_s%2);    // (?1) variable
	    
           double val=608.15;
           if (bc_val==1) { val=T2_in;}
            // Assemblying  Matrix & rhs --------------------------------------------------------
            if(mode == 1) FeM(lv_node) += Ipenalty*val; // T2_in read from reactor mesh
            KeM(lv_node,lv_node) += Ipenalty;  //  Dirichlet bc
            // ---------------------------------------------------------------------------------
          }// lb_node -end  local boundary loop -------------------------
        } // end if Dirichlet  boundary conditions

        // ======================================================================
        // Neumann boundary conditions
        else if(_bc_bd[sur_toply[NDOF_FEMB-1]] !=0) {


          // gaussian integration loop (n_gauss)
          // -----------------------------------------------
          for(int qp=0; qp<  elb_ngauss; qp++) {

            // quad/linear  [2]=quad [1]=linear------------------------------------
            det[2]     = _fe[2]->JacSur(qp,xxb_qnds);    // local coord _phi_g and jac
            JxW_g[2]=det[2]*_fe[2]->_weight1[ndim-2][qp];// weight
            _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);   // global coord _phi_g
#ifdef AXISYM   // axisymmetric  
            JxW_g  *=_ub_g[0];
#endif
            // ***********************************
            // local side loop (over the node face)
            for(int lsi_node=0; lsi_node< elb_ndof[2]; lsi_node++) {

              // set up row i
              const double phii_g=_phi_g[2][lsi_node]; // boundary test function
              int lei_node= sur_toply[lsi_node]; // local element index

              // flag setup: ---------------------------------------------------------------------
              // volume flag ->   \int_vol T_eqs dv *bc_v
              int bc_v=(int)_bc_vol[lei_node];
              double dtxJxW_g=JxW_g[2]*bc_v;      // Neumann bc flag and Jac
              // surface flag ->  \int (bc_var(?1)*T+bc_val(1?)*val)ds
              int bc_s=(int)_bc_bd[lei_node];
              int  bc_val = (int)((bc_s&2)>>1);  // (1?) non-homogeneous
              int  bc_var = (int)(bc_s%2);       // (?1) tau=A*variable
              // ---------------------------------------------------------------------------------

              // Assemblying rhs ----------------------------------------------------------------
//               if (mode == 1){

// 		if(xxb_qnds[2*NDOF_FEM+NDOF_FEM-1]<0.001) qn=1.;
//                 FeM(lei_node) += bc_val*dtxJxW_g*phii_g*(// non-homogeneous
//                                      0.
//                                    bc_var*qn // flusso
//                                    +(1 -bc_var)*qn                   // Neumann heat flux bc
//                                  );
// 	      } // -------------------------------------------------------------------------------
              // Assemblying Matrix --------------------------------------------------------------
//               for (int lsj_node=0; lsj_node< elb_ndof[2];  lsj_node++) {
//                 KeM(lei_node,sur_toply[lsj_node]) += dtxJxW_g*bc_var*alpha_eff*phii_g*_phi_g[2][lsj_node]; // Robin bc  (k*dt/dn = h*(-T))
//               }// end j  -------------------------------------------------------------------------

            }// i   +++++++++++++++++++++++++++++++
          } // end of the quadrature point qp-loop **********************


        } // Neumann non homog

      } //end if side
    } // ======================  end for boundary =======================================

    /// e) Global assemblying energy equation
    A[Level]->add_matrix(KeM,el_dof_indices);                  // global matrix
    if(mode == 1)   b[Level]->add_vector(FeM,el_dof_indices);  // global rhs

  } // end of element loop
  // clean
  el_dof_indices.clear();
  A[Level]->close();
  if(mode == 1)  b[Level]->close();   //A[Level]->print();
#ifdef PRINT_INFO
  std::cout<< " Matrix Assembled(T)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif
//   A[Level]->close();
//   double nnnorm=A[Level]->linfty_norm();
//  std::cout<< "norm " << nnnorm << "\n";


  return;
}


// =========================================================================================
/// This function compute average surface  values on  over mesh[0] and
///  transfer data to mesh[1]
void MGSolTC::print_ext_data(
  double mesh2mesh_data[]   // <- data vector on mesh[0]
) {
// ========================================================================================

  int flag_probe=_mgphys._flag_probe;
  const int Level=_NoLevels-1;
  const int offset= _mgmesh._NoNodes[_NoLevels-1];
  // Compute temperature on surface defined in _mgphys._node_glob_probe
  double sum_u=0.;
  for(int inode=0; inode< flag_probe; inode++) {
    sum_u += (*x_old[Level])(_node_dof[Level][_mgphys._node_glob_probe[inode]])*_refvalue[0];
  }
  sum_u /=flag_probe;
  // Writing average temperature from surface over mesh[0] to mesh[1]
  _mgmesh._mesh2_data[5+DIMENSION+1]=sum_u;  // writing on mesh[0] (reactor)
  mesh2mesh_data[5+DIMENSION+1]=sum_u;// writing on mesh[1]        (primary loop)

#ifdef PRINT_INFO
  std::cout << " probe TC computed  " << Level << " (Level) with "  << flag_probe  << " (nodes) "  << std::endl;

   std::cout << " 0n TC probe_count " << flag_probe << std::endl;
   std::cout << " Velocity = " << mesh2mesh_data[5+0] << "  " << 
        mesh2mesh_data[5+1] << "  "<< mesh2mesh_data[5+2] << "  "<<std::endl;
  std::cout << " Pressure = " <<  mesh2mesh_data[5+DIMENSION]  << std::endl;
  std::cout << " Temperature = " <<  mesh2mesh_data[5+DIMENSION+1]  << std::endl;
 
#endif
  return;
}

#include <iomanip>
// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolTC::MGTimeStep(
  const double time,  // time
  const int /*iter*/  // Number of max inter
) {
// =========================================================================================

/// A) Set up the time step
  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;

  /// B) Assemblying of the Matrix-Rhs
#if PRINT_TIME==1
  std::clock_t start_time=std::clock();
#endif
  GenMatRhs(time,_NoLevels-1,1);                                                // matrix and rhs
  for(int Level = 0 ; Level < _NoLevels-1; Level++)   GenMatRhs(time,Level,0);  // matrix
#if PRINT_TIME==1
  std::clock_t end_time=std::clock();
  std::cout << "  Assembly time -----> ="<< std::setprecision(7) << double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif

  /// C) Solution of the linear MGsystem (MGSolT::MGSolve).
  MGSolve(1.e-6,40);
#if PRINT_TIME==1
  end_time=std::clock();
  std::cout << " Assembly+solution time -----> ="<< std::setprecision(7) << double(end_time- start_time) / CLOCKS_PER_SEC
            << "s "<< std::endl;
#endif

  /// D) Update of the old solution at the top Level  (MGSolT::OldSol_update),
  x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);
  return;
}// =======================================================================================


































#endif
#endif // personal application


