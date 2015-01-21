// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#include "ReactData.h"

#ifdef NS_EQUATIONS

// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNSC.h"       // Navier-Stokes class header file


// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include -----------------------------------
#include "MGMesh.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsMap.h"  // Equation map class

// standard lib -----------------------------------------------
#include <string.h>          // string library
#include <sstream>
// local alg lib -----------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ==============================================================

// ==================================================================
/// This routine constructs the FSI class:
MGSolNSC::MGSolNSC(
  MGEquationsMap& mg_equations_map_in,
  int             nvars_in[],
  std::string     eqname_in,
  std::string     varname_in
):  MGSolDA(mg_equations_map_in,nvars_in,eqname_in,varname_in),
  /// A) reading parameters
  _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes (top level)
  _dt(_mgutils.get_par("dt")),            // parameter  dt
  _uref(_mgphys.get_par("Uref")),         // parameter  u reference
  _lref(_mgphys.get_par("Lref")),         // parameter  l reference
  _rhof(_mgphys.get_par("rho0")),         // parameter density
  _muf(_mgphys.get_par("mu0")) {          // parameter viscosity
  // ================================================================

  /// B) setting class variables
  // class variable names
  _dir=0;
  _mgphys._flag_probe=0; //
  int iname=0;
  for(; iname<nvars_in[2]; iname++) {    //quadratic
    std::ostringstream ostr;
    ostr << "qua" <<iname+1; //use the string stream just like cout,
    _var_names[iname]=ostr.str();
    _refvalue[iname]=_uref;
  }
  for(; iname<(nvars_in[1]+nvars_in[2]); iname++) {    //linear
    std::ostringstream ostr;
    ostr << "p_tube";//"lin"<<iname-nvars_in[2]+1; //use the string stream just like cout,
    _var_names[iname]=ostr.str();
    _refvalue[iname]=_rhof*_uref*_uref;
  }
  for(; iname<(nvars_in[0]+nvars_in[1]+nvars_in[2]); iname++) {    //piecewise
    std::ostringstream ostr;
    ostr << "pie"<<iname-nvars_in[2]-nvars_in[1]+1; //use the string stream just like cout,
    _var_names[iname]=ostr.str();
    _refvalue[iname]=1;
  }


  /// C ) setting solver type
  for(int l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVERNS);

  /// D) setting nondimensional parameters
  _IRe=_muf/(_rhof*_lref*_uref);          // Reynolds number
  _IFr=9.81*_lref/(_uref*_uref);          // Froud number
  _dirg[0] = _mgphys.get_par("dirgx");    // x-gravity
  _dirg[1] = _mgphys.get_par("dirgy");    // y-gravity
  _dirg[2] = _mgphys.get_par("dirgz");    // z-gravity

  return;
}//  =================================================================


///  ====================================================
/// This function assembles the matrix and the rhs:
///  ====================================================
void  MGSolNSC::GenMatRhs(
  const double time,   // time  <-
  const int    Level,  // Level <-
  const int    mode    // mode  <- (1=rhs+matrix) (0=only matrix)
) {  // ===============================================

  /// Set up
  // geometry and bc---------------------------------------------------------------------------------
  const int ndim = DIMENSION;                         // dimension
  const int offset = _mgmesh._NoNodes[_NoLevels-1];   // mesh nodes
  const int el_sides= _mgmesh._GeomEl._n_sides[0];    // element sides
  int       el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];   // element connectivity
  int       el_neigh[NDOF_FEM];                       // bd element connectivity
  int       sur_toply[NDOF_FEMB];                     // boundary topology
  double    xx_qnds[DIMENSION*NDOF_FEM];              // element node coords
  double    xxb_qnds[DIMENSION*NDOF_FEMB];            // boundary of the element node coords
  int       _bc_vol[NDOF_FEM*_n_vars];                // element  b.cond flags (Neu or Dir)
  int       _bc_bd[NDOF_FEM*_n_vars];                 // element  b.cond flags (different possibilities)

  // gauss integration  -----------------------------------------------------------------------------
  const int el_ngauss = _fe[2]->_NoGauss1[ndim-1];         // quadratic element gauss points
  const int elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];   // quadratic bd elem gauss points
  double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];   // determinant of Jacobean and det*gauss_weigth
  double dphijdx_g[3][DIMENSION];                          // global derivatives at gauss point
  double dphiidx_g[3][DIMENSION];                          // global derivatives at gauss point
  double normal[DIMENSION];  double h_eff[DIMENSION];
  h_eff[0]=0;  h_eff[1]=0;  h_eff[2]=0;
  const int t_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F]];  //Temperature quiiiiii
  double pump_work=_mgphys.get_par("pump_block");//take into account if user want to activate the pump
  
  double     xm[DIMENSION];                          // mid element point


  // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
  int elb_ndof[3];  elb_ndof[0]=NDOF_BK;
  elb_ndof[1]=NDOF_PB;  elb_ndof[2]=NDOF_FEMB;  // number of element boundary dofs
  int el_mat_nrows =0;                                              // number of matrix rows (dofs)
  for(int ideg=0; ideg<3; ideg++) el_mat_nrows +=_nvars[ideg]*_el_dof[ideg];
  int el_mat_ncols = el_mat_nrows;                    // square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector

  // coupling  fields -------------------------------------------------------------------------------
  int ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[DA_F]]; // Quadratic
  double _ub_g[3][14];                                     // values of external fields
  int probe_count=0;
  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix+rhs) ---------------------------
  A[Level]->zero();  if(mode ==1) b[Level]->zero();                 // global matrix A and rhs b
  DenseMatrixM KeM;  DenseVectorM FeM;                              // local  matrix KeM and rhs FeM
  KeM.resize(el_mat_nrows,el_mat_ncols);
  FeM.resize(el_mat_nrows); // resize local matrix and rhs

  // number of total elements for level
  int ndof_lev=0;
  for(int pr=0; pr <_mgmesh._iproc; pr++) {
    int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
    ndof_lev +=delta;
  }
  // From reactor
  double u2_in=-1.6; double v2_in=-1.6;
  double w2_in=-1.6; double T2_in=608.15; double p2_in=608.15;
  
// #ifdef COUPLED_MESH
//   FILE *fp;   fp = fopen("cmsh1.txt", "r"); // open the file
//   fscanf(fp,"%lf,%lf,%lf", &T2_in, &v2_in, &p2_in);     // read frm file
//   fclose(fp);                               // close the file
// #endif
  u2_in=_mgmesh._mesh2_data[0]; v2_in=_mgmesh._mesh2_data[1];
  w2_in=_mgmesh._mesh2_data[2]; p2_in=_mgmesh._mesh2_data[DIMENSION];
  T2_in=_mgmesh._mesh2_data[DIMENSION+1];
  double v_mono=-sqrt(u2_in*u2_in+v2_in*v2_in+w2_in*w2_in); //set velocity along z-axis

  /// b) Element  Loop over the volume (n_elem)
  const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
  const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
  for(int iel=0; iel < (nel_e - nel_b); iel++) {

    KeM.zero();    FeM.zero();      // set to zero matrix and rhs

    // geometry and element  fields ------------------------------------

    _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);   // Element Connectivity (el_conn) and coordinates (xx_qnds)
    _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh); // Neighbors of the element

    // set element-nodes variables  bc (bc_q_dofs)
    get_el_dof_bc(Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

      for(int idim=0; idim<DIMENSION; idim++) {
      double somma=0.;                                             // element middle point
      for(int d=0; d< NDOF_FEM; d++) {
        _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];// element nodes xxg (DIM)
        somma += xx_qnds[idim*NDOF_FEM+d];
      }
      xm[idim]=somma/NDOF_FEM;                                     // element middle point

      }
    
    // element field values
    for(int eq=0; eq<_data_eq[2].n_eqs; eq++) {
      _data_eq[2].mg_eqs[eq]->get_el_sol(
	0,_data_eq[2].indx_ub[eq+1]-_data_eq[2].indx_ub[eq],
        _el_dof[2],el_conn,offset,_data_eq[2].indx_ub[eq],_data_eq[2].ub);
    }
    //linear field
    _data_eq[1].mg_eqs[0]->get_el_sol(
      _nvars[2],_data_eq[1].indx_ub[0],_el_dof[1],el_conn,offset,0,_data_eq[1].ub);
    //  external node quantities -------------------------------------
    _data_eq[0].mg_eqs[0]->get_el_sol_piece(
      _nvars[2]+_nvars[1],_data_eq[0].indx_ub[0],_el_dof[0],iel+ndof_lev,offset,0,_data_eq[0].ub);


    // ======================================================================
    // Volume =============================================================
    // ======================================================================

    // ---------------------------------------------
    /// c) gaussian integration loop (n_gauss)
    // --------------------------------------------
    for(int qp=0; qp< el_ngauss; qp++) {

      // shape functions at gaussian points -----------------------------------
      for(int ideg=1; ideg<3; ideg++) {  // linear-quadratic  [1]=linear [2]=quad
        det[ideg]      = _fe[ideg]->Jac(qp,xx_qnds,InvJac[ideg]);     // Jacobian
        JxW_g[ideg] =det[ideg]*_fe[ideg]->_weight1[ndim-1][qp];       // weight
        _fe[ideg]->get_phi_gl_g(ndim,qp,_phi_g[ideg]);               // shape funct
        _fe[ideg]->get_dphi_gl_g(ndim,qp,InvJac[ideg],_dphi_g[ideg]); // global coord deriv
      }
      JxW_g[0]=JxW_g[2];
      _fe[0]->get_phi_gl_g(ndim,qp,_phi_g[0]);               // shape function piecewise

      //  fields -----------------------------------------------------------
      interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],
                    _phi_g[2],_el_dof[2],_ub_g[2]); // quadratic
      interp_el_sol(_data_eq[1].ub,0,_data_eq[1].indx_ub[0],_phi_g[1],_el_dof[1],_ub_g[1]); // linear
      interp_el_sol(_data_eq[0].ub,0,_data_eq[0].indx_ub[0],_phi_g[0],_el_dof[0],_ub_g[0]); //constant


      /// d) Local (element) assemblying DA equation
      // *********************** *******************************
      for(int i=0; i<_el_dof[2]; i++)     {    //  --- QUADRATIC ---
        // set up row i
        const double phii_g=_phi_g[2][i];
        for(int idim=0; idim< ndim; idim++)  dphiidx_g[2][idim]=_dphi_g[2][i+idim*_el_dof[2]];

        for(int ivar=0; ivar<_nvars[2]; ivar++) {
          double dtxJxW_g=JxW_g[2]*_bc_vol[i+ivar*_el_dof[2]];
          int index=i+ivar*_el_dof[2];

          // Rhs Assemblying  ----------------------------------------------------------
//           if(mode == 1) {
            // rhs
            FeM(index) = (v_mono
//                             +_ub_g[2][DIMENSION]*phii_g/_dt   // time
//                             +_ub_g[2][2]*phii_g                       // heat source
                          );
//           }
           KeM(index,i+ivar*_el_dof[2]) =1.;
          // Matrix Assemblying ---------------------------
//           for(int j=0; j<_el_dof[2]; j++) {
//             double phij_g= _phi_g[2][j];
//             double Lap=0.;
//             for(int idim=0; idim< ndim; idim++) {
//               dphijdx_g[2][idim]=_dphi_g[2][j+idim*_el_dof[2]];
//               Lap +=dphijdx_g[2][idim]*dphiidx_g[2][idim];            // Laplacian
//             }
// 
//             // energy-equation
//             KeM(index,j+ivar*_el_dof[2]) +=dtxJxW_g*(0.
// //                                              phii_g*phij_g// time term
// //                                                             +phii_g*dphijdx_g[2][2]
// //                                               + 0.001*Lap                   //diff
//                                            );
//           }
        } // quadratic variable cycle
      } // ----------------------------------------END QUADRATIC


//       for(int i=0; i<_el_dof[1]; i++)     {  //    --- LINEAR ---
//         // set up row i
//         const double phii_g=_phi_g[1][i];
//         for(int idim=0; idim< ndim; idim++)  dphiidx_g[1][idim]=_dphi_g[1][i+idim*_el_dof[1]];
//
//         for(int ivar=0; ivar<_nvars[1]; ivar++) {
//
//           double dtxJxW_g=JxW_g[1]*_bc_vol[i+(ivar+_nvars[2])*NDOF_FEM];
//           int index=i+ivar*_el_dof[1]+_el_dof[2]*_nvars[2];
//
//           // Rhs Assemblying  ----------------------------------------------------------
//           if(mode == 1) {
//             // rhs
//             FeM(index) = 0*dtxJxW_g*(
//                            _ub_g[1][ivar]*phii_g/_dt   // time
//                            +1.*phii_g                  // heat source
//                          );
//           }
//
//           // Matrix Assemblying ---------------------------
//           for(int j=0; j<_el_dof[1]; j++) {
//             double phij_g= _phi_g[1][j];
//             int jndex=j+_el_dof[2]*_nvars[2]+ivar*_el_dof[1];
//             double Lap=0.;
//             for(int idim=0; idim< ndim; idim++) {
//               dphijdx_g[1][idim]=_dphi_g[1][j+idim*_el_dof[1]];
//               Lap +=dphijdx_g[1][idim]*dphiidx_g[1][idim];            // Laplacian
//             }
//
//             // energy-equation
// //                         KeM(index,jndex) +=dtxJxW_g*(
// //                                                                phii_g*phij_g/_dt // time term
// //                                                                + Lap                   //diff
// //                                                            );
//             KeM(index,index) = 0+0*dtxJxW_g*(
//                                  phii_g*phij_g/_dt // time term
//                                  + 0*Lap                   //diff
//                                  +dphijdx_g[1][2]*phii_g
//                                );
//           }
//         } //// linear variable cycle
//       } // ----------------------------------END LINEAR


//       for(int i=0; i<_el_dof[0]; i++)     {  //    --- Piecewise ---
//         // set up row i
//         const double phii_g=_phi_g[0][i];
//         for(int ivar=0; ivar<_nvars[0]; ivar++) {
//
//           double dtxJxW_g=JxW_g[0]*_bc_vol[i+(ivar+_nvars[1]+_nvars[2])*NDOF_FEM];
//           int index=i+ivar*_el_dof[0]+_el_dof[2]*_nvars[2]+_el_dof[1]*_nvars[1];
//
//           // Rhs Assemblying  ----------------------------------------------------------
//           if(mode == 1) {
//             // rhs
//             FeM(index) += dtxJxW_g*(
// // 			                  _ub_g[2][0]*phii_g
//                             _ub_g[0][ivar]*phii_g/_dt   // time
//                             +_ub_g[2][0]*phii_g         //source x
//                             +_ub_g[2][1]*_ub_g[2][2]*phii_g*ivar   //source x*y*z
//                           );
//           }
//           // Matrix Assemblying ---------------------------
//           for(int j=0; j<_el_dof[0]; j++) {
//             double phij_g= _phi_g[0][j];
//             int jndex=j+ivar*_el_dof[0]+_el_dof[2]*_nvars[2]+_el_dof[1]*_nvars[1];
//
//             KeM(index,jndex) +=dtxJxW_g*(phii_g*phij_g/_dt);// time term
//           }
//         } ////  variable cycle
//       } // ----------------------------------END piecewise

    } // end of the quadrature point qp-loop ***********************
    double pos_ind[NDOF_FEMB/2];
    double neg_ind[NDOF_FEMB/2];
    for(int tmp=0; tmp<NDOF_FEMB/2; tmp++) {
      pos_ind[tmp]=-2;
      neg_ind[tmp]=-2;
    }

    // ======================================================================
    // ====================== boundary ======================================
    // ======================================================================
    int begin=0;
    for(int iside=0; iside< el_sides; iside++)  {
      int bc_tmp=0;

      if(el_neigh[iside] == -1) {

        for(int idof=0; idof<NDOF_FEMB; idof++) {                          //idof -> boundary node
          sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside]; //use map to find global node
          int idofb=sur_toply[idof];                                       //idofb -> element node
          elb_conn[idof]=el_conn[idofb];                                   //connectivity vector
          for(int idim=0; idim<DIMENSION; idim++) {
            xxb_qnds[idim*NDOF_FEMB+idof]=xx_qnds[idim*NDOF_FEM+idofb];    //get boundary coordinates
            _data_eq[2].ub[idim*NDOF_FEM+idof]=xxb_qnds[idim*NDOF_FEMB+idof];
          }
        }
        for(int iql=2; iql<3; iql++)           // --- QUADRATIC ---
          for(int ivar=0; ivar< _nvars[iql]; ivar++)    {
            // Dirichlet boundary conditions  ***********************************
            if(_bc_vol[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] ==0) {

              //[NDOF_FEMB-1] is the midpoint of a quadratic FEM element (HEX9 (2D) or HEX27 (3D))
              int bc_s=(int)_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM];     // b cond
              double det= _fe[iql]->JacSur(elb_ngauss-1,xxb_qnds);// jacobian
              double Ipenalty=det/_dt;                               // Dirichlet bc flag
              // local boundary loop   ---------------------------------------
              for(int lb_node=0; lb_node< elb_ndof[iql]; lb_node++) {
                int index= sur_toply[lb_node]+ivar*NDOF_FEM; // local vol index
                // flag setup (\int bc_var*T+bc_var*val)
                int  bc_val = (int)((bc_s&2)>>1);  // (1?)non homogeneous
                int  bc_var = (int)(bc_s%2);       // (?1) variable
                // Assemblying  Matrix & rhs
                if(mode == 1) {

//                   FeM(index) += bc_val*Ipenalty*_data_eq[iql].ub[ns_idx*NDOF_FEM+index];
                  FeM(index) += Ipenalty*(2.1);



                }
                KeM(index,index) += Ipenalty;  //  Dirichlet bc
              }// lb_node -end  local boundary loop -------------------------
            } // end if Dirichlet  boundary conditions
            // **********************************************************************

            else if(_bc_vol[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] !=0) {
                 int ss=1;
//               // Non homogenous Neumann boundary conditions  ***********************************
//
//               // gaussian integration loop (n_gauss)
//               // -----------------------------------------------
//               for(int qp=0; qp<  elb_ngauss; qp++) {
//
//                 // quad/linear  [2]=quad [1]=linear------------------------------------
//                 det[iql]  = _fe[iql]->JacSur(qp,xxb_qnds);    // local coord _phi_g and jac
//                 JxW_g[iql]=det[iql]*_fe[iql]->_weight1[ndim-2][qp];// weight
//                 _fe[iql]->get_phi_gl_g(ndim-1,qp,_phi_g[iql]);   // global coord _phi_g
//
//
// #ifdef AXISYM   // axisymmetric  (index ->0)
//                 interp_el_sol(_data_eq[iql].ub,0,DIMENSION,_phi_g[iql],elb_ndof[iql],_ub_g[iql]);
//                 JxW_g[iql]  *=_ub_g[iql][0];
//
// #endif
//                 // local side loop (over the node face)
//                 for(int lsi_node=0; lsi_node< elb_ndof[iql]; lsi_node++) {
//
//                   // set up row i
//                   const double phii_g=_phi_g[iql][lsi_node]; // boundary test function
//                   int index= sur_toply[lsi_node]+ivar*NDOF_FEM; // local element index
//                   int bc_s=(int)_bc_bd[index];
//                   int bc_v=(int)_bc_vol[index];
//                   double dtxJxW_g=JxW_g[iql]*bc_v; // Neumann bc flag and Jac
//
//                   // flag setup: +++++++++++++++++++++++++++++++++++++++
//                   //  \int (bc_var*T+bc_val*val)ds
//                   int  bc_val = (int)((bc_s&2)>>1);  // (1?) non-homogeneous
//                   int  bc_var = (int)(bc_s%2);       // (?1) tau=A*variable
//                   // ++++++++++++++++++++++++++++++++++++++++++++++++++++
//
//                   // Assemblying rhs ----------------------------
//                   if(mode == 1) {
// //                     FeM(index) += bc_val*(1-bc_var)*dtxJxW_g*phii_g*1.;
//                   }
//
//                   // Assemblying Matrix ---------------------------------
//                   for(int lsj_node=0; lsj_node< elb_ndof[iql];  lsj_node++) {
//                     int jndex=sur_toply[lsj_node]+ivar*NDOF_FEM;
//                     KeM(index,jndex) += dtxJxW_g*bc_var*phii_g*_phi_g[iql][lsj_node]; // Robin bc  (k*dt/dn = h*(-T))
//                   }// end j  ---------------------------------------
//
//                 }// i   +++++++++++++++++++++++++++++++
//               } // end of the quadrature point qp-loop **********************


            } // Neumann non homog

          }    // ---END QUADRATIC ---


//        } //end if side


        for(int iql=1; iql<2; iql++)  {      // --- LINEAR ---
          for(int ivar=0; ivar< _nvars[iql]; ivar++)    {
            // Dirichlet boundary conditions  ***********************************
//         int dn_flag=0;
//            for (int lb_node=0; lb_node< elb_ndof[iql]; lb_node++)
// 	     if (_bc_vol[sur_toply[lb_node]+ivar*_el_dof[iql]+_nvars[2]*NDOF_FEM] ==1) dn_flag++;


// 	     if (dn_flag< 2 ){ // Dirichlet  boundary conditions <-- ONLY THIS
            //[NDOF_FEMB-1] is the midpoint of a quadratic FEM element (HEX9 (2D) or HEX27 (3D))

            double det= _fe[iql]->JacSur(elb_ngauss-1,xxb_qnds);// jacobian
            double Ipenalty=det/_dt;                               // Dirichlet bc flag
            // local boundary loop   ---------------------------------------

            for(int lb_node=0; lb_node< elb_ndof[iql]; lb_node++) {


              if(_bc_vol[sur_toply[lb_node]+(ivar+_nvars[2])*NDOF_FEM] ==0) {
                begin=1;
                int index= sur_toply[lb_node]+ivar*_el_dof[iql]+_nvars[2]*NDOF_FEM; // local vol index
                // flag setup (\int bc_var*T+bc_var*val)
                int bc_s=(int)_bc_bd[index];     // b cond
                int  bc_val = (int)((bc_s&2)>>1);  // (1?)non homogeneous
                int  bc_var = (int)(bc_s%2);       // (?1) variable
                // Assemblying  Matrix & rhs
                if(bc_var == 0) {

//                               FeM(index) += (1-bc_var)*Ipenalty*_data_eq[iql].ub[ivar*NDOF_FEM+sur_toply[lb_node]];
                  FeM(index) += bc_var*Ipenalty*(0.);
                  if(neg_ind[3]<0) {
                    neg_ind[bc_tmp]=index; //pos_ind[bc_tmp]=index;
                    bc_tmp++;
                  }

                }
                if(bc_var == 1) {

//                               FeM(index) += (1-bc_var)*Ipenalty*_data_eq[iql].ub[ivar*NDOF_FEM+sur_toply[lb_node]];
                  FeM(index) += bc_var*Ipenalty*(0.);
                  if(pos_ind[3]<0) {
                    pos_ind[bc_tmp]=index; //pos_ind[bc_tmp]=index;
                    bc_tmp++;
                  }

                }
//                             KeM(index,index) += Ipenalty;  //  Dirichlet bc
              }//end if on bc_vol



              else if(_bc_vol[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] !=0 && xxb_qnds[3*9-1]<Z_LOW) {
                int index= sur_toply[lb_node]+ivar*_el_dof[iql]+_nvars[2]*NDOF_FEM;
                pos_ind[bc_tmp]=index;
                bc_tmp++;
              }
            }// lb_node -end  local boundary loop -------------------------
//          } // end if Dirichlet  boundary conditions
            // **********************************************************************



          } //end for variables
        }        // ---END LINEAR ---


//         for(int iql=0; iql<1; iql++)  {      // --- PIECEWISE ---
//           for(int ivar=0; ivar< _nvars[iql]; ivar++)    {
//             // Dirichlet boundary conditions  ***********************************
//             double det= _fe[2]->JacSur(elb_ngauss-1,xxb_qnds);// jacobian
//             double Ipenalty=det/_dt;                               // Dirichlet bc flag
//             if(_bc_vol[(ivar+_nvars[1]+_nvars[2])*NDOF_FEM]==0) {
//               int index= ivar*_el_dof[iql]+_nvars[1]*_el_dof[1]+_nvars[2]*_el_dof[2]; // local vol index
// 
//               // Assemblying  Matrix & rhs
//               if(mode == 1) {
//                 FeM(index) += Ipenalty*0.;//_data_eq[iql].ub[ivar*NDOF_FEM];
//               }
//               KeM(index,index) += Ipenalty;  //  Dirichlet bc
//             }//end if on bc_vol
//           } //end loop variables
//         } // -- END PIECEWISE --


      } //end if side

    } // ======================  end for boundary =======================================



    // ======================  p solver=======================================
    double hdir=0;

    for(int iql=1; iql<2; iql++)  {      // --- LINEAR ---
      for(int iside=0; iside< el_sides; iside++)  {
        if(el_neigh[iside] != -1) {

          for(int idof=0; idof<NDOF_FEMB; idof++) {                          //idof -> boundary node
            sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];//use map to find global node
            int idofb=sur_toply[idof];                                      //idofb -> element node
            elb_conn[idof]=el_conn[idofb];                                  //connectivity vector
            for(int idim=0; idim<DIMENSION; idim++) {
              xxb_qnds[idim*NDOF_FEMB+idof]=xx_qnds[idim*NDOF_FEM+idofb]; //get boundary coordinates
              _data_eq[2].ub[idim*NDOF_FEM+idof]=xxb_qnds[idim*NDOF_FEMB+idof];
            }
          }

          _fe[2]->normal_g(xxb_qnds,normal);
          double direction[DIMENSION];

          double loc_direction[DIMENSION];
          double dot=0;   double ddir=0;

          // 1D direction
          direction[0]=0.; direction[1]=0.;  direction[2]=-1.; //+direction-
          // flow direction ddir (+1 or -1)
          for(int idim=0; idim<ndim; idim++) {
            loc_direction[idim]= xx_qnds[(1+idim)*NDOF_FEM-1]-xxb_qnds[(1+idim)*NDOF_FEMB-1];
            ddir+=direction[idim]*loc_direction[idim];
          }
          if(ddir>0) ddir=1;   else ddir=-1;

          // element length    h_eff[]
          for(int idim=0; idim<ndim; idim++) {
            h_eff[idim]=fabs(xx_qnds[0+NDOF_FEM*idim]-xx_qnds[6+NDOF_FEM*idim]);//opposite nodes of the element
          }
          // Projection over the direction
          for(int idim=0; idim<ndim; idim++) {
            dot += direction[idim]*normal[idim];// direction projection
            hdir += direction[idim]*h_eff[idim];// element length
          }

// 	if (xxb_qnds[17]<0.9 && xxb_qnds[26]>1.){
// 	 direction[0]=0.; direction[1]=1.; direction[2]=0.;}
// 	 if (xxb_qnds[17]<0.9 && xxb_qnds[26]<1.){
// 	 direction[0]=0.; direction[1]=-1.; direction[2]=0.;}
// 	  if (xxb_qnds[26]>0.9 ){
// 	 direction[0]=0.; direction[1]=0; direction[2]=-1;}
//

          if(fabs(dot)>0.9) {
            for(int i=0; i< elb_ndof[1]; i++) {  // +++++++++++++++++++++++
//                             const int  indx_sol=sur_toply[i]+1*_el_dof[2];// index;
              const int  indx_row=sur_toply[i]+1*_el_dof[2];  //index ;
              // Assemblying Matrix ---------------------------------
// 	      FeM(indx_row)=-0.5;
              if(dot>0) { //KeM(indx_row,indx_row) = 1/fabs(hdir);
                pos_ind[i]=indx_row;
// 		  h_eff[2]+=xxb_qnds[2*NDOF_FEMB+i]/ elb_ndof[1];
              } else { // KeM(indx_row,indx_row) += -2*1/fabs(hdir);
                neg_ind[i]=indx_row;
// 		h_eff[2]+=-xxb_qnds[2*NDOF_FEMB+i]/ elb_ndof[1];
              }  // end dot>0
            }      // i-loop
          }         // end dot>0.9

        } //end if side

      } // ======================  end of p solver=======================================

    }
//         double rho1,rho2,t1,ta1,ta2,rho;
//     double t1=_data_eq[2].ub[t_idx*NDOF_FEM];
    double ta1=_data_eq[2].ub[t_idx*NDOF_FEM];
    double ta2=_data_eq[2].ub[t_idx*NDOF_FEM+6];//673.15;//_data_eq[2].ub[t_idx*NDOF_FEM];
    double dens=densityNS(608.15);
    double rho=0.5*(densityNS(ta1)+densityNS(ta2))/_rhof;
    int off_nod;//         int off_nod = 1000;
    double z_cord=xx_qnds[3*NDOF_FEM -1];

    for(int tmp=0; tmp<NDOF_FEMB/2; tmp++) {
// 	    off_nod=neg_ind[tmp]-NDOF_FEM; //position of node
// 	    t2=_data_eq[2].ub[t_idx*NDOF_FEM+off_nod];
// 	    rho2=(11367-1.1944*t1);rho1=(11367-1.1944*t1);
// 	   rho=((11367-1.1944*ta1)+(11367-1.1944*ta2))*0.5;
// 	    	    KeM(pos_ind[tmp],neg_ind[tmp])=-1/(0.5*(fabs(hdir)*rho2));//-(begin)/(fabs(hdir)*rho2);
      KeM(pos_ind[tmp],neg_ind[tmp])=-1/(0.5*(fabs(hdir)));//-(begin)/(fabs(hdir)*rho2);

// 	    	    	  printf("\n off diag   utili %f   nonut  %f \n    ", rho2,(11367-1.1944*ta1) );
// 	  }
      off_nod=pos_ind[tmp]-NDOF_FEM; //position of node
// 	      t1=_data_eq[2].ub[t_idx*NDOF_FEM+off_nod];
// 	 KeM(pos_ind[tmp],pos_ind[tmp])=1/(0.5*(fabs(hdir)*rho1));
      KeM(pos_ind[tmp],pos_ind[tmp])=1/(0.5*(fabs(hdir)));

// 	 printf("\n on diag   ut %f   nut  %f h  %f tut %f   tnut %f", rho1,(11367-1.1944*ta2),a ,t1, ta2);
// 	ta2=_data_eq[2].ub[t_idx*NDOF_FEM+t2_ind];
// 	rho1=(11367-1.1944*t1);
// 	rho2=(11367-1.1944*t2);

//====================================================================================================================================================
        double H_PUMP =  H_GV+LONG_GV - LONG_PUMP;
	double gravity_f = 0.;
            
	// gravity force h0
	

        // gravity force h1 and h2
        if(z_cord > H_GV+LONG_GV) gravity_f=-1.;
	if(z_cord > 13.485-0.4) gravity_f=0.;

        // Force in the porous media
        FeM(pos_ind[tmp])= ((begin==1)? 2:1)*(
	                  -9.81*(rho-1.)*_dirg[2]*gravity_f
                          +pump_work*
                          ((xm[2]>H_PUMP && xm[2]<H_PUMP+LONG_PUMP)?
                           1./(((int)(LONG_PUMP/(fabs(hdir))+0.5))*fabs(hdir)):0)
                         );
//====================================================================================================================================================


// 	    double dH;
// 	    double gravity_f;
// 	    double b_down=(13.485-4.485)*0.5+4.485 ;
//             double dim_pump =1.4;
//
// 	    dH = 0.;

//          if (a>13.485-dH && a<=13.485) gravity_f=1.;
// 	    else if (a<=C_OUT+dH) gravity_f=-1.;
// 	    else gravity_f=0.;

// 	    if( a>C_OUT) gravity_f=0.;
//
// 	    FeM(pos_ind[tmp]) = 9.81*_dirg[2]*gravity_f
//                                 +((a>=b_down && a<=b_down+dim_pump)? 36350.81975/(((int)(dim_pump/(fabs(hdir))+0.5))*fabs(hdir)*rho):0)
// 			        *pump_working;





    }


// #ifdef COUPLED_MESH
// // 	 printf("iproc    %d iel    %d offset %d zeta    %f\n", _iproc, iel , offset, xx_qnds[3*NDOF_FEM-1]);
// 
// // 	   	_data_eq[0].mg_eqs[0]->get_el_sol_piece(_nvars[2]+_nvars[1],_data_eq[0].indx_ub[0],_el_dof[0],iel+ndof_lev,offset,0,_data_eq[0].ub);
// // int aaa[1];
// // 	_data_eq[0].mg_eqs[0]->get_el_oldsol(_nvars[2]+_nvars[1],_nvars[0],_el_dof[0],iel+ndof_lev,offset,0,_data_eq[0].ub);
// // 	 get_el_sol_piece(_nvars[2]+_nvars[1],_data_eq[0].indx_ub[0],_el_dof[0],iel+ndof_lev,offset,0,_data_eq[0].ub);
// 
// 
// // 	kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+(ivar+ivar0)*offset]; // dof from top level
// //       uold[ id +(ivar+kvar0)*NDOF_FEM]= ((*x_old[_NoLevels-1])(kdof_top));
// // 	 printf("SCRITTO iproc    %d iel    %d offset %d     nel_e    %d      \n", _iproc, iel , offset, nel_b);
// // 	 ofstream myfile;
// // 	 myfile.open ("msh2.txt");
// // 	 myfile << "Temperature  Pressure \n"<<_data_eq[2].ub[t_idx*NDOF_FEM+NDOF_FEM-1]<<" "<<at;
// //          myfile.close();
// 
// 
//     if(xx_qnds[3*NDOF_FEM-1]<0.1) {
//       double pressure_out=((*x_old[_NoLevels-1])(_node_dof[_NoLevels-1][el_conn[6]+(1)*offset]));
// 
//       // print Temperature and Pressure into file DA_f
//       FILE *cmyfile;  cmyfile =fopen("cmsh2.txt","w");
//       fprintf(cmyfile, "%f,%f,%f\n Temperature Velocity Pressure",
//               _data_eq[2].ub[t_idx*NDOF_FEM+NDOF_FEM-1],
//               _data_eq[2].ub[ns_idx*NDOF_FEM+NDOF_FEM-1],pressure_out);
//       fclose(cmyfile);
// 
//     }
// #endif

















//             std::cout<<"\n\n\n KeM \n\n\n";
//     std::cout << FeM << endl;

    /// e) Global assemblying energy equation
    A[Level]->add_matrix(KeM,el_dof_indices);                  // global matrix
    if(mode == 1)   b[Level]->add_vector(FeM,el_dof_indices);  // global rhs
  } // end of element loop
  // clean
  el_dof_indices.clear();
  A[Level]->close(); if(mode == 1) b[Level]->close();
  //     A[Level]->print();  b[Level]->print();
#ifdef PRINT_INFO
  std::cout<< " Matrix Assembled(DA)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif

  return;
}



///  ====================================================
/// This function assembles the matrix and the rhs:
///  ====================================================
// ========================================================================
void  MGSolNSC::get_nod_probe(
) {  // ===============================================

  /// Set up
  int Level=_NoLevels-1;
  int mode=1;
  // geometry and bc---------------------------------------------------------------------------------
  const int ndim = DIMENSION;                         // dimension
  const int offset = _mgmesh._NoNodes[_NoLevels-1];   // mesh nodes
  const int el_sides= _mgmesh._GeomEl._n_sides[0];    // element sides
  int       el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];   // element connectivity
  int       el_neigh[NDOF_FEM];                       // bd element connectivity
  int       sur_toply[NDOF_FEMB];                     // boundary topology
  double    xx_qnds[DIMENSION*NDOF_FEM];              // element node coords
  double    xxb_qnds[DIMENSION*NDOF_FEMB];            // boundary of the element node coords


  // gauss integration  -----------------------------------------------------------------------------
  const int el_ngauss = _fe[2]->_NoGauss1[ndim-1];         // quadratic element gauss points
  const int elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];   // quadratic bd elem gauss points



  // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
  int elb_ndof[3];  elb_ndof[0]=NDOF_BK;
  elb_ndof[1]=NDOF_PB;  elb_ndof[2]=NDOF_FEMB;  // number of element boundary dofs
  int el_mat_nrows =0;                                              // number of matrix rows (dofs)
  for(int ideg=0; ideg<3; ideg++) el_mat_nrows +=_nvars[ideg]*_el_dof[ideg];
  int el_mat_ncols = el_mat_nrows;                    // square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector


  int probe_count=0;	      
  for(int kproc=0; kproc<_mgmesh._n_subdom; kproc++) {
    int nel_lev=0;
    for(int isubd=0; isubd<kproc; isubd++) {
      nel_lev += _mgmesh._off_el[0][isubd*_NoLevels+Level+1]-_mgmesh._off_el[0][isubd*_NoLevels+Level];//const
    }
    int ndof_lev=0;
    for(int pr=0; pr <kproc; pr++) {
      int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
      ndof_lev +=delta;
    }

    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*kproc+1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*kproc];   // stop element
    for(int iel=0; iel < (nel_e - nel_b); iel++) {

      // geometry and element  fields ------------------------------------

      _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds,kproc);   // Element Connectivity (el_conn) and coordinates (xx_qnds)

//       _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh); // Neighbors of the element
     for(int  iside=0; iside<el_sides; iside++)    {
    //get the global node number
    el_neigh[iside] = _mgmesh._el_neighbor[0][(iel+_mgmesh._off_el[0][Level+_NoLevels*kproc])*el_sides+iside];
  }
      // set element-nodes variables  bc (bc_q_dofs)
      get_el_dof_bc(Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);



      // ======================  p solver=======================================
//       double hdir=0;

//     for(int iql=1; iql<2; iql++)  {      // --- LINEAR ---
      for(int iside=0; iside< el_sides; iside++)  {
        if(el_neigh[iside] == -1) {

          for(int idof=0; idof<NDOF_FEMB; idof++) {                          //idof -> boundary node
            sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];//use map to find global node
            int idofb=sur_toply[idof];                                      //idofb -> element node
            elb_conn[idof]=el_conn[idofb];                                  //connectivity vector
            for(int idim=0; idim<DIMENSION; idim++) {
              xxb_qnds[idim*NDOF_FEMB+idof]=xx_qnds[idim*NDOF_FEM+idofb]; //get boundary coordinates
//               _data_eq[2].ub[idim*NDOF_FEM+idof]=xxb_qnds[idim*NDOF_FEMB+idof];
            } // idim
          } // idof

          double z_cord=xxb_qnds[3*NDOF_FEMB -1];
          if(z_cord<Z_LOW) {
            for(int tmp=0; tmp<NDOF_FEMB/2; tmp++) {
              _mgphys._node_glob_probe[probe_count]=elb_conn[tmp];
              probe_count++;
            } // tmp
          } // z_cord<

        } //end if side
      } //end iside
      // ======================  end of p solver=======================================


    } // end of element loop

  }
  // clean
  el_dof_indices.clear();



  // print on file --------------------------
  for(int ivar=5; ivar< 10; ivar++)  _mgmesh._mesh2_data[ivar]=0.;
//   printf(" NSC probe_count %d  \n",probe_count);
//   for(int   jvar=0; jvar< probe_count; jvar++)    {
//     printf(" %d  ",  _mgphys._node_glob_probe[jvar]);
//     if(jvar%20==19)  printf(" \n");
//   }
  printf("  \n");
  _mgphys._flag_probe=probe_count;
#ifdef PRINT_INFO
  std::cout << "  MGSolNSC::get_nod_probe: node defined  " << Level << " (Level) with "  << _mgphys._flag_probe  << " (nodes) "  << std::endl;
#endif


  return;
}

// ========================================================================================
/// This function compute average surface  values on  over mesh[0] and
///  transfer data to mesh[1]
void MGSolNSC::print_ext_data(
  double mesh2mesh_data[] // transfer vector from mesh[1]
) { // =====================================================================================
  int flag_probe=_mgphys._flag_probe;
  const int Level=_NoLevels-1;
  const int offset= _mgmesh._NoNodes[_NoLevels-1];

  // velocity (outlet from primary loop) data
  double sum_u=0.;
  for(int inode=0; inode< flag_probe; inode++) {
    sum_u  += (*x_old[Level])(_node_dof[Level][_mgphys._node_glob_probe[inode]])*_refvalue[0];
  }
  sum_u /=flag_probe;
  _mgmesh._mesh2_data[5+0]=0.;_mgmesh._mesh2_data[5+1]=0;// writing on mesh[1]
  _mgmesh._mesh2_data[5+2]=sum_u;// writing on mesh[1]
  mesh2mesh_data[5+0]=0.;mesh2mesh_data[5+1]=0.;     // writing on mesh[0]
  mesh2mesh_data[5+2]=sum_u;     // writing on mesh[0]

  // pressure (outlet from reactor) data
  double press=0.;
  for(int inode=0; inode< flag_probe; inode++) {
      press  += (*x_old[Level])(_node_dof[Level][_mgphys._node_glob_probe[inode]+1*offset]);
//          double pressure_out=((*x_old[_NoLevels-1])(_node_dof[_NoLevels-1][el_conn[6]+(1)*offset]));
  }
  press /=flag_probe;
  // writing data in the system vector
  _mgmesh._mesh2_data[5+DIMENSION]=press; // writing on mesh[1]
  mesh2mesh_data[5+DIMENSION]=press;      // writing on mesh[0]


#ifdef PRINT_INFO
  std::cout << " probe NSC computed  " << Level << " (Level) with "  << flag_probe  << " (nodes) "  << std::endl;
  
   std::cout << " 0n NSC probe_count " << flag_probe << std::endl;
  std::cout << " vel = " << mesh2mesh_data[5+0] << "  " << 
        mesh2mesh_data[1] << "  "<< mesh2mesh_data[5+2] << "  "<<std::endl;
  std::cout << " pressure = " <<  mesh2mesh_data[5+DIMENSION]  << std::endl;
  
#endif
 
  return;
}



/// ======================================================
/// This function controls the time step operations:
/// ======================================================
void MGSolNSC::MGTimeStep(const double time, const int) {

  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;
  if(_mgphys._flag_probe==0) {
    get_nod_probe();
  }

  /// [a] Assemblying of the rhs and matrix at the top level with GenMatRhs(time,top_level,1)
#if  PRINT_TIME==1
  std::clock_t start_time=std::clock();
#endif
  GenMatRhs(time,_NoLevels-1,1);

  /// [b] Assemblying of the other matrices with GenMatRhs(time,level,0) for all levels
  for(int Level = 0 ; Level < _NoLevels-1; Level++) GenMatRhs(time,Level,0);

#if    PRINT_TIME==1
  std::clock_t end_time=std::clock();
  std::cout << " Assembly time ="<< double(end_time- start_time) / CLOCKS_PER_SEC
            << " s "<< std::endl;
#endif
  /// [c] Solution of the linear system (MGSolverBase::MGSolve).
  MGSolve(1.e-6,40);
#if    PRINT_TIME==1
  std::clock_t end_timef=std::clock();
  std::cout << " Assembly+solution time ="<< double(end_timef- start_time) / CLOCKS_PER_SEC
            << "s "<< std::endl;
#endif
  /// [d] Update of the old solution at the top Level

  x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);
  return;
}




#endif  //ENDIF NS_EQUATIONS

