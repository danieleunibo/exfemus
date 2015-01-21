#include "Equations_conf.h"

// ============================================
#ifdef T_ADJ_EQUATIONS // 3D-2D Energy equation
// ============================================


// class local configuration -------
#include "MGSolverT_ADJ.h"
#include "MGSclass_conf.h"

// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include -----------------------------------
#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "EquationSystemsExtendedM.h"  // Equation map class

// standard lib -----------------------------------------------
#include <string.h>          // string library

// local alg lib -----------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers

#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "InterfaceFunctionM.h"
#endif

// ======================================================
/// This function constructs the 3d-2D MGSolT class
MGSolT_ADJ::MGSolT_ADJ(MGEquationsSystem& mg_equations_map_in,
               int             nvars_in[],         // KLQ number of variables
               std::string eqname_in,
               std::string varname_in):
    MGSolDA(mg_equations_map_in,nvars_in,eqname_in,varname_in),
    /// A) reading parameters
    // mesh params ------------
    _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes
    // phys params ------------
    _dt(_mgutils.get_par("dt")),       // parameter  dt
    _uref(mg_equations_map_in.get_par("Uref")),    // parameter  vel reference
    _lref(mg_equations_map_in.get_par("Lref")),    // parameter  length reference
    _Tref(mg_equations_map_in.get_par("Tref")),    // parameter  temperature reference
    _rhof(mg_equations_map_in.get_par("rho0")),    // parameter  density reference
    _muf(mg_equations_map_in.get_par("mu0")),       // parameter  viscosity reference
    _cp0(mg_equations_map_in.get_par("cp0")),      // parameter  Cp reference
    _kappa0(mg_equations_map_in.get_par("kappa0"))// parameter  conductivity reference
//   _h_conv(mg_equations_map_in.get_par("hconv")),  // parameter  heat transfer convection coefficient
//   _T_inf(mg_equations_map_in.get_par("T_inf"))
{
    //  =================================================
    /// B) setting class ariables
    _var_names[0]=varname_in;
    _refvalue[0]=_Tref;

    /// C ) setting solver type
    for (int l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVERT);

    /// D) setting nondimensional parameters
    _alpha=_kappa0/(_rhof*_cp0);
    _IPrdl=_rhof*_alpha/_muf;
    _IRe=_muf/(_rhof*_uref*_lref);
    _IPrdl_turb=1./PRT;
    _alpha_turb=0.;
//   _Nusselt=(_h_conv*_lref)/_kappa0;
    _qheat=mg_equations_map_in.get_par("qheat")*_lref/(_rhof*_cp0*_Tref*_uref);
    _qs=mg_equations_map_in.get_par("qs")/(_rhof*_cp0*_Tref*_uref);

    return;
}



//  ====================================================
/// This function assembles the matrix and the rhs:
//  ====================================================
void  MGSolT_ADJ::GenMatRhs(
    const double /* time*/, // time  <-
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
    double cont_zone=1.;
    
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
    for (int ideg=1; ideg<3; ideg++) {
        el_ndof[ideg]=_fe[ideg]->_NoShape[ndim-1];
        elb_ndof[ideg]=_fe[ideg]->_NoShape[ ndim-2];
        el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
    };
    int el_mat_ncols = el_mat_nrows;                    // square matrix
    std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector

    // coupling  fields -------------------------------------------------------------------------------
    int ns_NS=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NS_F]];  // Velocity
    int ns_T=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F]];    // Temperature
    int ns_TADJ=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[KTT_F]];    // Temperature adjoint
    double  val_tbg[3];                                        // turbulence
    double vel_g[DIMENSION];                                   // velocity
    for (int idim=0; idim<DIMENSION; idim++) vel_g[idim] =0.;  // velocity not coupled
    double Pe_h[3],f_upwind[3],h_eff[3];                       // local Peclet, upwind, h_eff
    double kappa_mg[2];
    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) ---------------------------
    A[Level]->zero();
    if (mode ==1) b[Level]->zero();  // global matrix+rhs
    DenseMatrixM KeM;
    DenseVectorM FeM;                // local  matrix+rhs
    KeM.resize(el_mat_nrows,el_mat_ncols);
    FeM.resize(el_mat_nrows);        // resize  local  matrix+rhs

    int ndof_lev=0;
        for (int pr=0;pr <_mgmesh._iproc; pr++) {
            int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
            ndof_lev +=delta;
        }
    
    
    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
    for (int iel=0; iel < (nel_e - nel_b); iel++) {
        // set to zero matrix and rhs and center
        KeM.zero();
        FeM.zero();

        // geometry and element  fields ------------------------------------
        // Element Connectivity (el_conn)  and coordinates (xx_qnds)
        _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);
        _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);

        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc(Level,iel+ndof_lev,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

        // grid and gaussian points
        for (int idim=0; idim<DIMENSION; idim++) {
            for (int d=0; d< NDOF_FEM; d++) _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];// element nodes xxg (DIM)
            // element grid distance
            h_eff[idim]=fabs(xx_qnds[idim*NDOF_FEM+2+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+0]);
            double h_1=fabs(xx_qnds[idim*NDOF_FEM+3+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+1]);
            if (h_eff[idim] <  h_1) h_eff[idim]=h_1; // Max dx diagonal term
            f_upwind[idim]=0.;                       // upwind
        }
        // element field values
        for (int deg=0; deg<3; deg++) {
            for (int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
                _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],
                                                     el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub);
            }
        }
        //  external node quantities -------------------------------------

        //  external cell properties -------------------------------------
        double phase=1.;
        double Qad =_qheat;
        double rhocp=1.;
        double rho=1.;
        double qn =_qs;
//   	_dt =0.5;
// 	vel_g[0]=4.*(xx_qnds[8]-xx_qnds[8]*xx_qnds[8]); vel_g[1]=0.; vel_g[2]=0.;
// 	vel_g[0]=-xx_qnds[8]*xx_qnds[8]+2.*xx_qnds[8]; vel_g[1]=0.; vel_g[2]=0.;  //bottom boundary layer
// 	vel_g[0]=-xx_qnds[8]*xx_qnds[8]+1.; vel_g[1]=0.; vel_g[2]=0.;             //top boundary layer
	
// 	vel_g[0]=1.; vel_g[1]=0.; vel_g[2]=0.;  
// 	cont_zone=1.;
// 	cont_zone=(xx_qnds[17]-xx_qnds[17]*xx_qnds[17])*
// 		  (8.*xx_qnds[8]-xx_qnds[8]*xx_qnds[8])/4.;  //for Paraview (coordsY-coordsY*coordsY)*(8*coordsX-coordsX*coordsX)/4
// 	if(xx_qnds[8]>0.1 && xx_qnds[8]<0.9 && xx_qnds[17]>0.3 && xx_qnds[17]<0.7) cont_zone=1;

	// ======================================================================
        // Volume =============================================================
        // ======================================================================

        // ---------------------------------------------
        /// c) gaussian integration loop (n_gauss)
        // --------------------------------------------
        for (int qp=0; qp< el_ngauss; qp++) {

            // shape functions at gaussian points -----------------------------------
            for (int ideg=1; ideg<3; ideg++) { // linear-quadratic  [1]=linear [2]=quad
                det[ideg]      = _fe[ideg]->Jac(qp,xx_qnds,InvJac[ideg]);     // Jacobian
                JxW_g[ideg] =det[ideg]*_fe[ideg]->_weight1[ndim-1][qp];       // weight
                _fe[ideg]->get_phi_gl_g(ndim,qp,_phi_g[ideg]);               // shape funct
                _fe[ideg]->get_dphi_gl_g(ndim,qp,InvJac[ideg],_dphi_g[ideg]); // global coord deriv
            }

            double alpha_eff = 1.;//_IPrdl*_IRe;

            //  fields -----------------------------------------------------------
            interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],
                          _phi_g[2],el_ndof[2],_ub_g[2]); // quadratic
#ifdef AXISYM   // axisymmetric (index -> 0)
            JxW_g[2]  *=_ub_g[2][0];
#endif
	    
// 	    cont_zone = (_ub_g[2][1]-_ub_g[2][1]*_ub_g[2][1])*
// 		  (8.*_ub_g[2][0]-_ub_g[2][0]*_ub_g[2][0])/4.;
            vel_g[0]=-_ub_g[2][1]*_ub_g[2][1]+2.*_ub_g[2][1];
//             if(_ub_g[2][ns_TADJ] < 1.e-8) _ub_g[2][ns_TADJ]=1.e-8;
            /// d) Local (element) assemblying energy equation
            // *********************** *******************************
            for (int i=0; i<el_ndof[2]; i++)     {
                // set up row i
	        double Adv=0.;
                const double phii_g=_phi_g[2][i];
                for (int idim=0; idim< ndim; idim++) {  
		  dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
// 		  Adv +=_ub_g[2][ns_NS+idim]*dphiidx_g[2][idim];
		  Adv +=vel_g[idim]*dphiidx_g[2][idim];
		}

                double dtxJxW_g=_dt*JxW_g[2]*_bc_vol[i];


                // Rhs Assemblying  ----------------------------------------------------------
                if (mode == 1) {
                    // rhs
                    FeM(i) += dtxJxW_g*(
                                  rhocp*_ub_g[2][ns_TADJ]*phii_g/_dt +  // time
                                  (_ub_g[2][ns_T]-T_CONTROL)*cont_zone*phii_g                       // source
                              );
                }


                // Matrix Assemblying ---------------------------
                for (int j=0; j<el_ndof[2]; j++) {
                    double phij_g= _phi_g[2][j];
                    double Lap=0.;
                    for (int idim=0; idim< ndim; idim++) {
                        dphijdx_g[2][idim]=_dphi_g[2][j+idim*el_ndof[2]];
                        Lap +=(alpha_eff) // upwind
                              *dphijdx_g[2][idim]*dphiidx_g[2][idim];            // Laplacian
//             Grad += dphijdx_g[2][idim]*phii_g;            // Gradient
                    }
                    // energy-equation
                    KeM(i,j) +=dtxJxW_g*(
                                   rhocp*phii_g*phij_g/_dt // time term
                                   + rhocp*Adv*phij_g      //adv
                                   + Lap                   //diff
                               );
                    // anisotropic term
//        for (int idim=0; idim< ndim; idim++)
//    for (int jdim=0; jdim< ndim; jdim++)
//            KeM(i,j) +=dtxJxW_g*ANYS_TTB*( // anisotropic term
//                        -1.*f_anis*_IRe*_nut_ratio*
//                        (vel_dxg[jdim+idim*DIMENSION]+vel_dxg[idim+jdim*DIMENSION])
//                      )*dphijdx_g[jdim]*dphiidx_g[idim];


                }
            } // ----------------------------------------
        } // end of the quadrature point qp-loop ***********************


        // ======================================================================
        // ====================== boundary ======================================
        // ======================================================================

        for (int iside=0; iside< el_sides; iside++)  {
            if (el_neigh[iside] == -1) {

                double alpha_eff = 1.;//_IPrdl*_IRe + _nut_ratio*_IRe/_IPrdl_turb;

                for (int idof=0; idof<NDOF_FEMB; idof++) {
                    sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
                    int idofb=sur_toply[idof];
                    elb_conn[idof]=el_conn[idofb];                 // connectivity vector
                    for (int idim=0; idim<DIMENSION; idim++) {
                        xxb_qnds[idim*NDOF_FEMB+idof]=xx_qnds[idim*NDOF_FEM+idofb]; // coordinates
                        _data_eq[2].ub[idim*NDOF_FEM+idof]=xxb_qnds[idim*NDOF_FEMB+idof];
                    }
                }

                // Dirichlet boundary conditions  ***********************************
//         printf(" \n  bc vol %d  bd %d ",_bc_vol[sur_toply[NDOF_FEMB-1]],_bc_bd[sur_toply[NDOF_FEMB-1]]);
                if (_bc_vol[sur_toply[NDOF_FEMB-1]] ==0) {

                    //[NDOF_FEMB-1] is the midpoint of a quadratic FEM element (HEX9 (2D) or HEX27 (3D))
                    int bc_s=(int)_bc_bd[sur_toply[NDOF_FEMB-1]];     // b cond
                    double det= _fe[2]->JacSur(elb_ngauss-1,xxb_qnds,InvJac[2]);// jacobian
                    double Ipenalty=det/_dt;                               // Dirichlet bc flag
                    // local boundary loop   ---------------------------------------
                    for (int lb_node=0; lb_node< elb_ndof[2]; lb_node++) {
                        int lv_node= sur_toply[lb_node]; // local vol index
                        // flag setup (\int bc_var*T+bc_var*val)
                        int  bc_val = (int)((bc_s&2)>>1);  // (1?)non homogeneous
                        int  bc_var = (int)(bc_s%2);       // (?1) variable
                        // Assemblying  Matrix & rhs
                        if (mode == 1) {

//                             FeM(lv_node) += bc_val*Ipenalty*_data_eq[2].ub[ns_T*NDOF_FEM+lv_node];
//	                       FeM(lv_node) += (1- bc_val)*(1-bc_var)*Ipenalty*700.;
//                             FeM(lv_node) +=(1- bc_val)*(1-bc_var)*Ipenalty*(573.15*0.+_ub_g[2][1]*0.);


                        }
                        KeM(lv_node,lv_node) += Ipenalty;  //  Dirichlet bc
                    }// lb_node -end  local boundary loop -------------------------
                } // end if Dirichlet  boundary conditions
                // **********************************************************************

//                else if (_bc_bd[sur_toply[NDOF_FEMB-1]] != 0) {
//                     // Non homogenous Neumann boundary conditions  ***********************************
// 
//                     // gaussian integration loop (n_gauss)
//                     // -----------------------------------------------
//                     for (int qp=0; qp<  elb_ngauss; qp++) {
// 
//                         // quad/linear  [2]=quad [1]=linear------------------------------------
//                         det[2]  = _fe[2]->JacSur(qp,xxb_qnds,InvJac[2]);    // local coord _phi_g and jac
//                         JxW_g[2]=det[2]*_fe[2]->_weight1[ndim-2][qp];// weight
//                         _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);   // global coord _phi_g
// 
// 
// #ifdef AXISYM   // axisymmetric  (index ->0)
//                         interp_el_sol(_data_eq[2].ub,0,DIMENSION,_phi_g[2],elb_ndof[2],_ub_g[2]);
//                         JxW_g[2]  *=_ub_g[2][0];
// 
// #endif
//                         double f_qn=1.;
// 
// // #ifdef TBK_EQUATIONS
// //             double kappa_w=kappa_mg[0]+1.e-20;
// //             double ystar=rho*0.547722557505*_y_dist*sqrt(kappa_w)/(_IRe);(ret*0.001841/(10340*0.03025))
// //             f_qn=(ystar/_IPrdl_turb< 11.66) ? 1.:2.5*_IPrdl_turb*log(9.*ystar/_IPrdl_turb)/ystar;
// //
// // #endif
// 
//                         f_qn=1.;
// 
// // 	                printf(" \n  ystar %lf kappa %lf y %lf f_qn %lf",ystar/_IPrdl_turb,kappa_w,_y_dist,f_qn);
//                         // ***********************************
//                         // local side loop (over the node face)
//                         for (int lsi_node=0; lsi_node< elb_ndof[2]; lsi_node++) {
// 
//                             // set up row i
//                             const double phii_g=_phi_g[2][lsi_node]; // boundary test function
//                             int lei_node= sur_toply[lsi_node]; // local element index
//                             int bc_s=(int)_bc_bd[lei_node];
//                             int bc_v=(int)_bc_vol[lei_node];
//                             double dtxJxW_g=JxW_g[2]*bc_v; // Neumann bc flag and Jac
// 
//                             // flag setup: +++++++++++++++++++++++++++++++++++++++
//                             //  \int (bc_var*T+bc_val*val)ds
//                             int  bc_val = (int)((bc_s&2)>>1);  // (1?) non-homogeneous
//                             int  bc_var = (int)(bc_s%2);       // (?1) tau=A*variable
// // 	       if(bc_val==1 && bc_var==0) std::cout<<_ub_g[2][0];
//                             // ++++++++++++++++++++++++++++++++++++++++++++++++++++
// 
//                             // Assemblying rhs ----------------------------
//                             if (mode == 1)
//                                 FeM(lei_node) += bc_val*dtxJxW_g*phii_g*(// non-homogeneous
//                                                      bc_var*(((_ub_g[2][1]>0.0000001)?1.:-1.)*alpha_eff*360000./(10340.*145.75*0.263659*0.002874)) // Robin bc (k*dt/dn = h*T0)
//                                                      +(1 -bc_var)*qn*f_qn                   // Neumann heat flux bc
//                                                  );
// 
//                             // Assemblying Matrix ---------------------------------
//                             for (int lsj_node=0; lsj_node< elb_ndof[2];  lsj_node++) {
//                                 KeM(lei_node,sur_toply[lsj_node]) += 0.*dtxJxW_g*bc_var*alpha_eff*phii_g*_phi_g[2][lsj_node]; // Robin bc  (k*dt/dn = h*(-T))
//                             }// end j  ---------------------------------------
// 
//                         }// i   +++++++++++++++++++++++++++++++
//                     } // end of the quadrature point qp-loop **********************
// 
// 
//                 } // Neumann non homog

            } //end if side
        } // ======================  end for boundary =======================================

        /// e) Global assemblying energy equation
        A[Level]->add_matrix(KeM,el_dof_indices);                  // global matrix
        if (mode == 1)   b[Level]->add_vector(FeM,el_dof_indices); // global rhs

    } // end of element loop
    // clean
    el_dof_indices.clear();
    A[Level]->close();
    if (mode == 1) b[Level]->close();
#ifdef PRINT_INFO
    std::cout<< " Matrix Assembled(T)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif

    return;
}


// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolT_ADJ::MGTimeStep(
    const double time,  // time
    const int /*iter*/  // Number of max inter
) {
// =========================================================================================

/// A) Set up the time step
#if PRINT_INFO==1
    std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;
#endif
    /// B) Assemblying of the Matrix-Rhs
#if PRINT_TIME==1
    std::clock_t start_time=std::clock();
#endif
    GenMatRhs(time,_NoLevels-1,1);                                                // matrix and rhs
    for (int Level = 0 ; Level < _NoLevels-1; Level++)   GenMatRhs(time,Level,0); // matrix
#if PRINT_TIME==1
    std::clock_t end_time=std::clock();
    std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif

    /// C) Solution of the linear MGsystem (MGSolT::MGSolve).
     MGSolve(1.e-8,60);
#if PRINT_TIME==1
    end_time=std::clock();
    std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif

    /// D) Update of the old solution at the top Level  (MGSolT::OldSol_update),
    x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);
    return;
}// =======================================================================================


#endif
// #endif // personal application

