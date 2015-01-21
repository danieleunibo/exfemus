#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS // 3D-2D Energy equation
// ============================================


// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverT.h"       // Navier-Stokes class header file


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
// ==================================================================

// ======================================================
/// This function constructs the 3d-2D MGSolT class
MGSolT::MGSolT(MGEquationsSystem& mg_equations_map_in,
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
void  MGSolT::GenMatRhs(
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
    int ns_T=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F]];    // Temperature
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
  	_dt =0.5;
#ifdef TBK_EQUATIONS     // turbulence field  

#ifdef DIST_FIX
        _y_dist=DIST_FIX;
#else
        _y_dist=_mgmesh._dist[ iel+nel_b];
#endif
#endif
        kappa_mg[0] =0.;
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

            double alpha_eff = _IPrdl*_IRe;

            //  fields -----------------------------------------------------------
            interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],
                          _phi_g[2],el_ndof[2],_ub_g[2]); // quadratic
#ifdef AXISYM   // axisymmetric (index -> 0)
            JxW_g[2]  *=_ub_g[2][0];
#endif

            // Temperature field -> [T_F] -> (quad, _indx_eqs[NS_F])
#ifndef CONST
            rhocp *=densityT(uold_g[0])*cipiT(uold_g[0]);// nondimensional density (rho/rhoref)
            alpha_eff *=kappaT(uold_g[0]);// nondimensional viscosity (mu/muref)
#endif

            // Velocity field -> [NS_F] -> (quad, _indx_eqs[NS_F])
            double mod2_vel=1.e-20;  // velocity modulus
            _sP=1.e-20;              // turbulence production
#ifdef NS_EQUATIONS
            const int ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NS_F]];  // Navier-Stokes

            // velocity field (->vel_g[])
            for (int idim=0; idim< DIMENSION; idim++) {
                vel_g[idim] =_ub_g[2][ns_idx+idim]; // velocity field
                mod2_vel +=vel_g[idim]*vel_g[idim]; // velocity modulus
            }
            mod2_vel =sqrt(mod2_vel); // velocity modulus

            // velocity field derivatives (->vel_dxg[][])
            double vel_dxg[DIMENSION*DIMENSION];
            interp_el_gdx(_data_eq[2].ub,ns_idx,DIMENSION,_dphi_g[2],el_ndof[2],vel_dxg);

            // turbulence production term _sP
            for (int idim=0; idim< DIMENSION; idim++) {
                for (int jdim=0; jdim< DIMENSION; jdim++) {
                   const double sss=0.5*(vel_dxg[idim+jdim*DIMENSION]+vel_dxg[jdim+idim*DIMENSION]);
                   _sP += 2.*sss*sss;
                }
                // upwind
                Pe_h[idim]=0.5*mod2_vel*h_eff[idim]/(_IRe*_IPrdl);                             // Peclet number
                f_upwind[idim]=UP_WIND_T*rhocp*0.5*(1./tanh(Pe_h[idim])-1./Pe_h[idim])*
                               h_eff[idim]/(mod2_vel);                            // upwind
            }

#ifdef TBK_EQUATIONS    // Turbulent viscosity  -> [K_F] -> (quad, _indx_eqs[K_F])  
            const int k_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[K_F]];  // kappa-equations
            _kappa_g[0]= _ub_g[2][k_idx];// val_tbg[0]=
            _kappa_g[1]= _ub_g[2][k_idx+1]; //  val_tbg[1]
            //      val_tbg[2]=y_dist;val_tbg[3]=sqrt(sP);
#ifdef TTBK_EQUATIONS    // Turbulent viscosity   -> [KTT_F] -> (quad, _indx_eqs[KTT_F])
            const int kTT_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[KTT_F]];  // kappa-equations
            _kappaT_g[0]= _ub_g[2][kTT_idx];// val_tbg[0]=
            _kappaT_g[1]= _ub_g[2][kTT_idx+1]; //  val_tbg[1]
#endif
            _nut_ratio=0.;
            f_mu(val_tbg);
            kappa_mg[0] +=  _kappa_g[0];
            alpha_eff = _IPrdl*_IRe +_nut_ratio*_IRe*_IPrdl_turb; 
	    //        std::cout << alpha_eff << "  " <<  1./_IRe << "  " << _nut_ratio*_IRe << " \n";
//                       *0.09*0.2*sqrt(sqrt(_nut_ratio/0.09))*(1.-ANYS_TTB) (1.-ANYS_TTB)+f_anis*0.666666666667*_kappa_g
#endif  // TBK_EQUATIONS  
#endif  // NS_EQUATIONS 
//       std::cout << " alfa_eff" << alpha_eff << " ni " << _nut_ratio << " up 0 " <<  f_upwind[0]*vel_g[0]*vel_g[0] << " up 1 " <<  f_upwind[1]*vel_g[1]*vel_g[1] << " h 0 " <<  h_eff[0] << " h 1 " <<  h_eff[1]<< "\n";
            /// d) Local (element) assemblying energy equation
            // *********************** *******************************
            for (int i=0; i<el_ndof[2]; i++)     {
                // set up row i
                const double phii_g=_phi_g[2][i];
                for (int idim=0; idim< ndim; idim++)  dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
                double dtxJxW_g=JxW_g[2]*_bc_vol[i];


                // Rhs Assemblying  ----------------------------------------------------------
                if (mode == 1) {

                    // rhs
                    FeM(i) += dtxJxW_g*(
                                  rhocp*_ub_g[2][ns_T]*phii_g/_dt   // time
                                  +Qad*phii_g                       // heat source
                                  +rhocp*vel_g[1]*phii_g
// //                                   *36000*2.*0.03025/(10340.0*145.75*.62*(0.03025*0.03025-0.004*0.004))
// // 		      -(alpha_eff+f_upwind[1]*vel_g[1]*vel_g[1])*dphiidx_g[2][1]
// 				  *360./(10340.*145.75*V_MID*0.03025)  //plane
				  *2.*360./(10340.*145.75*V_MID*0.03025)  //cylinder
                              );
                }


                // Matrix Assemblying ---------------------------
                for (int j=0; j<el_ndof[2]; j++) {
                    double phij_g= _phi_g[2][j];
                    double Adv=0.;
                    double Lap=0.;
                    for (int idim=0; idim< ndim; idim++) {
                        dphijdx_g[2][idim]=_dphi_g[2][j+idim*el_ndof[2]];
                        Adv +=vel_g[idim]*dphijdx_g[2][idim];                    // advection
                        Lap +=(alpha_eff+f_upwind[idim]*vel_g[idim]*vel_g[idim]) // upwind
                              *dphijdx_g[2][idim]*dphiidx_g[2][idim];            // Laplacian
//             Grad += dphijdx_g[2][idim]*phii_g;            // Gradient
                    }
                    // energy-equation
                    KeM(i,j) +=dtxJxW_g*(
                                   rhocp*phii_g*phij_g/_dt // time term
                                   + rhocp*Adv*phii_g      //adv
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

                double alpha_eff = _IPrdl*_IRe + _nut_ratio*_IRe/_IPrdl_turb;

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

                            FeM(lv_node) += bc_val*Ipenalty*_data_eq[2].ub[ns_T*NDOF_FEM+lv_node];
//	      FeM(lv_node) += (1- bc_val)*(1-bc_var)*Ipenalty*700.;
                            FeM(lv_node) +=(1- bc_val)*(1-bc_var)*Ipenalty*(573.15*0.+_ub_g[2][1]*0.);


                        }
                        KeM(lv_node,lv_node) += Ipenalty;  //  Dirichlet bc
                    }// lb_node -end  local boundary loop -------------------------
                } // end if Dirichlet  boundary conditions
                // **********************************************************************

                else if (_bc_bd[sur_toply[NDOF_FEMB-1]] !=0) {

                    // Non homogenous Neumann boundary conditions  ***********************************

                    // gaussian integration loop (n_gauss)
                    // -----------------------------------------------
                    for (int qp=0; qp<  elb_ngauss; qp++) {

                        // quad/linear  [2]=quad [1]=linear------------------------------------
                        det[2]  = _fe[2]->JacSur(qp,xxb_qnds,InvJac[2]);    // local coord _phi_g and jac
                        JxW_g[2]=det[2]*_fe[2]->_weight1[ndim-2][qp];// weight
                        _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);   // global coord _phi_g


#ifdef AXISYM   // axisymmetric  (index ->0)
                        interp_el_sol(_data_eq[2].ub,0,DIMENSION,_phi_g[2],elb_ndof[2],_ub_g[2]);
                        JxW_g[2]  *=_ub_g[2][0];

#endif
                        double f_qn=1.;

// #ifdef TBK_EQUATIONS
//             double kappa_w=kappa_mg[0]+1.e-20;
//             double ystar=rho*0.547722557505*_y_dist*sqrt(kappa_w)/(_IRe);(ret*0.001841/(10340*0.03025))
//             f_qn=(ystar/_IPrdl_turb< 11.66) ? 1.:2.5*_IPrdl_turb*log(9.*ystar/_IPrdl_turb)/ystar;
//
// #endif

                        f_qn=1.;

// 	                printf(" \n  ystar %lf kappa %lf y %lf f_qn %lf",ystar/_IPrdl_turb,kappa_w,_y_dist,f_qn);
                        // ***********************************
                        // local side loop (over the node face)
                        for (int lsi_node=0; lsi_node< elb_ndof[2]; lsi_node++) {

                            // set up row i
                            const double phii_g=_phi_g[2][lsi_node]; // boundary test function
                            int lei_node= sur_toply[lsi_node]; // local element index
                            int bc_s=(int)_bc_bd[lei_node];
                            int bc_v=(int)_bc_vol[lei_node];
                            double dtxJxW_g=JxW_g[2]*bc_v; // Neumann bc flag and Jac

                            // flag setup: +++++++++++++++++++++++++++++++++++++++
                            //  \int (bc_var*T+bc_val*val)ds
                            int  bc_val = (int)((bc_s&2)>>1);  // (1?) non-homogeneous
                            int  bc_var = (int)(bc_s%2);       // (?1) tau=A*variable
// 	       if(bc_val==1 && bc_var==0) std::cout<<_ub_g[2][0];
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++

                            // Assemblying rhs ----------------------------
                            if (mode == 1)
                                FeM(lei_node) += bc_val*dtxJxW_g*phii_g*(// non-homogeneous
                                                     bc_var*(((_ub_g[2][1]>0.0000001)?1.:-1.)*alpha_eff*360000./(10340.*145.75*0.263659*0.002874)) // Robin bc (k*dt/dn = h*T0)
                                                     +(1 -bc_var)*qn*f_qn                   // Neumann heat flux bc
                                                 );

                            // Assemblying Matrix ---------------------------------
                            for (int lsj_node=0; lsj_node< elb_ndof[2];  lsj_node++) {
                                KeM(lei_node,sur_toply[lsj_node]) += 0.*dtxJxW_g*bc_var*alpha_eff*phii_g*_phi_g[2][lsj_node]; // Robin bc  (k*dt/dn = h*(-T))
                            }// end j  ---------------------------------------

                        }// i   +++++++++++++++++++++++++++++++
                    } // end of the quadrature point qp-loop **********************


                } // Neumann non homog

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
    for (int Level = 0 ; Level < _NoLevels-1; Level++)   GenMatRhs(time,Level,0); // matrix
#if PRINT_TIME==1
    std::clock_t end_time=std::clock();
    std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif

    /// C) Solution of the linear MGsystem (MGSolT::MGSolve).
     MGSolve(1.e-6,40);
#if PRINT_TIME==1
    end_time=std::clock();
    std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif

    /// D) Update of the old solution at the top Level  (MGSolT::OldSol_update),
    x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);
    return;
}// =======================================================================================







// =====================================================
#ifdef TBK_EQUATIONS
// =====================================================
void  MGSolT::f_mu(double val[]) {

    if (_kappa_g[0]< 1.e-20) _kappa_g[0]= 1.e-20; // kappa
    if (_kappa_g[1]< 1.e-20) _kappa_g[1]= 1.e-20; // kappa
    double tau_k=1.;  // turbulent time

#if (TBK_EQUATIONS==0)   //kappa  (Prandtl length)--------------------------
    tau_k=_y_dist/sqrt(_kappa_g[0]);
#endif  // end kappa -------------------------------------------------------

#if (TBK_EQUATIONS/2==1)  // kappa-epsilon (TBK_EQUATIONS=2 or 3) ---------- 
    tau_k= CMU*_kappa_g[0]/_kappa_g[1];// tau=1/omega=CMU*kappa/epsilon
#endif   // kappa-epsilon --------------------------------------------------

#if (TBK_EQUATIONS/2==2)       // kappa-omega ------------------------------
    tau_k=1./_kappa_g[1]; // tau=1/omega

#endif  // -----------------------------------------------------------------
    if ( tau_k> MAX_TAU)  tau_k= MAX_TAU;
// turbulent viscosity
    double Re_t= _kappa_g[0]*tau_k/_IRe;
    if (Re_t > MU_TOP)  Re_t =MU_TOP;
    if (Re_t < MU_LOW)  Re_t =MU_LOW;
    _nut_ratio=Re_t;

//     // Boundary corrections
   double R_t=Re_t/CMU;                                             // turbulent Reynolds number
   double R_eps = _y_dist*sqrt(_kappa_g[0]/sqrt(R_t))/_IRe;
//    double R_eps =_y_dist/sqrt((_muf/_rhof)*sqrt((_muf/_rhof)/_kappa_g[1]));    // *sqrt(_kappa_g[0]/sqrt(R_t))/_IRe; 
//    double y_plus = _y_dist*0.547722557505*sqrt(_kappa_g[0])/_IRe;
   
    
    
#ifdef LOWRE     /// B) Low-Reynolds  correction model
      double f_mu =(1.-exp(-1.*R_eps/(14.)))*(1.-exp(-1.*R_eps/14.));
  _nut_ratio *= f_mu * (1.+5./pow(R_t,0.75)*exp(-1.*R_t*R_t/40000.));    // correction 
#endif
    
#ifdef SST       ///C) SST k-w model
    
    // F1 and F2 coeff calculation
      double F1,F2;
      double F1_first  = 500.*_IRe*tau_k/(_y_dist*_y_dist);
      double F2_first= 2.*sqrt(_kappa_g[0])*tau_k/(BETASTAR*_y_dist);
      if (F1_first > F2_first) F2_first=F1_first;
      F2=tanh(F2_first*F2_first);

      // nu_t calculation
//       double alpha_star        = 1.;//(0.024+ Re_t/6.)/(1.+ Re_t/6.);
      double alpha_starb = sqrt(_sP)*F2*tau_k/0.31; 
      if (alpha_starb < 1.) {alpha_starb=1.;} /*printf("alpha star is %f and and %f and %f\n",alpha_starb, 1./tau_k, F2);*/
      _nut_ratio /= alpha_starb;
 
      if (_nut_ratio > MU_TOP)  _nut_ratio =MU_TOP;
      if (_nut_ratio < MU_LOW)  _nut_ratio =MU_LOW;
//       printf("mu_turb is %f \n",_mu_turb);
      
#endif
    
    

#ifdef TTBK_EQUATIONS
    /// A) Energy  Turbulent viscosity
    // Energy  Turbulent viscosity
//     if(_kappaT_g[1]< 1.e-20) _kappaT_g[1]= 1.e-20;            // kappa
//     if(_kappaT_g[0]< 1.e-20) _kappaT_g[0]= 1.e-20;           // kappa

#if (TTBK_EQUATIONS/2==1)  // kappa-epsilon (TBK_EQUATIONS=2 or 3) ---------- 
 double tauT_k= CMU*_kappaT_g[0]/_kappaT_g[1]; // tau=1/omega=CMU*kappa/epsilon
  #endif   // kappa-epsilon --------------------------------------------------

  #if (TTBK_EQUATIONS/2==2)       // kappa-omega ------------------------------
 double tauT_k=1./_kappaT_g[1]; // tau=1/omega
  #endif
  if (tauT_k> MAX_TAU) tauT_k= MAX_TAU; // tau=1/omega=CMU*kappa/epsilon
    double rT=tauT_k/tau_k; 

  #if (TTBK_EQUATIONS/2==1)     // kappa_theta-epsilon_theta
      double f_alpha =(1.-exp(-1.*R_eps/(19.*sqrt(_IPrdl)) ))*(1.-exp(-1.*R_eps/14.));
      double f_d=exp(-1.*R_t*R_t/(200.*200.));
      double f_asym =(1.-exp(-1.*R_eps/(19.*sqrt(_IPrdl)) ))*(1.-exp(-1.*R_eps/14.));

       double a_wall = sqrt(2.*rT)*1.3*_IPrdl/pow(R_t,0.75)*f_d*f_alpha;
       double a_inter =2.*rT/(.3 +rT)*f_alpha*exp(-1.*R_t*R_t/(500.*500.));
       double asymp = 0.9*f_asym;

      double b_nu = f_mu*(1.+5./pow(R_t,0.75)*exp(-1.*R_t*R_t/40000.));
      
      _IPrdl_turb=1.11111*(a_inter+a_wall+asymp)/b_nu;
//        _IPrdl_turb=1./(0.85+0.7*_IPrdl/(_nut_ratio)); // kays model
//       _IPrdl_turb=1./2.3;                     // SED model
//      std::cout <<  _y_dist<<  " " << 1./_IPrdl_turb << "\n ";
  #endif

  #if (TTBK_EQUATIONS/2==2)       // kappa_theta-omega_theta ------------------------------
      double f_alpha =(1.-exp(-1.*R_eps/(15.*sqrt(_IPrdl)) ))*(1.-exp(-1.*R_eps/16.));
      double f_d=exp(-1.*R_t*R_t/(200.*200.));

       double a_wall = sqrt(2.*rT)*5.*_IPrdl/pow(R_t,0.75)*f_d;
       double a_inter =2.*rT/(.3+rT)*exp(-1.*R_t*R_t/(500.*500.));
       double asymp = 0.8;

//       double b_nu = f_mu*(1.+5./pow(R_t,0.75)*exp(-1.*R_t*R_t/40000.));
      
//       _IPrdl_turb=1.11111111111*(f_alpha*(a_inter+a_wall+asymp));
      
//       _IPrdl_turb=1./2.4;
      _IPrdl_turb=1./(0.85+0.7*_IPrdl/(_nut_ratio)); //Kays correlation
//     std::cout <<  _y_dist<<  " " << 1./_IPrdl_turb << "\n ";
  #endif

#endif
    return;
}
// end ---   f_mu -------------------------
#endif // ----------  end TBK_EQUATIONS  
#endif
// #endif // personal application

