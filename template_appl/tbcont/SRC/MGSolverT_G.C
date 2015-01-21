#include "Equations_conf.h"

// ============================================
#ifdef T_G_EQUATIONS // 3D-2D Energy equation
// ============================================

// class local configuration -------
#include "MGSolverT_G.h"
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
MGSolT_G::MGSolT_G(MGEquationsSystem& mg_equations_map_in,
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
    _beta = mg_equations_map_in.get_par("beta");
    _lambda = mg_equations_map_in.get_par("lambda");
    return;
}



//  ====================================================
/// This function assembles the matrix and the rhs:
//  ====================================================
void  MGSolT_G::GenMatRhs(
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
    double vel_g[DIMENSION];                                   // velocity
    for (int idim=0; idim<DIMENSION; idim++) vel_g[idim] =0.;  // velocity not coupled
    double     normal[DIMENSION];
    
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
        elb_ndof[ideg]=_fe[ideg]->_NoShape[ndim-2];
        el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
    };
    int el_mat_ncols = el_mat_nrows;                    // square matrix
    std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector

    // coupling  fields -------------------------------------------------------------------------------
    int ns_NS=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NS_F]];  // Velocity
    int ns_T=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F]];    // Temperature
    int ns_TADJ=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[KTT_F]];    // Temperature adjoint

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
        for (int idim=0; idim<DIMENSION; idim++)
            for (int d=0; d< NDOF_FEM; d++) 
	      _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];// element nodes xxg (DIM)

        // element field values
        for (int deg=0; deg<3; deg++) {
            for (int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
                _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],
                                                     el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub);
            }
        }

        double cont_zone=1.;
//         cont_zone=(xx_qnds[17]-xx_qnds[17]*xx_qnds[17])*
// 		  (8.*xx_qnds[8]-xx_qnds[8]*xx_qnds[8])/4.;  //for Paraview (coordsY-coordsY*coordsY)*(8*coordsX-coordsX*coordsX)/4
							     // (T-0.7)*(T-0.7)*cont
	// ======================================================================
        // Volume =============================================================
        // ======================================================================

        // volume gaussian integration loop 
	for (int qp=0; qp< el_ngauss; qp++) {
	      det[2] = _fe[2]->Jac(qp,xx_qnds,InvJac[2]);     // Jacobian
	      JxW_g[2] =det[2]*_fe[2]->_weight1[ndim-1][qp];       // weight
	      _fe[2]->get_phi_gl_g(ndim,qp,_phi_g[2]);               // shape funct
	      
		for (int i=0; i<el_ndof[2]; i++) {
		double dtxJxW_g= _dt*JxW_g[2]*_bc_vol[i];
                // Rhs Assemblying  -----------------------------
                if (mode == 1) FeM(i) += T_CONTROL*dtxJxW_g*_phi_g[2][i];

                // Matrix Assemblying ---------------------------
                 KeM(i,i) += dtxJxW_g*_phi_g[2][i];
		} // ----------------------------------------
	 } // end volume integration loop


        // ======================================================================
        // ====================== boundary ======================================
        // ======================================================================
	      double theta_dx[NDOF_FEM*(DIMENSION)];
	      for(int init=0; init< NDOF_FEM*DIMENSION; init++) theta_dx[init]=0.;
                
		for(int qp=0; qp< NDOF_FEM; qp++)  {
		  det[2] = _fe[2]->Jac(qp,xx_qnds,InvJac[2]);     // Jacobian
                  _fe[2]->get_dphi_gl_g(ndim,qp,InvJac[2],_dphi_g[2]); // global coord deriv
		  for (int idim=0; idim<DIMENSION; idim++) {
		    for(int jnode=0; jnode< NDOF_FEM; jnode++)  {
		    theta_dx[qp+idim*NDOF_FEM]+=_data_eq[2].ub[ns_TADJ*NDOF_FEM+jnode]*_dphi_g[2][jnode+idim*NDOF_FEM];
		    }
		  }
		}
        
        
        
        for (int iside=0; iside< el_sides; iside++)  {
            if (el_neigh[iside] == -1) {
	      
	      double Ipenalty=0.;

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

                if (_bc_bd[sur_toply[NDOF_FEMB-1]] == 2) { // Begin Equation defined on the boundary!!!  ***********************************

		    // volume gaussian integration loop 
		    for (int qp=0; qp< el_ngauss; qp++) {

			// shape functions at gaussian points -----------------------------------
                            det[2]      = _fe[2]->Jac(qp,xx_qnds,InvJac[2]);     // Jacobian
			    JxW_g[2] =det[2]*_fe[2]->_weight1[ndim-1][qp];       // weight
			    _fe[2]->get_phi_gl_g(ndim,qp,_phi_g[2]);               // shape funct
			    _fe[2]->get_dphi_gl_g(ndim,qp,InvJac[2],_dphi_g[2]); // global coord deriv
                            Ipenalty=det[2]/_dt; 
     
			 //  fields -----------------------------------------------------------
			  interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],
					_phi_g[2],el_ndof[2],_ub_g[2]); // quadratic

// 			  cont_zone = (_ub_g[2][1]-_ub_g[2][1]*_ub_g[2][1])*
// 		         (8.*_ub_g[2][0]-_ub_g[2][0]*_ub_g[2][0])/4.;
			 
			  vel_g[0]=-_ub_g[2][1]*_ub_g[2][1]+2.*_ub_g[2][1];
			 /// d) Local (element) assemblying equation
			   // *********************** *******************************
			  for (int i=0; i<el_ndof[2]; i++)     {
			  // set up row i
			  const double phii_g=_phi_g[2][i];
			  double Adv_theta=0.;
			  double Lap=0.;
			  for (int idim=0; idim< ndim; idim++) {  
			      dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
// 			      Adv_theta += _ub_g[2][ns_NS+idim]*dphiidx_g[2][idim];
			      Adv_theta += vel_g[idim]*dphiidx_g[2][idim];
			      Lap += alpha_eff*theta_dx[qp+idim*NDOF_FEM]*dphiidx_g[2][idim];
			  }
			  
			  int  bc_val = (int)((_bc_bd[i]&2)>>1);  // (1?) non-homogeneous
			  double dtxJxW_g= _dt*JxW_g[2]*bc_val;


			  if (mode == 1) FeM(i) += dtxJxW_g*(cont_zone*(T_CONTROL-_ub_g[2][ns_T])*phii_g
			                           + Adv_theta*_ub_g[2][ns_TADJ]
			                           + Lap
						   );
			  }
		    }//end volume integration

                    // surface gaussian integration loop
                    // -----------------------------------------------
                    for (int qp=0; qp<  elb_ngauss; qp++) {

                        // quad/linear  [2]=quad [1]=linear------------------------------------
                        det[2]  = _fe[2]->JacSur(qp,xxb_qnds,InvJac[2]);    // local coord _phi_g and jac
                        JxW_g[2]=det[2]*_fe[2]->_weight1[ndim-2][qp];// weight
                        _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);   // global coord _phi_g
// 			if (xxb_qnds[1]-xxb_qnds[0] > 1.e-8 ) InvJac[2][0]=2./(xxb_qnds[1]-xxb_qnds[0]); else InvJac[2][0]=0.;
			InvJac[2][0]=1.;
			InvJac[2][1]=0.;InvJac[2][2]=0.;InvJac[2][3]=0.;
			_fe[2]->get_dphi_gl_g(ndim-1,qp,InvJac[2],_dphi_g[2]);   // global coord d_phi_g

                        
#ifdef AXISYM   // axisymmetric  (index ->0)
                        interp_el_sol(_data_eq[2].ub,0,DIMENSION,_phi_g[2],elb_ndof[2],_ub_g[2]);
                        JxW_g[2]  *=_ub_g[2][0];
#endif

                        // ***********************************
                        // local side loop (over the node face)
                        for (int lsi_node=0; lsi_node< elb_ndof[2]; lsi_node++) {

                            // set up row i
                            const double phii_g=_phi_g[2][lsi_node]; // boundary test function
                            for (int idim=0; idim< ndim-1; idim++) 
			      dphiidx_g[2][idim*elb_ndof[2]]=_dphi_g[2][lsi_node+idim*elb_ndof[2]]; //];
                            int lei_node= sur_toply[lsi_node]; // local element index
                            int bc_s=(int)_bc_bd[lei_node];
                            int bc_v=(int)_bc_vol[lei_node];
                            

                            // flag setup: +++++++++++++++++++++++++++++++++++++++
                            //  \int (bc_var*T+bc_val*val)ds
                            int  bc_val = (int)((bc_s&2)>>1);  // (1?) non-homogeneous
                            int  bc_var = (int)(bc_s%2);       // (?1) tau=A*variable
                            // ++++++++++++++++++++++++++++++++++++++++++++++++++++
			    
			    double dtxJxW_g=_dt*JxW_g[2]*bc_val; // Neumann bc flag and Jac

                            // Assemblying Matrix ---------------------------------
                            for (int lsj_node=0; lsj_node< elb_ndof[2];  lsj_node++) {
			      double Lap=0.;
			      for (int idim=0; idim< ndim-1; idim++) {
				  dphijdx_g[2][idim*elb_ndof[2]]=_dphi_g[2][lsj_node+idim*elb_ndof[2]];//];
				  Lap += alpha_eff*dphijdx_g[2][idim*elb_ndof[2]]*dphiidx_g[2][idim*elb_ndof[2]];            // Laplacian
			      }
			      
                                KeM(lei_node,sur_toply[lsj_node]) += dtxJxW_g*(
									  _beta*phii_g*_phi_g[2][lsj_node]
									  +_lambda*Lap
									      ); 
                            }// end j  ---------------------------------------

                        }// i   +++++++++++++++++++++++++++++++
                    } // end of the quadrature point qp-loop **********************
                    
//                     if (xxb_qnds[0]< 1.+ BDRY_TOLL) {
// 		      FeM(sur_toply[0]) += Ipenalty*0.5;
// 		      KeM(sur_toply[0],sur_toply[0]) += Ipenalty; }
//                     if (xxb_qnds[1]> 2.- BDRY_TOLL) {
//                       FeM(sur_toply[1]) += Ipenalty*0.5;
// 		      KeM(sur_toply[1],sur_toply[1]) += Ipenalty;}

                } // End Equation on the boundary

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
//     A[Level]->print();
//     b[Level]->print();
#ifdef PRINT_INFO
    std::cout<< " Matrix Assembled(T)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif

    return;
}


// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolT_G::MGTimeStep(
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


#endif
// #endif // personal application

