#include "Equations_conf.h"

// ============================================
#ifdef T_COUP_EQUATIONS // 3D-2D Energy adjoint equation
#if T_COUP_EQUATIONS == 2  // Discontinuous Galerkin
// ============================================

// class local configuration -------
#include "MGSolverT_COUP.h"
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
MGSolT_COUP::MGSolT_COUP(MGEquationsSystem& mg_equations_map_in,
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
//   _h_conv(_mgphys.get_par("hconv")),  // parameter  heat transfer convection coefficient
//   _T_inf(_mgphys.get_par("T_inf"))
{
    //  =================================================
    /// B) setting class variables
    

    _var_names= new std::string[_n_vars];            // names
    _refvalue = new double[_n_vars];
    
    if (_n_vars == 1) {
      _var_names[0]=varname_in;
      _refvalue[0]=_Tref;
    }
    else {
      _var_names[0]="T_ad";
      _var_names[1]="T";
      _refvalue[0]=_Tref;
      _refvalue[1]=_Tref;
    }
    
    
    /// C ) setting solver type
    for (int l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVERT);

    /// D) setting nondimensional parameters
    _alpha=_kappa0/(_rhof*_cp0);
    _IPrdl=_rhof*_alpha/_muf;
    _IRe=_muf/(_rhof*_uref*_lref);
    _IPrdl_turb=1.;
    _alpha_turb=0.;
//   _Nusselt=(_h_conv*_lref)/_kappa0;
    _qheat=mg_equations_map_in.get_par("qheat")*_lref/(_rhof*_cp0*_Tref*_uref);
    _qs=   mg_equations_map_in.get_par("qs")/(_rhof*_cp0*_Tref*_uref);
    _beta= mg_equations_map_in.get_par("beta");
    
    return;
}



//  ====================================================
/// This function assembles the matrix and the rhs:
//  ====================================================
void  MGSolT_COUP::GenMatRhs(
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
    const int elb_ngauss = _fe[2]->_NoGauss1[ndim-2];             // bd elem gauss points
    double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];             // Jac, Jac*w Jacobean
    double dphijdx_g[3][DIMENSION];
    double dphiidx_g[3][DIMENSION];   // global derivatives at g point
    
    //sub-element connectivity
    double    coor_subdom[DIMENSION*NDOF_FEM];          // sub-element node coords
    int       sub_el_nod[NSUBDOM*NDOF_P];               // sub elements division
    int       int_sides_nod[NDOF_FEMB*12];              // internal sides nodes of sub elements
    int       ext_sides_nod[NDOF_FEMB*2*12];            // external sides nodes of sub elements
    int       mat_i[2*(NDOF_P+(DIMENSION-2)*4)];        // indeces for the local matrix when assemblying internal sides fluxes 
    int       bou_i[2*(NDOF_P+(DIMENSION-2)*4)];        // indeces for the local matrix when assemblying external sides fluxes (boundary)
    double    normal[DIMENSION];                        // normal to the boundary
//     double    vel[DIMENSION];                           // velocity on the boundary
    double    vxn;                                      // scalar product v dot n

  
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
    int ns_TADJ=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[TA_F]];    // Temperature Adjoint
//     int ns_T=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F]];    // Temperature
//     int ns_TADJ=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[KTT_F]];    // Temperature adjoint
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
	
	init_subelnod(sub_el_nod,int_sides_nod,ext_sides_nod,mat_i,bou_i); //sub element connectivity and topology

        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc(Level,iel+ndof_lev,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

        // grid and gaussian points
        for (int idim=0; idim<DIMENSION; idim++) {
            for (int d=0; d< NDOF_FEM; d++) _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];// element nodes xxg (DIM)
            // element grid distance
            h_eff[idim]=fabs(xx_qnds[idim*NDOF_FEM+2+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+0]);
            double h_1=fabs(xx_qnds[idim*NDOF_FEM+3+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+1]);
            if (h_eff[idim] <  h_1) h_eff[idim]=h_1; // Max dx diagonal term
            f_upwind[idim]=h_eff[idim];                       // upwind
        }
        
        // quadratic field values
            for (int eq=0; eq<_data_eq[2].n_eqs; eq++) {
                _data_eq[2].mg_eqs[eq]->get_el_sol(0,_data_eq[2].indx_ub[eq+1]-_data_eq[2].indx_ub[eq],
                                                     el_ndof[2],el_conn,offset,_data_eq[2].indx_ub[eq],_data_eq[2].ub);
            }
        // linear field values
// 	_data_eq[1].mg_eqs[1]->get_el_sol(1,1,_el_dof[1],el_conn,offset,1,_data_eq[1].ub);//take temperature with pressure (NS)
	_data_eq[1].mg_eqs[0]->get_el_sol(0,1,_el_dof[1],el_conn,offset,1,_data_eq[1].ub);//take temperature with no pressure (no NS)
        
        //  external node quantities -------------------------------------

        //  external cell properties -------------------------------------
        double phase=1.;
        double Qad =_qheat;
        double rhocp=1.;
        double rho=1.;
        double qn =_qs;
	double alpha_eff = 1.;//_IPrdl*_IRe;
// 	if (xx_qnds[8]< 3.) qn=0.; else qn=0.1;
//   	_dt =0.5;
// 	vel_g[0]=4.*(xx_qnds[17]-xx_qnds[17]*xx_qnds[17]); vel_g[1]=0.; vel_g[2]=0.; //double boundary layer
	vel_g[0]=-xx_qnds[8]*xx_qnds[8]+2.*xx_qnds[8]; vel_g[1]=0.; vel_g[2]=0.;  //bottom boundary layer
// 	vel_g[0]=-xx_qnds[8]*xx_qnds[8]+1.; vel_g[1]=0.; vel_g[2]=0.;             //top boundary layer
	
// 	vel_g[0]=1.; vel_g[1]=0.; vel_g[2]=0.;  
// 	cont_zone=0.;
// 	cont_zone=4.*(xx_qnds[17]-xx_qnds[17]*xx_qnds[17])*xx_qnds[8];  //for Paraview  4*(coordsY-coordsY*coordsY)*coordsX
// 	cont_zone=16.*(xx_qnds[17]-xx_qnds[17]*xx_qnds[17])*
// 		      (xx_qnds[8]-xx_qnds[8]*xx_qnds[8]);  //for Paraview 16*(coordsY-coordsY*coordsY)*(coordsX-coordsX*coordsX)
// 	cont_zone=(xx_qnds[17]-xx_qnds[17]*xx_qnds[17])*
// 		  (8.*xx_qnds[8]-xx_qnds[8]*xx_qnds[8])/4.;  //for Paraview 0.5*(T-0.7)*(T-0.7)*(coordsY-coordsY*coordsY)*(8*coordsX-coordsX*coordsX)/4
// 		  							    
	
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

            

            //  fields -----------------------------------------------------------
            interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],
                          _phi_g[2],el_ndof[2],_ub_g[2]); // quadratic
	    interp_el_sol(_data_eq[1].ub,0,2,_phi_g[1],el_ndof[1],_ub_g[1]); // linear
#ifdef AXISYM   // axisymmetric (index -> 0)
            JxW_g[2]  *=_ub_g[2][0];
#endif
// 	    cont_zone = (_ub_g[2][1]-_ub_g[2][1]*_ub_g[2][1])*
// 		  (8.*_ub_g[2][0]-_ub_g[2][0]*_ub_g[2][0])/4.;
            vel_g[0]=-_ub_g[2][1]*_ub_g[2][1]+2.*_ub_g[2][1];
            /// d) Local (element) assemblying energy equation
            // *********************** *******************************
            for (int i=0; i<el_ndof[2]; i++)     {
                // set up row i
	        double Adv_theta=0.;
                const double phii_g=_phi_g[2][i];
                for (int idim=0; idim< ndim; idim++) {  
		  dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
// 		  Adv_theta += _ub_g[2][ns_NS+idim]*dphiidx_g[2][idim];
		  Adv_theta += vel_g[idim]*dphiidx_g[2][idim];
		}

                double dtxJxW_g=_dt*JxW_g[2]*_bc_vol[i];


                // Rhs Assemblying  ----------------------------------------------------------
                if (mode == 1) {
//                     // rhs energy
//                     FeM(i) += dtxJxW_g
//                                   *(rhocp*_ub_g[2][ns_T]*phii_g/_dt)   // time
//                               ;
		    // rhs adjoint
		    FeM(i) += dtxJxW_g*(
                                  -rhocp*_ub_g[2][ns_TADJ]*phii_g/_dt   // time
                                   + (T_CONTROL-_ub_g[1][1])*cont_zone*phii_g            // source
                              );
                }


                // Matrix Assemblying ---------------------------
                for (int j=0; j<el_ndof[2]; j++) {
                    double phij_g= _phi_g[2][j];
//                     double Adv=0.;
                    double Lap=0.;
                    for (int idim=0; idim< ndim; idim++) {
                        dphijdx_g[2][idim]=_dphi_g[2][j+idim*el_ndof[2]];
// 			Adv += _ub_g[2][ns_NS+idim]*dphijdx_g[2][idim];          // advection
//                         Adv +=vel_g[idim]*dphijdx_g[2][idim];                    // advection
                        Lap +=(alpha_eff)      //+0.*f_upwind[idim]*vel_g[idim]*vel_g[idim]) // upwind
                              *dphijdx_g[2][idim]*dphiidx_g[2][idim];            // Laplacian
//             Grad += dphijdx_g[2][idim]*phii_g;            // Gradient
                    }
                    // energy-equation
//                     KeM(i,j) +=dtxJxW_g*(
//                                    rhocp*phii_g*phij_g/_dt // time term
//                                    +rhocp*Adv*phii_g      //adv
//                                    + Lap                   //diff
//                                );
		    // adjoint-equation
                    KeM(i,j) +=dtxJxW_g*(
                                                   - rhocp*phii_g*phij_g/_dt // time term
                                                   - rhocp*Adv_theta*phij_g      //adv
                                                   - Lap                   //diff
                               );
		  
		}
            } // ----------------------------------------
        } // end of the quadrature point qp-loop ***********************


        
        
        
        
        
// ---------------------------------- BEGIN T discontinous Galerkin assemblying

            // ================ diagonal terms ================
	
            for (int i=0; i<_el_dof[1]; i++)     { //    ----> piecewise in linear node--
              for(int sub=0;sub<NDOF_P;sub++) 
               for (int dim=0;dim<DIMENSION;dim++) 
		 coor_subdom[sub+dim*NDOF_FEM]=xx_qnds[sub_el_nod[i*NDOF_P+sub]+dim*NDOF_FEM];
		 
	      for (int qp=0; qp< el_ngauss; qp++) {  	       
                  det[1]  = _fe[1]->Jac(qp,coor_subdom,InvJac[1]);     // determinant Jacobian   
		  JxW_g[1]= _dt*det[1]*(_fe[1]->_weight1[ndim-1][qp]);
//                 
                  // Rhs Assemblying  ----------------------------------------------------------
                  if (mode == 1) FeM(i+el_ndof[2]) += JxW_g[1]*_data_eq[1].ub[i+NDOF_FEM]/_dt;  // time
                  
                         
                  // Matrix Assemblying ---------------------------
                  KeM(i+el_ndof[2],i+el_ndof[2]) +=  JxW_g[1]/_dt;         //time
            
	        }//gauss loop
            }//row-subelement loop
       
       
            // ================ off diagonal terms ================
            int plane=0;
	    int N_intside=4;
#if DIMENSION==3
          for(plane=0;plane<3;plane++) { //loop on planes of hex27: 0 half bottom volume, 1 center planes, 2 half top volume
#endif
            for(int side=0;side<N_intside;side++) { //loop on internal sides of sub-elements
              double Adv=0.; double Lap=0.;
	      int i=plane*N_intside*2+side*2;
	      int j=plane*N_intside*2+side*2+1;
	      for(int dot=0;dot<NDOF_FEMB;dot++) // coordinates of the internal side
		for (int dim=0;dim<DIMENSION;dim++)
	         xxb_qnds[dot+dim*NDOF_FEMB]=xx_qnds[int_sides_nod[plane*N_intside*NDOF_FEMB+side*NDOF_FEMB+dot]+dim*NDOF_FEM];  
		
		
		double h_vec[DIMENSION]; double h=0.;  vxn=0.; 
	       for(int idir=0;idir<DIMENSION;idir++) {
		 h_vec[idir]= _data_eq[2].ub[mat_i[j]+idir*NDOF_FEM]-_data_eq[2].ub[mat_i[i]+idir*NDOF_FEM];
		 normal[idir]=0.;
	       }
		
		_fe[1]->normal_g(xxb_qnds,normal); // normal
		for(int idir=0;idir<DIMENSION;idir++) { //velocity
// 		  if ((normal[idir]*_data_eq[2].ub[mat_i[i]+(idir+ns_NS)*NDOF_FEM])>0) vel_g[idir]=_data_eq[2].ub[mat_i[i]+(idir+ns_NS)*NDOF_FEM];
// 		  else vel_g[idir]=_data_eq[2].ub[mat_i[j]+(idir+ns_NS)*NDOF_FEM];
		}
		for(int idir=0;idir<DIMENSION;idir++) {vxn+=normal[idir]*vel_g[idir]; h += normal[idir]*h_vec[idir];}

	      for (int  gpb=0; gpb<  elb_ngauss; gpb++) { //gaussian integration loop on internal side 
		  det[1]  = _fe[1]->JacSur(gpb,xxb_qnds,InvJac[1]);   // local coord jac
		  JxW_g[1]= _dt*det[1]*_fe[1]->_weight1[ndim-2][gpb]; // weight
		  Adv += JxW_g[1]*vxn;
		  Lap += -JxW_g[1]*alpha_eff/h;
	      }
	      
	      if (vxn>-1.e-6) {
// 		if (mode == 1) FeM(mat_i[i]+el_ndof[2])-=Adv*_data_eq[1].ub[mat_i[i]];   //explicit method
// 		if (mode == 1) FeM(mat_i[j]+el_ndof[2])+=Adv*_data_eq[1].ub[mat_i[i]];   //explicit method
		KeM(mat_i[i]+el_ndof[2],mat_i[i]+el_ndof[2])-= Adv;                   //implicit method
		KeM(mat_i[i]+el_ndof[2],mat_i[j]+el_ndof[2])+= Adv;                   //implicit method
	      }
	      else {
// 		if (mode == 1) FeM(mat_i[i]+el_ndof[2])-=Adv*_data_eq[1].ub[mat_i[j]];   //explicit method
//   	        if (mode == 1) FeM(mat_i[j]+el_ndof[2])+=Adv*_data_eq[1].ub[mat_i[j]];   //explicit method
		KeM(mat_i[i]+el_ndof[2],mat_i[i]+el_ndof[2])-= Adv;                   //implicit method
		KeM(mat_i[i]+el_ndof[2],mat_i[j]+el_ndof[2])+= Adv;                   //implicit method
	      }
	      
	      //Laplacian term
	      KeM(mat_i[i]+el_ndof[2],mat_i[j]+el_ndof[2]) += Lap;
	      KeM(mat_i[i]+el_ndof[2],mat_i[i]+el_ndof[2]) -= Lap;
	      
	      KeM(mat_i[j]+el_ndof[2],mat_i[i]+el_ndof[2]) += Lap;
	      KeM(mat_i[j]+el_ndof[2],mat_i[j]+el_ndof[2]) -= Lap;
	      
	    }//loop on internal sides of sub-elements
#if DIMENSION==3
 	 }//plane loop
#endif	    
// 	    
// // 	 // ====================== boundary ======================================

	double theta_dx[NDOF_FEM*(DIMENSION)];
	for(int init=0; init< NDOF_FEM*DIMENSION; init++) theta_dx[init]=0.;
                
		for(int inode=0; inode< NDOF_FEM; inode++)  {
		  det[2] = _fe[2]->Jac_nodes(inode,xx_qnds,InvJac[2]);     // Jacobian
                  _fe[2]->get_dphi_node(ndim,inode,InvJac[2],_dphi_g[2]); // global coord deriv
		  for (int idim=0; idim<DIMENSION; idim++) {
		    for(int jnode=0; jnode< NDOF_FEM; jnode++)  {
		    theta_dx[inode+idim*NDOF_FEM]+=_data_eq[2].ub[ns_TADJ*NDOF_FEM+jnode]*_dphi_g[2][jnode+idim*NDOF_FEM];
		    }
		  }
		}


	 int N_extside=2;
#if DIMENSION==3
	     N_extside=4;
#endif	    
	     for (int iside=0; iside< el_sides; iside++)  {
                 if (el_neigh[iside] == -1) {
		   
		  for(int midside=0; midside<N_extside; midside++){ //external sides, divided in two(2D) (0-4),(4-1) <-iside=0 or in four (3D)
		   
		    for(int iiside=0;iiside<NDOF_FEMB;iiside++) // coordinates of the external side
		     for (int dim=0;dim<DIMENSION;dim++) {
		       int indx=iside*N_extside*NDOF_FEMB+midside*NDOF_FEMB+iiside;
	               xxb_qnds[iiside+dim*NDOF_FEMB]=
	               xx_qnds[ext_sides_nod[indx]+dim*NDOF_FEM]; 
		     }
		     
		     for (int idir=0;idir<DIMENSION;idir++) normal[idir]=0.;
		     
		      int i= bou_i[iside*N_extside+midside]; //dof to be filled
		       vxn=0.; double Adv=0.;
		       double Ipenalty= det[1]/_dt;
		       _fe[1]->normal_g(xxb_qnds,normal); // normal
// 		       for(int idir=0;idir<DIMENSION;idir++) vel_g[idir]=_data_eq[2].ub[i+(idir+ns_NS)*NDOF_FEM]; //velocity
		       for(int idir=0;idir<DIMENSION;idir++) vxn += vel_g[idir]*normal[idir];
// 			  
		      for (int  gpb=0; gpb<  elb_ngauss; gpb++) { //gaussian integration loop on internal side 
			  det[1] = _fe[1]->JacSur(gpb,xxb_qnds,InvJac[1]);   // local coord jac
			  JxW_g[1]=det[1]*_dt*_fe[1]->_weight1[ndim-2][gpb]; // weight
			  Adv+=JxW_g[1]*vxn;
		      }	       
		      if (_bc_vol[i+el_ndof[2]]!=0) { //bc_vol != 0 
			if (vxn>=-1.e-6) {  
// 			  if (mode == 1) FeM(i+el_ndof[2]) -= Adv*_data_eq[1].ub[i];  // explicit
			}
			else { //assign the value if the flux is entering from outside
			  KeM(i+el_ndof[2],i+el_ndof[2])   += Ipenalty;
			  if (mode == 1) FeM(i+el_ndof[2]) += Ipenalty*0.7;
			}   
		      }
		      else { //bc_vol = 0 
			
			if(_bc_bd[i+el_ndof[2]]==0) {
			  double thetaxn=0.;
			  for (int idim=0; idim<DIMENSION; idim++) thetaxn += theta_dx[i+idim*NDOF_FEM]*normal[idim];
			    KeM(i+el_ndof[2],i+el_ndof[2])  += Ipenalty;
			    if (mode == 1) FeM(i+el_ndof[2])+= Ipenalty*alpha_eff*thetaxn/_beta;
			  }
			  else{
			    KeM(i+el_ndof[2],i+el_ndof[2])  += Ipenalty;
			    if (mode == 1) FeM(i+el_ndof[2])+= Ipenalty*0.5; 
			  }
		      } 
		  }// midside loop
	        } //if neighbor
	    }
// // ----------------------------------END T discontinous Galerkin assemblying
        

        
        // ======================================================================
        // ====================== boundary ======================================
        // ======================================================================
	
        for (int iside=0; iside< el_sides; iside++)  {
            if (el_neigh[iside] == -1) {


                for (int idof=0; idof<NDOF_FEMB; idof++) {
                    sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
                    int idofb=sur_toply[idof];
                    elb_conn[idof]=el_conn[idofb];                 // connectivity vector
                    for (int idim=0; idim<DIMENSION; idim++) {
                        xxb_qnds[idim*NDOF_FEMB+idof]=xx_qnds[idim*NDOF_FEM+idofb]; // coordinates
                        _data_eq[2].ub[idim*NDOF_FEM+idof]=xxb_qnds[idim*NDOF_FEMB+idof];
                    }
                }

                              
                // Dirichlet boundary conditions for theta ***********************************
//         printf(" \n  bc vol %d  bd %d ",_bc_vol[sur_toply[NDOF_FEMB-1]],_bc_bd[sur_toply[NDOF_FEMB-1]]);
                if (_bc_vol[sur_toply[NDOF_FEMB-1]] == 0) {
		     //[NDOF_FEMB-1] is the midpoint of a quadratic FEM element (HEX9 (2D) or HEX27 (3D))
                    double det= _fe[2]->JacSur(elb_ngauss-1,xxb_qnds,InvJac[2]);// jacobian
                    double Ipenalty= det/_dt;                               // Dirichlet bc flag
                    // local boundary loop   ---------------------------------------
                    for (int lb_node=0; lb_node< elb_ndof[2]; lb_node++) {
                        int lv_node= sur_toply[lb_node]; // local vol index
                        int  bc_s_theta = (int) _bc_bd[lv_node];
			int  bc_var_theta = (int)(bc_s_theta%2);   
                        if (mode == 1) {

			    FeM(lv_node) += bc_var_theta*Ipenalty*_data_eq[2].ub[(ns_TADJ)*NDOF_FEM+lv_node];
			    FeM(lv_node) += (1.-bc_var_theta)*Ipenalty*0.; //theta value is 0
                        }

                        KeM(lv_node,lv_node) += Ipenalty;  //  Dirichlet bc
                    }// lb_node -end  local boundary loop -------------------------
		} // end if Dirichlet  boundary conditions
                // **********************************************************************
            } //end if side
        } // ======================  end for boundary =======================================

        /// e) Global assemblying energy equation
        A[Level]->add_matrix(KeM,el_dof_indices);                  // global matrix
        if (mode == 1)   b[Level]->add_vector(FeM,el_dof_indices); // global rhs

    } // end of element loop
    // clean
    el_dof_indices.clear();
    A[Level]->close();
//     A[Level]->print();
    if (mode == 1) {b[Level]->close(); /*b[Level]->print();*/}
#ifdef PRINT_INFO
    std::cout<< " Matrix Assembled(T)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif

    return;
}


// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolT_COUP::MGTimeStep(
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


void MGSolT_COUP::init_subelnod(
  int sub_el_nod[],
  int int_sides_nod[],
  int ext_sides_nod[],
  int mat_index[],
  int boun_index[]
) {
int N_intside, N_extside, N_sub, N_mat, N_bou;
#if DIMENSION==2
N_sub=NDOF_P*NSUBDOM;
N_intside=NDOF_FEMB*NSUBDOM;
N_extside=NDOF_FEMB*NSUBDOM*2;
N_mat=2*NDOF_P;
N_bou=2*NDOF_P;

int map[NDOF_P*NSUBDOM]={
  0, 4, 8, 7,
  4, 1, 5, 8,
  8, 5, 2, 6,
  7, 8, 6, 3
};
int int_side[NDOF_FEMB*NSUBDOM]={
  4,8,0, 5,8,0, 6,8,0, 7,8,0                               // internal sides
};
int mat_ind[2*NDOF_P]={
  0,1, 1,2, 2,3, 3,0
};
int ext_side[NDOF_FEMB*NSUBDOM*2]={
  0,4,0, 4,1,0, 1,5,0, 5,2,0, 2,6,0, 6,3,0, 3,7,0, 7,0,0    // external sides
};
int bou_ind[2*NDOF_P]={
  0,1, 1,2, 2,3, 3,0
};

#else
N_sub=NDOF_P*NSUBDOM;
N_intside=NDOF_FEMB*12;
N_extside=NDOF_FEMB*4*6;
N_mat=2*12;
N_bou=4*6;

int map[NDOF_P*NSUBDOM]={
  0, 8, 20, 11, 12, 21, 26, 24,
  8,  1, 9, 20, 21, 13, 22, 26,
  20, 9, 2, 10, 26, 22, 14, 23,
  11, 20, 10, 3, 24, 26, 23, 15,
  
  12, 21, 26, 24, 4, 16, 25, 19,
  21, 13, 22, 26, 16, 5, 17, 25,
  26, 22, 14, 23, 25, 17, 6, 18,
  24, 26, 23, 15, 19, 25, 18, 7
};
int int_side[NDOF_FEMB*12]={
  8,20,26,21,0,0,0,0,0,   9,20,26,22,0,0,0,0,0, 10,20,26,23,0,0,0,0,0, 11,20,26,24,0,0,0,0,0,  //bottom volume
  12,21,26,24,0,0,0,0,0, 21,13,22,26,0,0,0,0,0, 22,14,23,26,0,0,0,0,0, 23,15,24,26,0,0,0,0,0,  //central horizontal surface
  21,26,25,16,0,0,0,0,0, 22,26,25,17,0,0,0,0,0, 23,26,25,18,0,0,0,0,0, 24,26,25,19,0,0,0,0,0   //top volume
};
int mat_ind[2*12]={
  0,1, 1,2, 2,3, 3,0,           //bottom volume
  0,4, 1,5, 2,6, 3,7,           //central horizontal surface
  4,5, 5,6, 6,7, 7,4            //top volume
};
int ext_side[NDOF_FEMB*4*6]={
  0,11,20,8,0,0,0,0,0,   8,20,9,1,0,0,0,0,0,  20,10,2,9,0,0,0,0,0,  11,3,10,20,0,0,0,0,0, // bottom surface
  0,8,21,12,0,0,0,0,0,   8,1,13,21,0,0,0,0,0, 21,13,5,16,0,0,0,0,0, 12,21,16,4,0,0,0,0,0, // front surface
  1,9,22,13,0,0,0,0,0,   9,2,14,22,0,0,0,0,0, 22,14,6,17,0,0,0,0,0, 13,22,17,5,0,0,0,0,0, // right surface
  2,10,23,14,0,0,0,0,0, 10,3,15,23,0,0,0,0,0, 23,15,7,18,0,0,0,0,0, 14,23,18,6,0,0,0,0,0, // behind surface
  3,11,24,15,0,0,0,0,0, 11,0,12,24,0,0,0,0,0, 24,12,4,19,0,0,0,0,0, 15,24,19,7,0,0,0,0,0, // left surface
  4,16,25,19,0,0,0,0,0, 16,5,17,25,0,0,0,0,0, 25,17,6,18,0,0,0,0,0, 19,25,18,7,0,0,0,0,0, // top surface
};
int bou_ind[4*6]={
  0,1,2,3, 0,1,5,4, 1,2,6,5, 2,3,7,6, 3,0,4,7, 4,5,6,7
};
#endif

//assignment
for(int i=0;i<N_sub;i++)     sub_el_nod[i]=map[i];
for(int i=0;i<N_intside;i++) int_sides_nod[i]=int_side[i];
for(int i=0;i<N_extside;i++) ext_sides_nod[i]=ext_side[i]; 
for(int i=0;i<N_mat;i++)     mat_index[i]=mat_ind[i];
for(int i=0;i<N_bou;i++)     boun_index[i]=bou_ind[i];

return;
}


#endif
#endif