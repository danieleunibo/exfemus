#include "Equations_conf.h"
// ============================================
#ifdef COLOR_EQUATIONS // 3D-2D Energy equation
// ============================================

#include "MGSolverCOL.h"
// #include "MGSFSIconf.h"

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



// ======================================================
/// This function constructs the 3d-2D MGSolCOLX_Y_Z class 
MGSolCOL::MGSolCOL(MGEquationsSystem& mg_equations_map_in,
		       const int vars_in[],
                       std::string eqname_in,
                       std::string varname_in):
    MGSolDA(mg_equations_map_in,vars_in,eqname_in,varname_in),
    // mesh params ------------
    _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes
    // phys params ------------
    _dt(_mgutils.get_par("dt")),       // parameter  dt
    _rhof(mg_equations_map_in.get_par("rho0")),    // parameter  density reference
    _uref(mg_equations_map_in.get_par("Uref")),    // parameter  vel reference
    _lref(mg_equations_map_in.get_par("Lref")),    // parameter  length reference
    _Tref(mg_equations_map_in.get_par("Tref")),    // parameter  temperature reference
    _rhos(mg_equations_map_in.get_par("rhos"))     // parameter solid density
{ //  =================================================
  // class variable names
  _dir=0;
  
  if(!varname_in.compare("dx")) _dir=0;
  if(!varname_in.compare("dy")) _dir=1;
  if(!varname_in.compare("dz")) _dir=2;
#if DIMENSION ==2
  if (_dir==2) std::cout<<"Too many Dimension!!\n";
#endif
  _var_names[0]=varname_in;     _refvalue[0]=_lref;
  
  for (int l=0;l<NDOF_FEM;l++) {_bc_vol[l]=-1; _bc_bd[l]=-1;}
  
  // class solver type (SOLVERT  in MGSTconf.h)
//   for (int l=0;l<_NoLevels;l++) _solver[l]->set_solver_type(CGM);
     for (int l=0;l<_NoLevels;l++) _solver[l]->set_solver_type(GMRESM);
    
  return;
}


//  ====================================================
/// This function assembles the matrix and the rhs:
//  ====================================================
void  MGSolCOL::GenMatRhs(
  const double /* time*/, // time  <-
  const int Level,  // Level <-
  const  int mode    // mode  <- (1=rhs) (0=only matrix)
) {  // ===============================================
// not nedeed here
  
  return;}
















//  ====================================================
/// This function assembles the matrix and the rhs:
//  ====================================================
void  MGSolCOL::ColorFunc(
  const double /* time*/, // time  <-
  const int Level,  // Level <-
  const  int mode    // mode  <- (1=rhs) (0=only matrix)
) {  // ===============================================

  /// Set up
  // geometry -----
  const int  ndim = DIMENSION;                                           //dimension
  double     xx_qnds[DIMENSION*NDOF_FEM], xxb_qnds[DIMENSION*NDOF_FEMB]; // element node coords
  int        el_conn[NDOF_FEM];//, elb_conn[NDOF_FEMB];                     // element connectivity
  int        el_neigh[NDOF_FEM];                                         // element connectivity
  const int  offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
  const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
  int        sur_toply[NDOF_FEMB]; // boundary topology
//   double     vel_g[DIMENSION];  for(int idim=0;idim<DIMENSION;idim++) vel_g[idim] =0.; 
//   double  val_tbg[3];
  int fl_int;
  int phase;
  int interface;
//   int interface;
//   double rho;
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
   // gauss integration  -----------------------------------------------------------------------------
    const int el_ngauss = _fe[2]->_NoGauss1[ndim-1];         // quadratic element gauss points
    const int elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];   // quadratic bd elem gauss points
    
    // element dofs (costant[0]-linear[1]-quadratic[2]) -----------------------------------------------
    int elb_ndof[3];
    elb_ndof[0]=NDOF_BK; elb_ndof[1]=NDOF_PB; elb_ndof[2]=NDOF_FEMB;  // number of element boundary dofs 
    int el_mat_nrows =0;                                              // number of matrix rows (dofs)
    for (int ideg=0; ideg<3; ideg++) el_mat_nrows +=_nvars[ideg]*_el_dof[ideg];
    int el_mat_ncols = el_mat_nrows;                    // square matrix
    std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector

       
    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix+rhs) ---------------------------
    DenseVectorM XeM;
    XeM.resize(el_mat_nrows);
    // number of total elements for level
    int ndof_lev=0;
        for (int pr=0;pr <_mgmesh._iproc; pr++) {
            int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
            ndof_lev +=delta;
        }
        
        
    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
    for (int iel=0; iel < (nel_e - nel_b); iel++) {
    
	XeM.zero();     // set to zero matrix and rhs

        // geometry and element  fields ------------------------------------
        
        _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);   // Element Connectivity (el_conn) and coordinates (xx_qnds)
        _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh); // Neighbors of the element

        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc(Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

        // ======================================================================
        // Volume =============================================================
        // ======================================================================
for (int n_dof_count =0; n_dof_count< el_mat_nrows; n_dof_count++)
  XeM(n_dof_count)=10.;
 // *********************** *******************************
    phase=(_bc_vol[NDOF_FEM-1]<2)?0:1; 
    for(int i=1;i<el_sides+1;i++)  if(_bc_vol[NDOF_FEM-i] > 4.5 ) interface=1.;
            for (int i=0; i<_el_dof[2]; i++)     {   //  --- QUADRATIC ---
                for (int ivar=0;ivar<_nvars[2];ivar++) {
		   int index=i+ivar*_el_dof[2];
                        XeM(index) = 2.*phase;
		}
	    }
            for (int i=0; i<_el_dof[1]; i++)     { //    --- LINEAR ---
                for (int ivar=0;ivar<_nvars[1];ivar++) {
                    int index=i+ivar*_el_dof[1]+_el_dof[2]*_nvars[2];
                        XeM(index) =2.*phase ;
                    }
	    }
            for (int i=0; i<_el_dof[0]; i++)     { //    --- Piecewise ---
                for (int ivar=0;ivar<_nvars[0];ivar++) {
                    int index=i+ivar*_el_dof[0]+_el_dof[2]*_nvars[2]+_el_dof[1]*_nvars[1];
                        XeM(index) = 6.*phase;
		}
	    }
for (int a=0; a < el_mat_nrows ; a++)
  x[Level]->set(el_dof_indices[a],XeM(a));



  } // end of element loop
  // clean
  el_dof_indices.clear();
#ifdef PRINT_INFO
  std::cout<< " Matrix Assembled(T)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif


  
  return;
}
















// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolCOL::MGTimeStep(const double time, const int /*iter*/) { // ------------------------------------

// #ifdef AXISYM
// if(_dir==1) return;
// #endif

  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;
  /// B) Assemblying of the rhs (top level in surface and volume with MGSolNS::GenRhs,GenRhsB),
#if PRINT_TIME==1
  std::clock_t start_time=std::clock();
#endif
  ColorFunc(time,_NoLevels-1,1);
//   GenMatRhs(time,_NoLevels-1,1);
//   A[_NoLevels-1]->close();
  /// C) Assemblying of the  MGmatrices (MGSolNS::GenMatrix),
  for(int Level = 0 ; Level < _NoLevels-1; Level++) {
//     GenMatRhs(time,Level,0); // matrix
//     A[Level]->close();
  }
#if PRINT_TIME==1
  std::clock_t end_time=std::clock();
  std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif

  /// E) Solution of the linear MGsystem (MGSolNS::MGSolve).
//   MGSolve(1.e-6,40);
#if PRINT_TIME==1
  end_time=std::clock();
  std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
            << "s "<< std::endl;
#endif
  /// A) Update of the old solution at the top Level  (MGSolNS::OldSol_update),
//   x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);
 
    //update x_oold
//   x_oold[_NoLevels-1]->zero();
//   x_oold[_NoLevels-1]->add(1,*x_old[_NoLevels-1]); 

//   //update x_old
   x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);
// 
  
//   const int flag_moving_mesh = _mgutils.get_par("moving_mesh");
//   
// //    //update dx_old ========================================
//    if(flag_moving_mesh) {
//     /// E) mesh update    
//    const int n_nodes=_mgmesh._NoNodes[_NoLevels-1];
//    const int offsetp=_dir*n_nodes;
//    for (int inode=0;inode<n_nodes;inode++) { 
//        double disp=(*x_old[_NoLevels-1])(inode)-(*x_oold[_NoLevels-1])(inode);
//       _mgmesh._dxdydz[inode+offsetp] += disp;
//      }
//    }
   // ==============================================================

 return;
}

#endif

