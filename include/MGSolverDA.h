#ifndef __mgsolverda_h__
#define __mgsolverda_h__

#include "Equations_conf.h"
// #ifdef DA_EQUATIONS

// conf includes
#include "MGFE_conf.h"
// algebra class 
#include "dense_vectorM.h"
#include "dense_matrixM.h"
#include "numeric_vectorM.h"

// local class
#include "MGSolverBase.h"   //for the inherited class
#include "MGFEMap.h"
// #include "MGFEMap.h"

// Forward declarations --------------------------
class MGUtils;
class MGSystem;
class MGMesh;
class MGEquationsMap;


// =====================================================
// Table external fields
// ----------------------------------------------------
 
// struct external fields

 /// Class containing all the possible external fields for the system
 class  external_field{
 public: 
   ///@{ \name EXTERNAL FIELDS INFORMATION 
  int n_eqs;                         ///< number of equations 
  int tab_eqs[20];                   ///<  map external system and index
  
  MGSolBase* mg_eqs[20];             ///< external system pointer
  int indx_ub[20];                   ///< index of external equations(const/linear/quad):
  double  ub[20*NDOF_FEM];           ///< element external field old solution
  ///@}
  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  external_field(){}
  ~external_field(){}
  ///@}
  };
// =====================================================


// ====================================================
/// Class for MG Diffusion-Advection equation solvers
// ====================================================

class MGSolDA : public MGSolBase {
  
public:

 // -------------------------------------------------------------
  // -------------------- class field ----------------------------
 
  // ------------------ integration -----------------------
  // shape and derivative function at a  gaussian point
  double _phi_g[3][NDOF_FEM];              ///< shape field (0-1-2 degree)
  double _dphi_g[3][NDOF_FEM*DIMENSION];   ///< shape derivative  (0-1-2 degree)
    
  //  external field 
  external_field _data_eq[3];              ///< external data structure

   // data -----------------------------------------   
  int   _el_dof[3];   ///< number of dof in each element  (piecewise, linear, quadratic)        
  MGFE* _fe[3];       ///< fem  (piecewise linear,piecewise quadratic)   

 // ========================================================================
///@{ \name  Constructor destructor memory allocation
// ========================================================================
  /// Level constructor  I
  MGSolDA(MGEquationsSystem& mg_equations_map,
	  const int nvars_in[],  ///< number of piecewise[0], linear[1], quadratic[2] variables
          std::string eq_name_in="DA",
          std::string varname_in="u");
	  

  /// Destructor (level structure)
  ~MGSolDA();   
 /// Clean all substructures
  void clean();   
  // Setting --------------------------------------------------
  virtual void init_dof( ///< Setting dof
    const int Level          // MG Level
  ); 
  virtual void init(     ///< Setting Memmory alloc
    const int Level         // MG Level
  );   
  ///@} 
  // ========================================================================
 ///@{ \name BOUNDARY/INITIAL CONDITIONS
  // ========================================================================
  void GenBc();         ///< Setting boundary conditions
  void GenBc_loop(
    const int vb,
    const int ndof_femv, 
    const int n_ub_dofs, 
    const int n_pb_dofs,
    const int n_kb_dofs,
    int bc_id[],
    int mat_id[]
  );
  virtual void GenIc(); ///< Element Initial conditions
//   virtual void GenIc0(  ///< Setting Initial conditions
//     const std::string ib_name,
//     const int node_dof_top[],
//     NumericVectorM &sol_top,
//     NumericVectorM &old_sol_top
//   );
///@} 
 
// ========================================================================
 ///@{ \name MULTILEVEL OPERATORS
// ========================================================================  
  // Reading/Writing MG operators ------------------------------------------
  void ReadProl(                    ///< Read Prolongation Op
    const int          Level,       // MG Level <-
    const std::string& name,        // file name for P <-
    SparseMMatrixM     &Mat,        //  Matrix for P <-
    const int          nvars_in[],
    int                node_dof_c[],
    int                node_dof_f[]
  ); 
  
void ReadRest (   ///< Restriction Op.
  const int Level,         // MG Level
  const std::string& name, // file name (reading from)
  SparseMMatrixM &Rest,    // Restriction Matrix
  const int nvars_in[],    // # cost,linear quad variables
  int  node_dof_f[],       // dof map fine mesh
  int  node_dof_c[],       // dof map coarse
  int  _node_dof_top[]     // dof map top level
);    // ----  
  
 
  void ReadMatrix(             ///< Reading Matrix
    const int Level,           // MG Level <-
    const std::string& name,   // file name for M <-
    SparseMatrixM &Mat,
    const int* nvars_in
  );  
///@} 
  
// ========================================================================  
 ///@{ \name    MultiGrid Solver
// ========================================================================  
    virtual void GenMatRhs( ///< Volume Assemblying  matrix-rhs 
    const double time,        // time
    const int Lev,            // Level
    const  int m              // rhs assembly control  
  );   
  /// MG time step solver (backward Euler).
  virtual void MGTimeStep(
    const double time,               // time               <-
    const int    mode                // rhs assembler flag <-   
  ); 
///@}  
  
  // ========================================================================
   ///@{ \name SOLUTION WRITE AND READ 
   // ========================================================================
  // Print
  virtual void print_bc(      ///< Print boundary conditions to a xdmf file.
    std::string namefile,     // filename <-
    const int   Level         // MGLevel  <-
  ); 
  virtual void print_u(      /// Print solution to a xdmf file. 
    std::string namefile,    // filename <-
    const int Level
  );        // MGLevel  <-
  
#ifdef HAVE_MED
  virtual void print_u_med(  /// Print solution to a med file. 
    std::string namefile,    // filename <-
    const int Level);        // MGLevel  <-
#endif
   virtual void print_xml_attrib(
    std::ofstream &out,    // file stream to print ->
    int nodes,
     int nelems,
    std::string file_name
  ) const;
  // Read 
  virtual void read_u(      /// Read solution from a xdmf file
    std::string name,       // filename <-
    int Level              // restart level <-
  );

// ============================================================================
/// This function  defines the boundary conditions for the system:
virtual void bc_intern_read(
  int /*face_id_node*/,  ///<  face identity           (in)
  int  /*mat_flag*/,     ///<  volume identity         (in)
  double /*xp*/[],       ///< xp[] node coordinates    (in) 
  int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
  int bc_flag[]          ///< boundary condition flag  (out)
);
// ============================================================================
/// This function reads Boundary conditions  from function
virtual void bc_read(
   int /*face_id_node*/,  ///<  face identity          (in)
  int  /*mat_flag*/,     ///<  volume identity         (in)
  double /*xp*/[],       ///< xp[] node coordinates    (in) 
  int bc_Neum[],         ///< Neuman (1)/Dirichlet(0)  (out)
  int bc_flag[]          ///< boundary condition flag  (out)
);

virtual void ic_read( /// initial conditions  from function
    int    k,           // bc from gambit    <-
    int    m,           // material from gambit    <-
    double xp[],        // point coordinates <-
    double value[]      //  point values     ->
  );    
  void MGFunctional(double, double){}
///@} 
 
 // Return ----------------------------------------------------
  // ========================================================================
   ///@{ \name INTERPOLATIONS FUNCTIONS
   // ========================================================================
// void get_el_lq(const uint Level,const int nvars[],const int el_nds[],
// 			 const uint el_conn[],const uint offset,
// 			 std::vector<uint>   & el_dof_indices,
//                           std::vector<uint>   bc_bd[],  std::vector<uint>   bc_vol[],		 
//                          std::vector<double> uold[])  const;
// void  get_el_lq(
//   const uint Level,       // level
//   const int nvars[],       // # of variables to get  <-
//   const int el_nds[],      // # of element nodes for this variable  <-
//   const uint el_conn[],   // connectivity <-
//   const uint offset,      // offset for connectivity <-
//   std::vector<uint>   & el_dof_indices, // element connectivity ->
//   uint  bc_vol[],        // element boundary cond flags ->
//   uint bc_bd[],        // element boundary cond flags ->			 
//   double uold[]           // element node values ->
// )  const;	
 
void  interp_el_sol (
  const double uold_b[],     // node values <-
  const int ivar0,          // init variable  <-		      
  const int nvars,          // # of variables  <-
 const double phi[],        // shape functions  <-
  const int n_shape,        // # of shape functions  <-
  double uold[]            // interpolated function ->
)  const;

void  interp_el_bd_sol (
  const double uold_b[],     // node values <-
  const int sur_tpgly[],     //surface nodes topology <-
  const int el_ndof,     //surface nodes topology <-
  const int ivar0,          // init variable  <-		      
  const int nvars,          // # of variables  <-
  const double phi[],        // shape functions  <-
  const int n_shape,        // # of shape functions  <-
  double uold[]             // interpolated function ->
)  const;  // =======================================



void  interp_el_gdx (
  double uold_b[], // node values <-
  const int ivar0,      // init variable  <-		     
  const int nvars,      // # of variables  <-
//  const double phi[],    // shape functions  <-
  const double dphi[],   // derivatives of the shape functions  <-			     
  const int n_shape,    // # of shape functions  <-
  double uold_dx[]       // interpolated derivatives ->
)  const;    

// void  interp_el_dx (
//   double uold_b[], // node values <-
//   const uint ivar0,      // init variable  <-		     
//   const uint nvars,      // # of variables  <-
//   const double phi[],    // shape functions  <-
//   const double dphi[],   // derivatives of the shape functions  <-
//   const uint n_shape,    // # of shape functions  <-
//   double uold[],         // interpolated function ->
//   double uold_dx[]       // interpolated derivatives ->
// )  const;
 
void  compute_jac (
  const int j,
  const int idim,			    
  double uold_b[], // node values <-
  const int nvars,      // # of variables  <-
  const double phi[],    // shape functions  <-
  const double dphi[],   // derivatives of the shape functions  <-			     
  const int n_shape,    // # of shape functions  <-
  double u_forw[],         // interpolated function ->
  double u_back[],         // interpolated function ->			    
  double u_forw_dx[],       // interpolated derivatives ->
  double u_back_dx[]       // interpolated derivatives ->		   

)  const;
///@}

   // ========================================================================
   ///@{ \name RETURN FUNCTIONS
   // ========================================================================
void  get_el_sol (
  const int ivar0,      // initial variable  <-
  const int nvars,      // # of variables to get  <-
  const int el_nds,     // # of element nodes for this variable  <-
  const int el_conn[],  // connectivity <-
  const int offset,     // offset for connectivity <-
  const int kvar0,      // offset  variable for  uold <-
  double  uold[]            // element node values ->
)  const;
// virtual void  f_aux (double a1[]);


//  void  get_el_q(const uint Level,const uint nvars, const uint el_nds,   
// 		const uint el_conn[], const uint offset,
// 		std::vector<uint>& el_dof_indices,std::vector<uint> & bc_vol,
// 		std::vector<uint>   & bc_bd,std::vector<double> & uold) const;
		
/// Dof , the bc and the solution  vector at the nodes of  an element  
  void  get_el(
    const int Level,// level <-
    const int nvar0, // intial varables number based on NDOF_FEM <-
	       const int nvar,// final varables number based on NDOF_FEM <-
	       const int el_nds,// number of nodes for element
               const int el_conn[], // connectivity <-
	       const int offset, // offset <-
               std::vector<int>  &el_dof_indices, // DOF indices ->
	       int   bc_dofs[][NDOF_FEM],// boudary conditions -> 
	       double uold[]//solution ->
	      )  const; 
	       
/// Dof , the bc and the solution  vector at the nodes of  an element  
/*  void  get_el(const uint Level,const uint el_off0,const uint nvar0,
	       const uint nvar,const uint el_nds,
               const uint el_conn[],const uint offset,
               std::vector<uint>  &el_dof_indices,
	       std::vector<uint>   bc_dofs[],
	       std::vector<double>  & uold)  const;*/	       
void  get_el_dof_bc(
  const int Level,       // level
  const int iel,
//   const int nvars[],       // # of variables to get  <-
  const int el_nds[],      // # of element nodes for this variable  <-
  const int el_conn[],   // connectivity <-
  const int offset,      // offset for connectivity <-
  std::vector<int>   & el_dof_indices, // element connectivity ->
 int bc_vol[],        // element boundary cond flags ->
 int  bc_bd[]        // element boundary cond flags ->			 
)  const ;  // ==============================================================			
///@}

// ========================================================================
   ///@{ \name EXTERNAL FIELDS
   // ========================================================================
  void set_ext_fields(const std::vector<FIELDS> & pbName);
 ///@} 
};

#endif

// #endif
