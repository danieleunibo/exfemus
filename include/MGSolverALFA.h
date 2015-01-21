#ifndef __mgsolverALFA_h__
#define __mgsolverALFA_h__

#include "Equations_conf.h"
// ===================================
#ifdef ALFA_EQUATIONS
// ==================================
// ===================================
// #ifdef ALFA_EQUATIONS_LEGENDRE
// ==================================

// config files ------------------
#include "MGFE_conf.h"

// class files ---------------------
#include "MGSclass_conf.h"

// Local Includes -----------------
#include "MGSolverDA.h"

// Forward declarations ----------
class MGSystem;



// =================================================
/// Class for mg energy equation solvers
class MGSolALFA: public MGSolDA {

private:
  // ==================================================
  /// a) MGSolNS element data
  // mesh
  const int  _offset;  ///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
  // constant reference parameters
  double _dt;     ///< =_mgutils.get_par("dt");
  int _step_levset; ///< =_mgutils.get_par("step") or calculated;
  int _mid_step; ///< integer to discriminate for split advection

  // -------------------- class field ----------------------------
  // element boundary conditions
 int   _bc_vol[NDOF_FEM]; ///< element  b.cond flags
 int   _bc_bd[NDOF_FEM]; ///< element  b.cond flags  

  // ------------------ integration -----------------------
  // shape and derivative function at a  gaussian point
 // double _phi_g[3][NDOF_FEM];              ///< shape field (0-1-2 degree)
 // double _dphi_g[3][NDOF_FEM*DIMENSION];   ///< shape derivative  (0-1-2 degree)
  //  fields at gaussian points
 double  _ub_g[3][14];             ///< external field  (0-1-2 degree)
 double  _ub_dxg[2*DIMENSION];     ///< external field derivative  (0-1-2 degree) 

  // =============================================================
public:

  /// b)   Init MGSolNS functions (constructor,destructor,external fields)
  MGSolALFA(MGEquationsMap& mg_equations_map_in,   ///< Constructor
	 int             nvars_in[],         // KLQ number of variables
         std::string eqname_in,
         std::string varname_in
        );
  ~MGSolALFA() {}                                  ///< Destructor
//   void set_ext_fields();                        ///< External field set up

  /// c) Read MGSolNS functions
  void bc_read(int bc_gam, int mat_gam , double xp[],int bc_Neu[], 
	       int bc_value[]);                  ///< Read bc
  void ic_read(int bc_gam, int mat_gam ,double xp[],double u_value[]);    ///< Read ic

  /// d)  Assemblying MGSolNS Operators
  void GenMatRhs(const double time,              ///< Volume Assemblying
                 const int Level,const int mode);
//   void GenMatRhsB(const double time,             ///< Surface Assemblying
//                   const int Level,const int mode){};
  void MGTimeStep(const double time,             ///< Time-step manager function
                  const int /*iter*/);
  void init_subelnod(int sub_el_nod[],           ///< Initialize topology of sub elements
		     int int_sides_nodes[],
		     int ext_sides_nodes[],
		     int mat_index[],
		     int boun_index[]
		    );
  
  void level_set(std::vector<NumericVectorM *> x);

};
#endif  // define ALFA_EQUATIONS
#endif //__mgsolveralfa_h__
