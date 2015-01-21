#ifndef __mgsolverddds_h__
#define __mgsolverddds_h__

#include "Equations_conf.h"

// =================================
 #ifdef DS_EQUATIONS
// =================================

// classe include ---------
#include "MGSolverDA.h"

// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGFEMap;

// =================================================
/// Class for mg displ solvers
class MGSolDS: public MGSolDA {

private:
  // -----------------------------------------------
  // auxilyary element data ------------------------
  // mesh
  const int  _offset;  ///< offset mesh nodes (top level)
  // physics
  const double _dt;     ///< time step
  const double _rhof;   ///< fluid density
  const double _uref;   ///< ref velocity
  const double _lref;   ///< ref length
  const double _Tref;   ///< ref temp
  const double _rhos;   ///< solid density
  int    _dir;          ///< direction

  // -------------------- class field ----------------------------
  // element boundary conditions
 int   _bc_vol[NDOF_FEM]; ///< element  b.cond flags
 int   _bc_bd[NDOF_FEM]; ///< element  b.cond flags  
 
  // ------------------ integration -----------------------
  // shape and derivative function at a  gaussian point
 // double _phi_g[3][NDOF_FEM];              ///< shape field (0-1-2 degree)
 // double _dphi_g[3][NDOF_FEM*DIMENSION];   ///< shape derivative  (0-1-2 degree)
  //  fields at gaussian points
 double  _ub_g[3][12];                     ///< external field  (0-1-2 degree)
 double  _ub_dxg[DIMENSION*DIMENSION];     ///< external field derivative  (0-1-2 degree) 
 double  _ub_old[12*NDOF_FEM];   
  // =============================================================

  // -----------------------------------------------------------
public:
  // Constructor-destructor --------------------------
  /// Level constructor
  MGSolDS(MGEquationsSystem& mg_equations_map_in,
	  const int vars_in[],
            std::string eqname_in="DS",
            std::string varname_in="d"
           );

  /// Destructor (level structure)
  ~MGSolDS() {}

  
   // Setting ---------------------------------------
  /// Boundary conditions
  void bc_read(int bc_gam,int mat_gam, double xp[],int bc_Neu[], int bc_value[]);
//  void bc_intern_read(double x[],int bc_Neu[],int u[]); ///< Read ic.
  /// Initial conditions
  void ic_read(int bc_gam,int mat_gam, double xp[],double u_value[]);

  void bc_intern_read (int bc_gam,int mat_gam,double xp[], int normal[],int bc[]);
   
  // Assemblying ------------------------------------
  /// Volume Assemblying.
  void GenMatRhs(const double time,const int Level,const int mode);
  
  void MoveMesh(const int Level);
//   /// Surface Assemblying
//   void GenMatRhsB(const double time, const int Level,const int mode);

  // Multigrid function -----------------------------
  /// Timestep function
  void MGTimeStep(const double time, const int /*iter*/);
  
  inline void print_ext_data(double /* table_data*/[]) { };

};

#endif  //  #endif DS_EQUATIONS
#endif

// ====================================================================================