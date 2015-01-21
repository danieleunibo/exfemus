#ifndef __mgsolversmP_h__
#define __mgsolversmP_h__

#include "Equations_conf.h"

// ==================================
#ifdef SMP_EQUATIONS
// ==================================

// This class files -----------------
#include "MGSSMconf.h"
// Local Includes -----------------
#include "MGSolverDA.h"
// Forward declarations ----------
class MGSystem;

// =================================================
/// Class for MG Pressure equation solvers
// =================================================
// config files ------------------
#include "MGFEconf.h"


// =================================================
/// Class for mg energy equation solvers
class MGSolSMP: public MGSolDA {
  
private:
  // -----------------------------------------------
  // auxilyary element data ------------------------
  // mesh
  const uint  _offset;  // = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
  // physics
  const double _dt;     // =_mgutils.get_par("dt");
  const double _rhof;   // =_mgphys.get_par("rho0");
  const double _uref;   // =_mgphys.get_par("Uref");
  const double _lref;   // =_mgphys.get_par("Lref");
  const double _Tref;   // =_mgphys.get_par("Tref");
  const double _cp0;    // =_mgphys.get_par("cp0");
  const double _kappa0; // =_mgphys.get_par("kappa0");
  const double _muf;    // =_mgphys.get_par("mu0");
  
  double       _alpha;  // nondimensional conducibility 
  double       _IPrdl;  // nondimensional Prandl number
  double       _IRe;    // nondimensional Reynolds number
  double       _qheat;  // nondimensional thermal source
  double       _qs;     // nondimensional thermal surface flux
  
  // field --------- -----------------------------------
  // element system field                                   
  double _u_old[NDOF_FEM];      // element old solution
  uint   _bc_dofs[2][NDOF_FEM]; // element  b.cond flags  
  
  // external fields -----------------------------------
  int _nvar_aux;     // index to external field
#ifdef SM_EQUATIONS  // Navier-Stokes system
  int _ns_idx;       // index to external field
#endif
#ifdef TBK_EQUATIONS // Turbulence
  int _tb_idx;       // index to external Turbulence eq
#endif 
  double _ub[(2*DIMENSION+4)*NDOF_FEM]; // element old solution

  // integration -----------------------------------------
  // shapes
  double _phi_g[NDOF_FEM],_dphi_g[NDOF_FEM*DIMENSION];//  shape and derivative function
    double _phi2_g[NDOF_FEM],_dphi2_g[NDOF_FEM*DIMENSION];//  shape and derivative function
  // field at a  gaussian point
  double _uold_g[1],_uold_dxg[NDOF_FEM*DIMENSION];   //  old shape and derivative values
  double _ub_g[2*DIMENSION+4];                       // xxg(DIM)+vel(DIM)+mut

  // -----------------------------------------------------------
public:
  // Constructor-destructor --------------------------
    /// Level constructor 
    MGSolSMP(MGEquationsMap& mg_equations_map_in,
           std::string eqname_in="SMP",   
           std::string varname_in="P"
	  );

   /// Destructor (level structure)
    ~MGSolSMP() {}  
    
    // Setting ---------------------------------------
    /// Boundary conditions
    void bc_read(double xp[],int bc_Neu[], int bc_value[]); 
    /// Initial conditions
    void ic_read(double xp[],double u_value[]);           
    
    // Assemblying ------------------------------------
    /// Volume Assemblying.
    void GenMatRhs(const double time,const uint Level,const int mode); 
    /// Surface Assemblying
    void GenMatRhsB(const double time, const uint Level,const int mode);
    
    // Multigrid function -----------------------------
    /// Timestep function 
    void MGTimeStep(const double time, const uint /*iter*/); 
#ifdef TBK_EQUATIONS
     void  f_mu(double val[]); 
#endif 
};




#endif  // define P_EQUATIONS
#endif //__mgsolverT_h__
