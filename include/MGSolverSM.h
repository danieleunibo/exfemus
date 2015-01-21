#ifndef __mgsolversm_h__
#define __mgsolversm_h__

#include "Equations_conf.h"
// =================================
 #ifdef SM_EQUATIONS
// =================================

// classe include ---------
#include "MGSolverDA.h"

// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGEquationsMap;
class MGFEMap;



/// Class MGSolSM:
class MGSolSM : public MGSolDA {
private:
  // =====================  start data =======================================
  /// a) MGSolSM element data
  // mesh
  const int  _offset; ///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
  int   _dir;
  
  // constant reference parameters
  const double _dt;    ///< =_mgutils.get_par("dt");
  const double _uref;  ///< =_mgphys.get_par("Uref");
  const double _lref;  ///< =_mgphys.get_par("Lref");
  double _dirg[3];     ///< =gravity
  // constant fluid properties
  const double _rhof;  ///< =_mgphys.get_par("rho0");
  const double _muf;   ///< =_mgphys.get_par("mu0");
   // constant solid propertie
  const double _rhos;  ///< =_mgphys.get_par("mu0");
  const double _ni;    ///< =_mgphys.get_par("mu0");
  const double _Emod;  ///< =_mgphys.get_par("mu0");
  const double _hs;    ///< =_mgphys.get_par("mu0");
  double _lambda;
  double _mus;
  
  // nondimensional numbers
  double _IRe;         ///< =Reynolds number
  double _IFr;         ///< =Froud number

  double _kappa_g[2];    ///< reference kappa
//   double _omega_g;    ///< reference omega
  double _y_dist;     ///< distance from the wall
  double _sP;         ///< turbulent tensor modulus
  double _mu_turb;     ///< effective turbulent viscosity 
  
//   // -------------------- class field ----------------------------
//   // element class field  (u,v,w,p)
// //   double _u1_old[NDOF_FEM*(DIMENSION+SM_EQUATIONS%2)]; ///<  old solution (element)
   int   _bc_vol[NDOF_FEM*(DIMENSION)]; ///<  element  b.cond flags (vol int)
   int   _bc_bd[NDOF_FEM*(DIMENSION)];  ///<  element  b.cond flags  (bd int)
// 
//   // ------------------ integration -----------------------
//   // shape and derivative function at a  gaussian point
//   double _phi_g[3][NDOF_FEM];            ///< shape field (0-1-2 degree)
//   double _dphi_g[3][NDOF_FEM*DIMENSION]; ///< shape derivative  (0-1-2 degree)
//   //  fields at gaussian points
   double  _ub_g[3][10];                  ///< external field  (0-1-2 degree)
   double  _ub_dxg[10];  ///< external field derivative  (0-1-2 degree)
  
  // ======================  end data ==========================================
  
  // ======================  start functions =================================== 
public:

  /// b)   Init MGSolSM functions (constructor,destructor,external fields)
  MGSolSM(MGEquationsMap& mg_equations_map,  ///< Constructor
          std::string eqname_in="SM0",std::string varname_in="u");
  ~MGSolSM() {};                            ///< Destructor

  /// c) Read MGSolSM functions
  void ic_read(double xp[],double u_value[]);    ///< Read bc
  void bc_read(double x[],int bc_Neu[],int u[]); ///< Read ic

  /// d)  Assemblying MGSolSM Operators
  void GenMatRhs(const double time,    ///< Volume Assemblying
                 const int Lev,const  int m);
 // Multigrid function -----------------------------
    /// Timestep function 
  void MGTimeStep(const double time,  ///< Time-step manager function
                  const int /*iter*/);

#ifdef TBK_EQUATIONS // Turbulence only (4-5) --------------------------
  void  f_mu(double val[]);
#endif
};


// =======================================================================
// =======================================================================
// Projection method (Poisson Pressure equation)
// =======================================================================
#if (SMP_EQUATIONS==1)
// =======================================================================

// =================================================
/// Class for MG Pressure equation solvers
// =================================================
class MGSolSMP: public MGSolDA {
// =================================================  
private:
    // =======================================================================
  /// a) MGSolSMP element data
  // mesh
  const int   _offset;  // = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
 
  // physics
  const double _dt;                          ///< =_mgutils.get_par("dt");
  const double _rhof;                        ///< =_mgphys.get_par("rho0");
  const double _uref;                        ///< =_mgphys.get_par("Uref");
  const double _lref;                        ///< =_mgphys.get_par("Lref");

//   // -------------------- class field ----------------------------
//   // element boundary conditions
   int    _bc_vol[NDOF_FEM];///<  element  b.cond flags (vol int)
   int    _bc_bd[NDOF_FEM]; ///<  element  b.cond flags  (bd int)
// 
//   // ------------------ integration -----------------------
//   // shape and derivative function at a  gaussian point
//   double _phi_g[3][NDOF_FEM];              ///< shape field (0-1-2 degree)
//   double _dphi_g[3][NDOF_FEM*DIMENSION];   ///< shape derivative  (0-1-2 degree)
//   //  fields at gaussian points
   double  _ub_g[3][10];                     ///< external field  (0-1-2 degree)
   double  _ub_dxg[DIMENSION*DIMENSION];     ///< external field derivative  (0-1-2 degree)
  // -----------------------------------------------------------
public:
  // Constructor-destructor --------------------------
    /// Level constructor 
    MGSolSMP(MGEquationsMap& mg_equations_map_in,
           std::string eqname_in="SM1P",   
           std::string varname_in="p"
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
    void GenMatRhs(const double time,const int Level,const int mode); 
    
    // Multigrid function -----------------------------
    /// Timestep function 
    void MGTimeStep(const double time, const int  /*iter*/); 

};

#endif  // endif SM_EQUATIONS==0
#endif // endif SM_EQUATIONS
#endif // endif _mgsolverns_h
