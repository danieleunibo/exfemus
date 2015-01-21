#ifndef __mgsolvernsA_h__
#define __mgsolvernsA_h__

#include "Equations_conf.h"
// =================================
 #ifdef NSA_EQUATIONS
// =================================
// classe include ---------
#include "MGSolverDA.h"

// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGEquationsSystem;
class MGFEMap;



/// Class MGSolNS:
class MGSolNSA : public MGSolDA {
private:
  // =====================  start data =======================================
  /// a) MGSolNS element data
  // mesh
  const int  _offset; ///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
  int   _dir;         ///< x-direction for segrgated mode
  
  // constant reference parameters
  const double _dt;    ///< =_mgutils.get_par("dt");
  const double _uref;  ///< =_mgphys.get_par("Uref");
  const double _lref;  ///< =_mgphys.get_par("Lref");
  double _dirg[3];     ///< =gravity
  // constant fluid properties
  const double _rhof;  ///< =_mgphys.get_par("rho0");
  const double _muf;   ///< =_mgphys.get_par("mu0");
  // nondimensional numbers
  double _IRe;         ///< =Reynolds number
  double _IFr;         ///< =Froud number

  double _kappa_g[2];  ///< reference kappa
//   double _omega_g;    ///< reference omega
  double _y_dist;     ///< distance from the wall
  double _sP;         ///< turbulent tensor modulus
  double _mu_turb;    ///< effective turbulent viscosity 
  double _beta;       ///< parameter for control 
  
//   // -------------------- class field ----------------------------
//   // element class field  (u,v,w,p)
   int   _bc_vol[NDOF_FEM*(DIMENSION+1)]; ///<  element  b.cond flags (vol int)
   int   _bc_bd[NDOF_FEM*(DIMENSION+1)];  ///<  element  b.cond flags  (bd int)
// 
//   // ------------------ integration -----------------------
//   // shape and derivative function at a  gaussian point
//   //  fields at gaussian points
   double  _ub_g[3][10];                  ///< external field  (0-1-2 degree)
   double  _ub_dxg[DIMENSION*DIMENSION];  ///< external field derivative  (0-1-2 degree)
  
  // ======================  end data ==========================================
  
  // ======================  start functions =================================== 
public:

  /// b)   Init MGSolNS functions (constructor,destructor,external fields):
  MGSolNSA(                               ///< Constructor
    MGEquationsSystem& mg_equations_map,    // equation map class (Mesh and parameters)
    int             nvars_in[],         // KLQ number of variables
    std::string     eqname_in="NSA0",     // base name system
    std::string     varname_in="ua"       // base name variable
  );
  ~MGSolNSA() {};                        ///< Destructor
         
  /// c)          Read MGSolNS functions:
  void ic_read(             ///< Read bc
    int bc_gam,
    int bc_mat,
    double xp[],            // point coordinates 
    double u_value[]        // point field values
  );   
  void bc_read(             ///< Read ic
   int bc_gam,
   int bc_mat,
    double x[],             // point coordinates 
    int bc_Neu[],           //  bc volume integral flag  
    int bc_bd[]             //  bc surface integral flag 
  ); 

  /// d)  Assemblying MGSolNS Operators
  void GenMatRhs(          ///< Volume and surface Assemblying
    const double time,        // time
    const int Lev,            // Level
    const  int m              // rhs assembly control
  );
 // Multigrid function -----------------------------
 
  void MGTimeStep(        ///< Time-step manager function
    const double time,         // time
    const int /*iter*/         // number max of iterations
  );

#ifdef TBK_EQUATIONS // Turbulence only (4-5) --------------------------
  void  f_mu(double val[]);
#endif
};


#endif // endif NS_EQUATIONS
#endif // endif _mgsolverns_h
