#ifndef __mgfsisolverns_h__
#define __mgfsisolverns_h__

#include "Equations_conf.h"
// =================================
 #ifdef FSI_EQUATIONS
// =================================
// classe include ---------
#include "MGSolverDA.h"

// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGEquationsMap;
class MGFEMap;



/// Class MGSolFSI:
class MGSolFSI : public MGSolDA {
private:
  // =====================  start data =======================================
  /// a) MGSolFSI element data
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
   int   _bc_vol[NDOF_FEM*(DIMENSION+1)]; ///<  element  b.cond flags (vol int)
   int   _bc_bd[NDOF_FEM*(DIMENSION+1)];  ///<  element  b.cond flags  (bd int)

// 
//   // ------------------ integration -----------------------
//   // shape and derivative function at a  gaussian point
//   //  fields at gaussian points
   double  _ub_g[3][10];                  ///< external field  (0-1-2 degree)
   double  _disp_g[3],_disp_old_g[3],_ale_vel_g[3];                  ///< external field  (0-1-2 degree)
   double  _disp_el[DIMENSION*NDOF_FEM];
   double  _old_solg[8*DIMENSION];
   double  _ub_dxg[DIMENSION*DIMENSION];  ///< external field derivative  (0-1-2 degree)
   double  _old_sol[20*NDOF_FEM];
  // ======================  end data ==========================================
  
  // ======================  start functions =================================== 
public:
       
  /// b)   Init MGSolFSI functions (constructor,destructor,external fields):
  MGSolFSI(                               ///< Constructor
    MGEquationsSystem& mg_equations_map, ///< equation map class (Mesh and parameters)
    int             nvars_in[],          ///< KLQ number of variables
    std::string     eqname_in="FSI0",     ///< base name system
    std::string     varname_in="u"       ///< base name variable
  );
  ~MGSolFSI() {};                        ///< Destructor
         
  void bc_intern_read(
  int bc_gam,
  int mat_gam,
  double xp[],   // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], // normal
  int bc_flag[]  // boundary condition flag
) ; // ===================================
  /// c)          Read MGSolFSI functions:
  /// This function  read bc
  void ic_read(            
    int bc_gam,
    int bc_mat,
    double xp[],            // point coordinates 
    double u_value[]        // point field values
  );   
   /// This function reads ic
  void bc_read(            
   int bc_gam,
   int bc_mat,
    double x[],             // point coordinates 
    int bc_Neu[],           //  bc volume integral flag  
    int bc_bd[]             //  bc surface integral flag 
  ); 

  /// d)  Assemblying MGSolFSI Operators
  /// This function computes volume and surface assemblying
  void GenMatRhs(          
    const double time,        // time
    const int Lev,            // Level
    const  int m              // rhs assembly control
  );
 // Multigrid function ------------------------------------------------------------------
  /// This function is the  time-step manager function 
  void MGTimeStep(       
    const double time,         ///< time
    const int /*iter*/         ///< number max of iterations
  );
  // ====================================================================================
  /// This fucntion compute functional for optimal control
//   void MGFunctional(           
//     const double time,            
//     double starting_distance
//   ) {};
    
    
#ifdef TBK_EQUATIONS // Turbulence only (4-5) --------------------------
  void  f_mu(
    double val[]
  );
#endif
};


#endif // endif NS_EQUATIONS
#endif // endif _mgsolverns_h
