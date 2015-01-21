#ifndef __mgsolverTBKA_h__
#define __mgsolverTBKA_h__

#include "Equations_conf.h"
// ===================================
#ifdef TBKA_EQUATIONS
// ==================================
#include "MGSolverDA.h"

// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGEquationsSystem;
class MGFEMap;

// =================================================
/// Class for mg energy equation solvers
class MGSolTBKA: public MGSolDA {

private:
  // ==================================================
  /// a) MGSolNS element data
  // mesh
  const int  _offset; ///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
  // constant reference parameters
  const double _dt;    ///< =_mgutils.get_par("dt");
  const double _uref;  ///< =_mgphys.get_par("Uref");
  const double _lref;  ///< =_mgphys.get_par("Lref"); 
  // constant fluid properties 
  const double _rhof;  ///< _mgphys.get_par("rho0");
  const double _muf;   ///< =_mgphys.get_par("mu0");
  // nondimensional numbers
  int _dir;
  double _IRe;        ///< nondimensional conducibility
  double _c_source[2];   ///< nondimensional turb source coeff
  double _c_diss[2];     ///< nondimensional turb diss coeff
  double _chi_k;      ///< cross term coeff
 
  double _kappa_g[2];    ///< reference kappa-omega
//   double _omega_g;    ///< reference omega
  double _y_dist;     ///< distance from the wall
  double _sP;         ///< turbulent tensor modulus
  double _mu_turb;     ///< effective turbulent viscosity 
  double _cross;      ///< SST term for cross coefficient
  double _sigma[2];      ///< diffusion coefficient for k and e/w equations
  
  // -------------------- class field ----------------------------
  // element boundary conditions
  int   _bc_vol[2*NDOF_FEM]; ///< element  b.cond flags
  int   _bc_bd[2*NDOF_FEM]; ///< element  b.cond flags  

  // ------------------ integration -----------------------
  // shape and derivative function at a  gaussian point
 // double _phi_g[3][NDOF_FEM];              ///< shape field (0-1-2 degree)
 // double _dphi_g[3][NDOF_FEM*DIMENSION];   ///< shape derivative  (0-1-2 degree)
  //  fields at gaussian points
  double  _ub_g[3][12];                     ///< external field  (0-1-2 degree)
  double  _ub_dxg[2*DIMENSION];     ///< external field derivative  (0-1-2 degree) 

  // =============================================================
public:

  /// b)   Init MGSolNS functions (constructor,destructor,external fields)
  MGSolTBKA(MGEquationsSystem& mg_equations_map_in,   ///< Constructor
	 int             nvars_in[],    
         std::string eqname_in="KA",
         std::string varname_in="ka"
        );
  ~MGSolTBKA() {}                                  ///< Destructor
//   void set_ext_fields();                        ///< External field set up

  /// c) Read MGSolNS functions
  void bc_read(int bc_gam,int bc_mat,double xp[],int bc_Neu[], 
	       int bc_value[]);                  ///< Read bc
  void ic_read(int bc_gam,int bc_mat,double xp[],double u_value[]);    ///< Read ic

  /// d)  Assemblying MGSolNS Operators
  void GenMatRhs(const double time,              ///< Volume Assemblying
                 const int Level,const int mode);
//   void GenMatRhsB(const double time,             ///< Surface Assemblying
//                   const int Level,const int mode){};
  void MGTimeStep(const double time,             ///< Time-step manager function
                  const int /*iter*/);

  void  f_mu(double val[]);
};
#endif  // define TBK_EQUATIONS
#endif //__mgsolverT_h__
