#ifndef __mgsolverT_h__
#define __mgsolverT_h__

#include "Equations_conf.h"
// ===================================
#ifdef T_EQUATIONS
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
class MGSolTC: public MGSolDA {

private:
  // ==================================================
  /// a) MGSolNS element data
  // mesh
  const int  _offset;  ///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
  // constant reference parameters
  double _dt;     ///< =_mgutils.get_par("dt");
  const double _uref;   ///< =_mgphys.get_par("Uref");
  const double _lref;   ///< =_mgphys.get_par("Lref");
  const double _Tref;   ///< =_mgphys.get_par("Tref");
  // constant fluid properties
  const double _rhof;   ///< =_mgphys.get_par("rho0");
  const double _muf;    ///< =_mgphys.get_par("mu0");
  const double _cp0;    ///< =_mgphys.get_par("cp0");
  const double _kappa0; ///< =_mgphys.get_par("kappa0");
//   const double _h_conv; ///< =_mgphys.get_par("hconv");
//   const double _T_inf;  ///< =_mgphys.get_par("T_inf");
  // nondimensional numbers
  double       _alpha;  ///< nondimensional conducibility
  double       _IPrdl;  ///< nondimensional Prandl number
  double     _alpha_turb;  ///< nondimensional conducibility
  double     _IPrdl_turb;  ///< nondimensional Prandl number
  double       _IRe;    ///< nondimensional Reynolds number
//   double       _Nusselt;///< nondimensional Nusselt number
   double       _qheat;  ///< nondimensional thermal source
   double       _qs;     ///< nondimensional thermal surface flux

  // turbulence
  double _kappa_g[2];    ///< reference kappa
  double _kappaT_g[2];    ///< reference omega
  double _y_dist;     ///< distance from the wall
  double _sP;         ///< turbulent tensor modulus
  double _nut_ratio;     ///< effective turbulent viscosity 
  
  // -------------------- class field ----------------------------
  // element boundary conditions
 int   _bc_vol[NDOF_FEM]; ///< element  b.cond flags
 int   _bc_bd[NDOF_FEM]; ///< element  b.cond flags  

  // ------------------ integration -----------------------
  // shape and derivative function at a  gaussian point
 // double _phi_g[3][NDOF_FEM];              ///< shape field (0-1-2 degree)
 // double _dphi_g[3][NDOF_FEM*DIMENSION];   ///< shape derivative  (0-1-2 degree)
  //  fields at gaussian points
 double  _ub_g[3][13];             ///< external field  (0-1-2 degree)
 double  _ub_dxg[2*DIMENSION];     ///< external field derivative  (0-1-2 degree) 

  // =============================================================
public:

  /// b)   Init MGSolNS functions (constructor,destructor,external fields)
  MGSolTC(MGEquationsMap& mg_equations_map_in,   ///< Constructor
	 int             nvars_in[],         // KLQ number of variables
         std::string eqname_in="T",
         std::string varname_in="T"
        );
  ~MGSolTC() {}                                  ///< Destructor
//   void set_ext_fields();                        ///< External field set up

  /// c) Read MGSolNS functions
  void bc_read(int bc_gam, double xp[],int bc_Neu[], 
	       int bc_value[]);                  ///< Read bc
  void ic_read(int bc_gam, double xp[],double u_value[]);    ///< Read ic

  /// d)  Assemblying MGSolNS Operators
  void GenMatRhs(const double time,              ///< Volume Assemblying
                 const int Level,const int mode);
//   void GenMatRhsB(const double time,             ///< Surface Assemblying
//                   const int Level,const int mode){};
  void MGTimeStep(const double time,             ///< Time-step manager function
                  const int /*iter*/);
#ifdef TBK_EQUATIONS
  void  f_mu(double val[]);
#endif
};





// =================================================
/// Class for mg energy equation solvers
class MGSolTC: public MGSolDA {

private:
  // ==================================================
  /// a) MGSolNS element data
  // mesh
  const int  _offset;  ///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
  // constant reference parameters
  double _dt;     ///< =_mgutils.get_par("dt");
  const double _uref;   ///< =_mgphys.get_par("Uref");
  const double _lref;   ///< =_mgphys.get_par("Lref");
  const double _Tref;   ///< =_mgphys.get_par("Tref");
  // constant fluid properties
  const double _rhof;   ///< =_mgphys.get_par("rho0");
  const double _muf;    ///< =_mgphys.get_par("mu0");
  const double _cp0;    ///< =_mgphys.get_par("cp0");
  const double _kappa0; ///< =_mgphys.get_par("kappa0");
//   const double _h_conv; ///< =_mgphys.get_par("hconv");
//   const double _T_inf;  ///< =_mgphys.get_par("T_inf");
  // nondimensional numbers
  double       _alpha;  ///< nondimensional conducibility
  double       _IPrdl;  ///< nondimensional Prandl number
  double     _alpha_turb;  ///< nondimensional conducibility
  double     _IPrdl_turb;  ///< nondimensional Prandl number
  double       _IRe;    ///< nondimensional Reynolds number
//   double       _Nusselt;///< nondimensional Nusselt number
   double       _qheat;  ///< nondimensional thermal source
   double       _qs;     ///< nondimensional thermal surface flux

  // turbulence
  double _kappa_g[2];    ///< reference kappa
  double _kappaT_g[2];    ///< reference omega
  double _y_dist;     ///< distance from the wall
  double _sP;         ///< turbulent tensor modulus
  double _nut_ratio;     ///< effective turbulent viscosity 
  
  // -------------------- class field ----------------------------
  // element boundary conditions
 int   _bc_vol[NDOF_FEM]; ///< element  b.cond flags
 int   _bc_bd[NDOF_FEM]; ///< element  b.cond flags  

  // ------------------ integration -----------------------
  // shape and derivative function at a  gaussian point
 // double _phi_g[3][NDOF_FEM];              ///< shape field (0-1-2 degree)
 // double _dphi_g[3][NDOF_FEM*DIMENSION];   ///< shape derivative  (0-1-2 degree)
  //  fields at gaussian points
 double  _ub_g[3][13];             ///< external field  (0-1-2 degree)
 double  _ub_dxg[2*DIMENSION];     ///< external field derivative  (0-1-2 degree) 

  // =============================================================
public:

  /// b)   Init MGSolNS functions (constructor,destructor,external fields)
  MGSolTC(MGEquationsMap& mg_equations_map_in,   ///< Constructor
	 int             nvars_in[],         // KLQ number of variables
         std::string eqname_in="T",
         std::string varname_in="T"
        );
  ~MGSolTC() {}                                  ///< Destructor
//   void set_ext_fields();                        ///< External field set up

  /// c) Read MGSolNS functions
  void bc_read(int bc_gam, double xp[],int bc_Neu[], 
	       int bc_value[]);                  ///< Read bc
  void ic_read(int bc_gam, double xp[],double u_value[]);    ///< Read ic

  /// d)  Assemblying MGSolNS Operators
  void GenMatRhs(const double time,              ///< Volume Assemblying
                 const int Level,const int mode);
//   void GenMatRhsB(const double time,             ///< Surface Assemblying
//                   const int Level,const int mode){};
  void MGTimeStep(const double time,             ///< Time-step manager function
                  const int /*iter*/);
#ifdef TBK_EQUATIONS
  void  f_mu(double val[]);
#endif
};



#endif  // define T_EQUATIONS
#endif //__mgsolverT_h__
