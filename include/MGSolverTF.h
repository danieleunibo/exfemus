#ifndef __mgsolver2f_h__
#define __mgsolver2f_h__

#include "Equations_conf.h"
// =================================
#ifdef TF_EQUATIONS
// =================================
// classe include ---------
#include "MGSolverDA.h"

// Forwarded classes ------
class MGUtils;
class MGSystem;
class MGEquationsMap;
class MGFEMap;



/// Class MGSolTF:
class MGSolTF : public MGSolDA {
private:
  // =====================  start data =======================================
  /// a) MGSolTF element data
  // mesh
  const int  _offset; ///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
  int   _dir;         ///< x-direction for segrgated mode
  
  // constant reference parameters
  const double _dt;    ///< =_mgutils.get_par("dt");
  const double _uref;  ///< =_mgphys.get_par("Uref");
  const double _lref;  ///< =_mgphys.get_par("Lref");
  double _dirg[3];     ///< =gravity
  int _phase;	       ///< phase 0=liquid, 1=gas
  // constant fluid properties
  const double _rhof;  ///< =_mgphys.get_par("rho0");
  const double _muf;   ///< =_mgphys.get_par("mu0");
  const double _rhog;  ///< =_mgphys.get_par("rhog");
  const double _mug;   ///< =_mgphys.get_par("mug");
  // nondimensional numbers
  double _IRe;         ///< =Reynolds number
  double _IFr;         ///< =Froud number

  double _kappa_g[2];  ///< reference kappa
//   double _omega_g;    ///< reference omega
  double _y_dist;     ///< distance from the wall
  double _sP;         ///< turbulent tensor modulus
  double _mu_turb;    ///< effective turbulent viscosity 
  
//   // -------------------- class field ----------------------------
//   // element class field  (u,v,w,p)
   int   _bc_vol[NDOF_FEM*(DIMENSION+1)]; ///<  element  b.cond flags (vol int)
   int   _bc_bd[NDOF_FEM*(DIMENSION+1)];  ///<  element  b.cond flags  (bd int)
// 
//   // ------------------ integration -----------------------
//   // shape and derivative function at a  gaussian point
//   //  fields at gaussian points
   double  _ub_g[3][14];                  ///< external field  (0-1-2 degree)
   double  _ub_dxg[DIMENSION*DIMENSION];  ///< external field derivative  (0-1-2 degree)
  
  // ======================  end data ==========================================
  
  // ======================  start functions =================================== 
public:

  /// b)   Init MGSolTF functions (constructor,destructor,external fields):
  MGSolTF(                               ///< Constructor
    MGEquationsMap& mg_equations_map,    // equation map class (Mesh and parameters)
    int             nvars_in[],          // KLQ number of variables
    std::string     eqname_in,           // base name system
    std::string     varname_in           // base name variable
  );
  ~MGSolTF() {};                        ///< Destructor
         
  /// c)          Read MGSol2F functions:
  void ic_read(             ///< Read ic
    double xp[],            // point coordinates 
    double u_value[]        // point field values
  );   
  void bc_read(             ///< Read bc
   int bc_gam,              // gambit interface
    double x[],             // point coordinates 
    int bc_Neu[],           //  bc volume integral flag  
    int bc_bd[]             //  bc surface integral flag 
  ); 

  /// d)  Assemblying MGSol2F Operators
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
  
  
  void init_subelnod(int sub_el_nod[],           ///< Initialize topology of sub elements
		     int int_sides_nodes[],
		     int ext_sides_nodes[],
		     int mat_index[],
		     int boun_index[]
		    );
  void alpha_check(
  std::vector<NumericVectorM *> x
);

};

// =================================================
/// Class for mg energy equation solvers
class MGSolPTF: public MGSolDA {

private:
  // ==================================================
  /// a) MGSolNS element data
  // mesh
  const int  _offset;  ///< = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
  // constant reference parameters
  double _dt;     ///< =_mgutils.get_par("dt");
  // density
  const double _rhof;  ///< =_mgphys.get_par("rho0");
  const double _rhog;  ///< =_mgphys.get_par("rhog");

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
  MGSolPTF(MGEquationsMap& mg_equations_map_in,   ///< Constructor
	 int             nvars_in[],         // KLQ number of variables
         std::string     eqname_in,
         std::string     varname_in
        );
  ~MGSolPTF() {};                                  ///< Destructor
//   void set_ext_fields();                        ///< External field set up

  /// c) Read MGSolNS functions
  void bc_read(int bc_gam, double xp[],int bc_Neu[], 
	       int bc_value[]);                  ///< Read bc
  void ic_read(double xp[],double u_value[]);    ///< Read ic

  /// d)  Assemblying MGSolNS Operators
  void GenMatRhs(const double time,              ///< Volume Assemblying
                 const int Level,const int mode);
//   void GenMatRhsB(const double time,             ///< Surface Assemblying
//                   const int Level,const int mode){};
  void MGTimeStep(const double time,             ///< Time-step manager function
                  const int /*iter*/);

};

#endif // endif 2F_EQUATIONS
#endif // endif _mgsolverns_h
