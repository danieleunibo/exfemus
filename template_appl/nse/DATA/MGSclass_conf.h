#include <Equations_conf.h>



// =============================================
///    NS_EQUATIONS   MGSNSconf.h
#ifdef NS_EQUATIONS
// =============================================

#ifndef __mgsnsc_h__
#define __mgsnsc_h__


//  ---------------------
// MULTIGRID PARAMETERS

// #define SOLVERNS VANKANSM
#define SOLVERNS GMRESM
// #define SOLVERNS BICGSTABM
#define SOLVERNSP CGM
// #define SOLVERNSP GMRESM


#define NLGEOM  1

//  ---------------------------------
// 3D NAVIER-STOKES ADVECTION term
//
// A) Stokes flow    ADVPIC_NS 0. ADVNEW_NS=0
// B)  Navier-Stokes ADVPIC_NS 1. ADVNEW_NS={0,1}
#define ADVPIC_NS 0.
//  Navier-Stokes  nonlinear iterations
// A) Picard iteration  ADVPIC_NS=1, ADVNEW_NS=0 
// B) Newton iteration  ADVPIC_NS=1, ADVNEW_NS=1
#define ADVNEW_NS 0.


// -----------  stabilization ------------------
// --------------------------------------------
//  Navier-Stokes  stab=STAB_NS*0.5*(div u,u.v)
#define STAB_NS 0.
//  Navier-Stokes  compressibility=KOMP_NS*dp/dt
#define KOMP_NS 1.e-20
//  Navier-Stokes  compressibility  upwind=Re+UP_WIND_NS*v^2
#define UP_WIND_NS (0.1)
//  Navier-Stokes  penalty =LAMBDA grad div u (solid)
#define LAMBDA (0.)
// Navier-Stokes supg
#define SUPG_NS (0.)
// Navier-Stokes  c -> antisymmetric (0.5) 
#define ADV_ASYM 0.





// Crank-Nicolson first order 0. 2nd order 0.5
#define CN_TIME 0.

//  boundary integral
// #define P_0 (0.)

// turbulent 
#define MU_LOW (1.e-12)
#define MU_TOP (1.e+12)


#define SQCMU (0.3)
#define YPLUS (1.)
#define ALPHA0 (1.)


#define CMU (0.09)



#endif


#endif // ============== end NS_equations ===================




// =============================================
///   T_EQUATIONS    MGSTconf.h
#ifdef T_EQUATIONS // 3D-2D Energy equation
// =============================================


#ifndef  __mgstconf_h__
#define  __mgstconf_h__

// turbulence  ----------------------
#define ADVE 1.

// Turbulence Prandl number 
#define PRT (0.85)

// #define SOLVERT VANKATM -----------------
#define SOLVERT GMRESM



// temperature lows -------------------------
 #define  CONST 1
// #define densityT(x) (1.)
// #define cipiT(x) (1.)
// #define kappa(x)  (1.)
#define  UP_WIND_T (1.)
// ------------------------
// #define  LINEAR 1 quad (2)
#define LQ_T (2)

// boundary flux
// #define Q0 (0.)

#endif

#endif   // ======== end Energy ===========



// =============================================
///    DA_EQUATIONS  MGSDAconf.h
#ifdef DA_EQUATIONS 
// =============================================
#ifndef __mgsnsc_h__
#define __mgsnsc_h__



//  -----------------------------
// NAVIER-STOKES ADVECTION term
// --------------------------------
// A) Stokes flow   ADV 0. B)  Navier-Stokes ADV 1.
#define ADV 1.
//  Navier-Stokes  A) Picard iteration  ADV1 (0.) B) Newton ADV 1.
#define ADV1 0.
//  Navier-Stokes  stab  +0.5*(div u,u.v)
#define STAB 1.
// 


#endif

#endif  // ============= end DA_EQUATIONS  ============





// =============================================
///    FSI_EQUATIONS  MGSFSIconf.h
#ifdef FSI_EQUATIONS
// =============================================

#ifndef __mgsnscfsi_h__
#define __mgsnscfsi_h__


//  
// MULTIGRID PARAMETERS
// -------------------------------
// Navier-Stokes solver type
#define SOLVER_FSI GMRESM  // options -> GMRESM BICGSTABM

// Pressure solver type (projection method only)
#define SOLVER_FSIP CGM  // options -> GMRESM CGM BICGSTABM

 
// ------------------------
//   SOLID MODEL
// ----------------------------
 
// geometric non-linearity
 #define NL_GEOM  (1)
 
// define penalty  (only with FSIP_EQUATIONS==1) 
#define  PENALTY_FSI (100.)

 
 
 // ---------------------------
//   SOLID-FLUID REGIONS
// -----------------------------
//  #define SOLID 0
//  #define STIFF 10
 

 
//  --------------------------
// 3D NAVIER-STOKES ADVECTION term
// ---------------------------------
// A) Stokes flow    ADVPIC_SM 0. ADVNEW_SM=0
// B)  Navier-Stokes ADVPIC_SM 1. ADVNEW_SM={0,1}
#define ADVPIC_FSI 0.
//  Navier-Stokes  nonlinear iterations
//  Navier-Stokes  nonlinear iterations
// A) Picard iteration  ADVPIC_SM=1, ADVNEW_SM=0 
// B) Newton iteration  ADVPIC_SM=1, ADVNEW_SM=1
#define ADVNEW_FSI 0.


// -----------  stabilization ------------------
// --------------------------------------------
//  Navier-Stokes  stab=STAB_SM*0.5*(div u,u.v)
// #define STAB_FSI 0.
//  Navier-Stokes  compressibility=KOMP_SM*dp/dt
#define KOMP_FSI 1.e-20
//  Navier-Stokes  compressibility  upwind=Re+UP_WIND_SM*v^2
#define UP_WIND_FSI (0.)
//  Navier-Stokes  penalty =LAMBDA grad div u (solid)
// #define LAMBDA (2100)
// Navier-Stokes supg
// #define SUPG_FSI (0.)
// Navier-Stokes  c -> antisymmetric (0.5) 
// #define ADV_ASYM 0.


// #define LQ (2)


// Crank-Nicolson first order 0. 2nd order 0.5
// #define CN_TIME 0.

//  boundary integral
// #define P_0 (0.)

// turbulent 
// #define MU_LOW (1.e-12)
// #define MU_TOP (1.e+12)


// #define SQCMU (0.3)
// #define YPLUS (1.)
// #define ALPHA0 (1.)



// P solver ------------------------------------

// #define SOLVERT VANKATM ---------------------
// #define SOLVERFSIP GMRESM


// temperature lows -----------------------------
//  #define  CONST 1
// // #define densityT(x) (1.)
// // #define cipiT(x) (1.)
// // #define kappa(x)  (1.)
// 
// // quadratic QL=2 linear QL=1 -------------------
// #define LQ (2)
// #define LQ_P (1)
// // boundary pressure --------------------------
// #define P_BD (1.)


#endif
#endif  // ============= end FSI_EQUATIONS  ============

// =============================================
///   DS_EQUATIONS   MGSDSconf.h
#ifdef DS_EQUATIONS
// =============================================


#endif // ============= end DS_EQUATIONS  ============



// =============================================
///   TBK_EQUATIONS   MGSTBKconf.h
#ifdef TBK_EQUATIONS
// =============================================
#ifndef  __mgstbkconf_h__
#define  __mgstbkconf_h__

// turbulence  -----------------------------
#define ADVE 1.

// Turbulence Prandl number 
#define PRT (0.85)

// #define SOLVERT VANKATM-----------------------------
#define SOLVER_TBK GMRESM


// temperature lows-----------------------------
 #define  CONST 1
// #define densityT(x) (1.)
// #define cipiT(x) (1.)
// #define kappa(x)  (1.)

// #define  LINEAR 1 quad (2)
//  #define LQ_TB  (2)

#define LES (100.)
// boundary flux
#define Q0 (1.)

// turbulence model k-e=0 k-t=1 k-w=2
#define KEWT_MODEL (2)

#define MU_LOW (1.e-12)
#define MU_TOP (1.e+12)

#define SUPG_K (0.)
#define SUPG_E (0.)
#define SUPG_W (0.)
#define UP_WIND_W (1.)
#define UP_WIND_K (1.)

//boundary
#define Y_WALL (0.00005)
#define SQCMU (0.3)
#define KAPPA_VK (0.4)

// kappa - epsilon model -----------------------------
#if ((TBK_EQUATIONS/2)==1)
// kappa omega constant
#define SIGMAO (2.)
#define SIGMAK (2.)
#define C2O (1.83)
#define C1O (1.55)
#define CMU (0.09)
#define BETA (0.09)
#define BETASTAR (0.09)
#endif

// kappa -omega model ----------------------------------------
#if ((TBK_EQUATIONS/2)==2)

//    #define SST (1)
   #define A1 (1000.)

//   #define LOWRE (1)

// #define ALFA_W (0.555555555556)
#define ALFA_W (0.52)
// #define BETA_W (0.075)
#define BETA_W (0.072)
#define SIGMA_W (.5)

#define BETASTAR (0.09)
#define SIGMA_K (.5)

#endif
#endif
#endif  // ============= end MGSTB_EQUATIONS  ============




// =============================================
///   MGSTTBK_EQUATIONS    MGSTTBKconf.h
#ifdef TTBK_EQUATIONS
// =============================================
#ifndef  __mgsttbkconf_h__
#define  __mgsttbkconf_h__

// turbulence  ==========================
#define ADVE 1.

// Turbulence Prandl number 
#define PRT (0.85)

// #define SOLVERT VANKATM =======================
#define SOLVERTBK GMRESM


// temperature lows ================================
 #define  CONST 1
// #define densityT(x) (1.)
// #define cipiT(x) (1.)
// #define kappa(x)  (1.)

// #define  LINEAR 1 quad (2)
//  #define LQ_TB  (2)

#define LES (100.)
// boundary flux
#define Q0 (1.)

// turbulence model k-e=0 k-t=1 k-w=2
#define KEWT_MODEL (2)

#define MU_LOW (1.e-12)
#define MU_TOP (1.e+12)

#define SUPG_K (0.)
#define SUPG_E (0.)
#define SUPG_W (0.)
#define UP_WIND_TW (1.)
#define UP_WIND_TK (1.)

//boundary
#define Y_WALL (0.00005)
#define SQCMU (0.3)
#define KAPPA_VK (0.4)

// kappa - epsilon model ================================
#if ((TBK_EQUATIONS/2)==1)
// kappa omega constant
#define SIGMAO (2.)
#define SIGMAK (2.)
#define CT2O (1.)
#define CT1O (1.)
#define CMU (0.09)
#define BETA (0.09)
#define BETASTAR (0.09)
#endif

// kappa -omega model ================================
#if ((TBK_EQUATIONS/2)==2)

//    #define SST (1)
   #define A1 (1000.)

//   #define LOWRE (1)

// #define ALFA_W (0.555555555556)
#define ALFA_W (0.52)
// #define BETA_W (0.075)
#define BETA_W (0.072)
#define SIGMA_W (.5)

#define BETASTAR (0.09)
#define SIGMA_K (.5)

#endif
#endif
#endif  // ============= end MGSTTBK_EQUATIONS  ============

