// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

// class files --------------------------------------------------
#include "MGSolverNS.h"       // Navier-Stokes class header file


// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation   
#include "Printinfo_conf.h"  // Print options

// local Femus class include -----------------------------------
#include "MGMesh.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsMap.h"  // Equation map class

// standard lib -----------------------------------------------
#include <string.h>          // string library

// local alg lib -----------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ==============================================================


// =========================================================
/// This function generates the initial conditions for the NS system:
void MGSolNS::ic_read(
  int bc_gam,
  int mat_gam, 
  double xp[],
  double u_value[]
) {// =======================================
//xp[] is the NON-DIMENSIONAL node coordinate
  double ILref = 1./_lref;
  u_value[0]=0.; u_value[1]=0.; u_value[2]=0.; u_value[3]=0.; // DEFAULT IC

  std::string msh_name=_mgutils.get_file("MESHNAME");

  if(msh_name=="elsy") {

//     double modulo = 0.85;// modulo is the value of inlet velocity

//   if(xp[2]>3.5){  u_value[0]=3.; u_value[1]=2.; u_value[2]=1.; u_value[3]=0;}

    if(bc_gam==4000) {
      u_value[2]=-0.7; 

    }

//       if(bc_gam==3100){
// 	u_value[0]=0.; u_value[1]=-modulo; u_value[2]=0.; u_value[3]=0.;}
//
//       if(bc_gam==3200){
// 	u_value[0]=-modulo/sqrt(5.); u_value[1]=-2*(modulo/sqrt(5.)); u_value[2]=0.; u_value[3]=0.;}
//
//       if(bc_gam==3300){
// 	u_value[0]=-modulo/sqrt(2.); u_value[1]=-modulo/sqrt(2.); u_value[2]=0.; u_value[3]=0.;}
//
//       if(bc_gam==3400){
// 	u_value[0]=-2*(modulo/sqrt(5.)); u_value[1]=-modulo/sqrt(5.); u_value[2]=0.; u_value[3]=0.;}
//
//       if(bc_gam==3500){
// 	u_value[0]=-modulo; u_value[1]=0.; u_value[2]=0.; u_value[3]=0.;}

  }

  return;
}

// ========================================
// This function  defines the boundary conditions for the NS system:
void MGSolNS::bc_read(
  int bc_gam,
  int mat_gam, 
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[],        // normal
  int bc_flag[]         // boundary condition flag
) {

  const double Lref = _mgphys.get_par("Lref");
  double ILref = 1./Lref;


  bc_flag[0]=0; bc_flag[1]=0; bc_flag[2]=0;  bc_flag[3]=0;  // DEFAULT BC
  bc_Neum[0]=1; bc_Neum[1]=1; bc_Neum[2]=1;  bc_Neum[3]=1;  //     "


  std::string msh_name=_mgutils.get_file("MESHNAME");

  if(msh_name=="box") {
    if(xp[0]<LXB*ILref+BDRY_TOLL) {
      bc_flag[0]=0; bc_flag[2]=0;
      bc_Neum[0]=0; bc_Neum[2]=1;
    }
    if(xp[1]<LYB*ILref+BDRY_TOLL) {
      bc_flag[1]=0; bc_flag[2]=0;
      bc_Neum[1]=0; bc_Neum[2]=1;
    }
    if(xp[0]>LXE*ILref-BDRY_TOLL) {
      bc_flag[0]=0; bc_flag[1]=0; bc_flag[2]=0;
      bc_Neum[0]=0; bc_Neum[1]=0; bc_Neum[2]=1;
    }
    if(xp[1]>LYE*ILref-BDRY_TOLL) {
      bc_flag[0]=0; bc_flag[1]=0; bc_flag[2]=0;
      bc_Neum[0]=0; bc_Neum[1]=0; bc_Neum[2]=1;
    }
    if(xp[2]<LZB*ILref+BDRY_TOLL) {
      bc_flag[0]=0; bc_flag[1]=0; bc_flag[2]=4;
      bc_Neum[0]=0; bc_Neum[1]=0; bc_Neum[2]=0;
    }
  }

  if(msh_name=="elsy_core") { // ================================ elsy core ===============

    if(bc_gam==50) {
      bc_flag[0]=0; bc_flag[1]=0;
      bc_Neum[0]=0; bc_Neum[1]=0;
    }

    if(bc_gam==60) {
      bc_flag[0]=0; bc_flag[1]=0;
      bc_Neum[0]=0; bc_Neum[1]=0;
    }
   // ------------------------------------------------------------------
    if(bc_gam==70) { 
      bc_flag[0]=0; bc_flag[1]=0; bc_flag[2]=0;
      bc_Neum[0]=0; bc_Neum[1]=0; bc_Neum[2]=0;
    }
    // SIMMETRY X -------------------------------------------------------------------------
    if(xp[0]<0.*ILref+BDRY_TOLL) {  
      bc_flag[0]=0;
      bc_Neum[0]=0;
    }
   // SIMMETRY Y --------------------------------------------------------------------------
   if(xp[1]<0.*ILref+BDRY_TOLL) {  
      bc_flag[1]=0;
      bc_Neum[1]=0;
    }
   // INLET
   if(xp[2]<0.*ILref+BDRY_TOLL) {                     
      bc_flag[0]=0; bc_flag[1]=0; bc_flag[2]=4;
      bc_Neum[0]=0; bc_Neum[1]=0; bc_Neum[2]=1;
    }

  }
  // ================================ elsy ================================================
  if(msh_name=="elsy") {
    // X SIMMETRY PLANE -------------------------------------------------------------------
    if(xp[0] < BDRY_TOLL) {           
      bc_flag[0]=0; /*bc_flag[1]=0;*/
      bc_Neum[0]=0; /*bc_Neum[1]=0;*/
    }
    // Y SIMMETRY PLANE -------------------------------------------------------------------
    if(xp[1] < BDRY_TOLL) {           
      /*bc_flag[0]=0;*/ bc_flag[1]=0;
      /*bc_Neum[0]=0;*/ bc_Neum[1]=0;
    }
   // CONTROL RODS & CORE LATERAL SURFACES ------------------------------------------------
    if((bc_gam==200)||(bc_gam==300)) {   
      bc_flag[0]=0; bc_flag[1]=0;
      bc_Neum[0]=0; bc_Neum[1]=0;
    }
    // CONTROL RODS BOTTOM SURFACE --------------------------------------------------------
    if(bc_gam==400) {                              
      bc_flag[0]=0; bc_flag[1]=0; bc_flag[2]=0;
      bc_Neum[0]=0; bc_Neum[1]=0; bc_Neum[2]=0;
    }
    // TOP SURFACE ------------------------------------------------------------------------
    if (bc_gam==2100){
/*       bc_flag[0]=0; bc_flag[1]=0;*/  bc_flag[2]=0;
/*       bc_Neum[0]=0; bc_Neum[1]=0;*/  bc_Neum[2]=0;
    }
    // INLET (bottom old plenum) ----------------------------------------------------------
    // if ((bc_gam==3100)||(bc_gam==3200)||(bc_gam==3300)||(bc_gam==3400)||(bc_gam==3500)){
    //     bc_flag[0]=4; bc_flag[1]=4; bc_flag[2]=0;
    //     bc_Neum[0]=1; bc_Neum[1]=1; bc_Neum[2]=0;
    //  }
    // INLET (bottom new plenum) ----------------------------------------------------------
    if(bc_gam==4000) {
      bc_flag[0]=0;  bc_flag[1]=0;  bc_flag[2]=4;
      bc_Neum[0]=0;  bc_Neum[1]=0;  bc_Neum[2]=0;
    }

    // WALLS (top plenum) -----------------------------------------------------------------
    if((bc_gam==1300)||(bc_gam==1400)||(bc_gam==1500)||(bc_gam==4100)||(bc_gam==4200)) {
      bc_flag[0]=0; bc_flag[1]=0; bc_flag[2]=0;
      bc_Neum[0]=0; bc_Neum[1]=0; bc_Neum[2]=0;
    }
    // OUTLET (top plenum= 1100 (right)  and 1200 (left) ----------------------------------
    if((bc_gam==1100)||(bc_gam==1200)) {
      bc_flag[0]=0; bc_flag[1]=0; bc_flag[2]=0;
      bc_Neum[0]=1; bc_Neum[1]=1; bc_Neum[2]=0;
      bc_Neum[3]=10;    bc_flag[3]=0;
    }
    // TRANSVERSAL SURFACES (z=cost) ------------------------------------------------------
    if((bc_gam==2400)||(bc_gam==2200)||(bc_gam==2300)/*||(bc_gam==2100)*/) {
      bc_flag[0]=0; bc_flag[1]=0; bc_flag[2]=0;
      bc_Neum[0]=0; bc_Neum[1]=0; bc_Neum[2]=0;
    }
    // OUTLET VELOCITY NORMAL -------------------------------------------------------------
    // if ((bc_gam==1100)){   // OUTLET VELOCITY NORMAL
    //   bc_flag[0]=2; bc_flag[1]=0; bc_flag[2]=0;
    //   bc_Neum[0]=0; bc_Neum[1]=1; bc_Neum[2]=0;
    //  }
    // OUTLET VELOCITY NORMAL -------------------------------------------------------------
    // if ((bc_gam==1200)){   // OUTLET VELOCITY NORMAL
    //   bc_flag[0]=1; bc_flag[1]=0; bc_flag[2]=0;
    //   bc_Neum[0]=0; bc_Neum[1]=1; bc_Neum[2]=0;
    //   }


  }

  return;
}



// ===============================================================
// --------------   PRESSURE EQUATION [P_F] ------------------
// ===============================================================
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================
#if (NS_EQUATIONS%2==0)


// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
// =================================================
void MGSolNSP::ic_read(double xp[],double u_value[]) {
// xp[]=(xp,yp) u_value[]=(u,v,p)
  u_value[0] =0.;
}

// ========================================================
/// This function  defines the boundary conditions for the energy equation:
void MGSolNSP::bc_read(double xp[],int bc_Neum[],int bc_flag[]) {
  // =================================================
  const double Lref = _mgphys.get_par("Lref");
  double ILref = 1./Lref;
#if DIMENSION==2
// xp[]=(xp,yp) bc_flag[T]=0-Dirichlet 1-Neumann
  //    boundary conditions box

  if(xp[0]< LXB*ILref+BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }   // left
  if(xp[0]> LXE*ILref-BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }   // right
  if(xp[1]< LYB*ILref+BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }   // bottom
  if(xp[1]> LYE*ILref-BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=0;  }   // top
#endif

#if DIMENSION==3
// =================================================
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
  //    boundary conditions box
  if(xp[0] < LXB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
  if(xp[0] > LXE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
  if(xp[1] < LYB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
  if(xp[2] < LZB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
  if(xp[2] > LZE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
  if(xp[1] > LYE*ILref-BDRY_TOLL) {    bc_Neum[0]=0;    bc_flag[0]=0;  }
#endif
  return;
} // end boundary conditions ==========================





#endif // ENDIF NS_EQUATIONS==0
#endif  //ENDIF NS_EQUATIONS

