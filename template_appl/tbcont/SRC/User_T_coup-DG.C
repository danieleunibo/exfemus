#include "Equations_conf.h"

// ============================================
#ifdef T_COUP_EQUATIONS // 3D-2D Energy adjoint equation
#if T_COUP_EQUATIONS == 2  // Discontinuous Galerkin
// ============================================

// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverT_COUP.h"       // Navier-Stokes class header file


// config file --------------------------------------------------------------------------
#include "MGGeomEl.h"        // Geometrical element
#include "MGFE_conf.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include ------------------------------------------------------------
#include "MGMesh.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsSystem.h"  // Equation map class
#include "MGUtils.h"
// standard lib -------------------------------------------------------------------------
#include <string.h>          // string library

// local alg lib ------------------------------------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ======================================================================================



// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
// =================================================
void MGSolT_COUP::ic_read(int bc_gam,int bc_mat, double xp[],double u_value[]) {
// =================================================
#if DIMENSION==2
// xp[]=(xp,yp) u_value[]=(T,T_ad)
 
  u_value[0] =  0.;
  u_value[1] =  0.7;
  
//   if (((xp[1]-0.2)*(xp[1]-0.2)+(xp[0]-0.3)*(xp[0]-0.3)) < 0.05) u_value[1] = 1.;    //sphere
//   if (xp[0]<0.5) u_value[1] = 1.;    //wave

#else
// =================================================
  // xp[]=(xp,yp,zp) u_value[]=(T,T_ad)
  
  u_value[0] =  1.;
  u_value[1] =  0.;

#endif
}

// ========================================================
/// This function  defines the boundary conditions for the energy equation:
void MGSolT_COUP::bc_read(int bc_gam,int bc_mat, double xp[],int bc_Neum[],int bc_flag[]) {
  // =================================================
  double ILref = 1./_lref;
#if DIMENSION==2
// xp[]=(xp,yp) bc_flag[]=0-Dirichlet 1-Neumann
  
// //   //    boundary conditions box
  if (xp[1]> LYE*ILref-BDRY_TOLL) { bc_Neum[0]=1;       bc_flag[0]=0;  bc_Neum[1]=1;       bc_flag[1]=0; }   // top
  if (xp[0]< LXB*ILref+BDRY_TOLL) { bc_Neum[0]=0;       bc_flag[0]=0;  bc_Neum[1]=1;       bc_flag[1]=0; }   // left
  if (xp[0]> LXE*ILref-BDRY_TOLL) { bc_Neum[0]=1;       bc_flag[0]=0;  bc_Neum[1]=1;       bc_flag[1]=0; }   // right 
  if (xp[1]< LYB*ILref+BDRY_TOLL) { bc_Neum[0]=0;       bc_flag[0]=0;  bc_Neum[1]=0;       bc_flag[1]=0; }   // bottom

//   if (xp[1] < LYB*ILref + BDRY_TOLL && xp[0] > 1.-BDRY_TOLL && xp[0] < 2.+ BDRY_TOLL)
//       { bc_Neum[0]=0; bc_flag[0]=0; bc_Neum[1]=0; bc_flag[1]=0; } // injection 
// //   channel with injection
//   if (bc_gam==10) {bc_Neum[0]=0;       bc_flag[0]=0;  bc_Neum[1]=0;       bc_flag[1]=1;} //little injection
//   if (bc_gam==20) {bc_Neum[0]=0;       bc_flag[0]=0;  bc_Neum[1]=1;       bc_flag[1]=0;} //big injection (left)
//   if (bc_gam==30) {bc_Neum[0]=0;       bc_flag[0]=0;  bc_Neum[1]=0;       bc_flag[1]=1;} //wall
//   if (bc_gam==40) {bc_Neum[0]=1;       bc_flag[0]=0;  bc_Neum[1]=1;       bc_flag[1]=0;} //outlet
//   if (xp[0] > 5. - BDRY_TOLL && (xp[1] <0. + BDRY_TOLL || xp[1] >1. - BDRY_TOLL))  //corners
//    {bc_Neum[0]=0;       bc_flag[0]=0;  bc_Neum[1]=0;       bc_flag[1]=1;} //wall
  

#endif

#if DIMENSION==3
// =================================================
  // xp[]=(xp,yp) bc_flag[]=0-Dirichlet 1-Neumann
  //    boundary conditions box

  if (xp[1]> LYE*ILref-BDRY_TOLL)  {
    bc_flag[0]=0;    // top
    bc_Neum[0]=1;
  }
  if (xp[1]< LYB*ILref+BDRY_TOLL)  {
    bc_flag[0]=2;    // bottom
    bc_Neum[0]=1;
  }
  if (xp[2]< LZB*ILref+BDRY_TOLL)  {
    bc_flag[0]=0;
    bc_Neum[0]=1;
  }
  if (xp[2]> LZE*ILref-BDRY_TOLL)  {
    bc_flag[0]=0;
    bc_Neum[0]=1;
  }
  if (xp[0]< LXB*ILref+BDRY_TOLL)  {
    bc_flag[0]=0;    // left
    bc_Neum[0]=1;
  }
  if (xp[0]> LXE*ILref-BDRY_TOLL)  {
    bc_flag[0]=0;    // right
    bc_Neum[0]=1;
  }

  
#endif

  return;
} // end boundary conditions ==========================


#endif
#endif
