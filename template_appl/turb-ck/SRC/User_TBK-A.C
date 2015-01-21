#include "Equations_conf.h"  // <--- Equations configure

#ifdef TBKA_EQUATIONS

// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverTBK-A.h"       // Navier-Stokes class header file


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



// =========================================================
/// This function generates the initial conditions for the NS system:
void MGSolTBKA::ic_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  double u_value[]
) {// =======================================
//xp[] is the NON-DIMENSIONAL node coordinate
  double ILref = 1./_lref;

#if DIMENSION==2 // --------------  2D //--------------------------
  // xp[]=(xp,yp) u_value[]=(u,v,p)

#if (TBKA_EQUATIONS%2==0)                // splitting 
  if(_dir==0)  { //  kappa
    u_value[0] = .1;
  }
  if(_dir==1)   { // omega
    u_value[0] = .1;
//     if(xp[1]<0.01) {
//        if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) { u_value[0] =4.; }
//      }
  }
#else                          // --------- coupled ---------------------------
  //  kappa -> 0  epsilon -> 1
  u_value[0] = 1.e-2;
//   u_value[1] = 1.e-3;
//   if(xp[1]<0.01) {
//     u_value[1] =4.;
//     if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) {  u_value[0] =.1; }
//   }
//    if (xp[1]<0.001 && xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL)
//     u_value[1] =6852.79932546*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref;
//      if (xp[1]<0.001 && xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL)
//   if(xp[1]<0.01)
//     if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) u_value[1] =0.777+ 0*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref;
#endif


#endif  // //-----------------------//------------------------------------

// ===================================================================
#if DIMENSION==3  // --------------  3D //--------------------------

#if (TBKA_EQUATIONS%2==0)    // ------- splitting ------------------- 
  if(_dir==0)  {//  kappa +++++++++++++++++++++++++
    u_value[0] = 0.;
  }
  if(_dir==1)   { // omega ++++++++++++++++++++++++
    u_value[0] = 0.;
    if(xp[1]<0.01) {
      if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) { u_value[0] =0.777+ 0*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref; }
    }
  }
#else               // --------- coupled ---------------------------
  // xp[]=(xp,yp,zp)   kappa -> 0  omega -> 1
  u_value[0]= 0.;
  u_value[1] =0.;
  if(xp[1]<0.01) {
    if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) { u_value[0] =0.777+ 0*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref; }
  }
#endif

#endif // //-------------------------------------------
  return;
}

// ========================================
/// This function  defines the boundary conditions for the NS system:
void MGSolTBKA::bc_read(
  int bc_gam,
  int bc_mat,
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], // normal
  int bc_flag[]         // boundary condition flag
) {// ===================================
  //     0 ->  single component
  //     +4 ->  nonhomogeneous
  //     +2 -> tg
  //     +1 ->  normal
// ===================================
  double ILref = 1./_lref;

#if DIMENSION==2  // ----- 2D boundary conditions ----------
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann

#if (TBKA_EQUATIONS%2==0)        // splitting 
  if(_dir==0) { //  kappa
    
    if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  // top 
    if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  //  left
    if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  //  right
    if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  // bottom
    
          // cyl med
//     if (bc_gam==300)  {bc_flag[0]=0; bc_Neum[0]=1;} //  top         
//     if (bc_gam==200)  {bc_flag[0]=0; bc_Neum[0]=1;} //  bottom      
//     if (bc_gam==400)  {bc_flag[0]=0; bc_Neum[0]=1;} //  left        
//     if (bc_gam==100)  {bc_flag[0]=0; bc_Neum[0]=0;} //  right 
//     
  }
  if(_dir==1) {  // omega
    if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  // top
    if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=1; bc_Neum[0]=1; }  //  left
    if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=1; bc_Neum[0]=1; }  //  right
    if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=2; bc_Neum[0]=0; }  // bottom
  }

#else           // --------- coupled ---------------------------
  //  kappa -> 0  omega -> 1
//   // cyl
//   if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=1; bc_Neum[1]=1;} //  top
//   if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=1; bc_Neum[1]=1;} //  bottom
//   if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=1; bc_Neum[1]=1;} //  left
//   if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=1; bc_flag[1]=0; bc_Neum[0]=1; bc_Neum[1]=0;} //  right
// 
//   
//   // cyl med
//  if (bc_gam==300)  {bc_flag[0]=0; bc_Neum[0]=1;} //  top         
//  if (bc_gam==200)  {bc_flag[0]=0; bc_Neum[0]=1;} //  bottom      
//  if (bc_gam==400)  {bc_flag[0]=0; bc_Neum[0]=1;} //  left        
//  if (bc_gam==100)  {bc_flag[0]=1; bc_Neum[0]=1;} //  right 
  
  
#endif // //---------------------------------------------

#endif // //---------------------------------------------

#if DIMENSION==3 // -------- 3D boundary conditions ---------
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
  //    boundary conditions box
#if (TBKA_EQUATIONS%2==0)          // splitting 
  if(_dir==0) {  //  kappa

    if(xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  //   OUTLET
    if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  // 
    if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  // 
    if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  // 
    if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  //  
    if(xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  //   INLET 
    
//       bc_flag[0]=0; bc_Neum[0]=0; // wall
//   if (bc_gam==400) bc_flag[0]=1; bc_Neum[0]=0;
//   if (bc_gam==300) bc_flag[0]=1; bc_Neum[0]=0;
  
//   if (xp[2]< LZB*ILref + BDRY_TOLL){bc_flag[0]=0; bc_Neum[0]=0;} //  INLET
//   if (xp[2]> LZE*ILref - BDRY_TOLL){bc_flag[0]=0; bc_Neum[0]=1;} //  OUTLET
    
  }
  if(_dir==1) { // omega
    if(xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  //  left
    if(xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  //  right
    if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; }  // top
    if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=2; bc_Neum[0]=0; }  // bottom
    if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  //  left
    if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; }  //  right
  }
#else            // --------- coupled ---------------------------
  if(xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=2; bc_Neum[0]=0; bc_Neum[1]=0;} //  INLET
  if(xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=0; bc_Neum[1]=1;} //  OUTLET
  if(xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=0; bc_Neum[1]=0;}
  if(xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=0; bc_Neum[1]=0;}
  if(xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=1; bc_Neum[1]=1;} //current
  if(xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=1; bc_Neum[1]=1;} //current
#endif
#endif  // //-----------------------//---------------------------

    return;
  }

#endif // ----------  endif TBKA_EQUATIONS  