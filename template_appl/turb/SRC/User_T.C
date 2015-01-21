

#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS // 3D-2D Energy equation
// ============================================

// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverT.h"       // Navier-Stokes class header file


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
void MGSolT::ic_read(int bc_gam,int bc_mat, double xp[],double u_value[]) {
// =================================================
#if DIMENSION==2
// xp[]=(xp,yp) u_value[]=(u,v,p)
//   u_value[0] =573.15;
  
  u_value[0] =  1.;//*(1.-xp[1]);
//    if (xp[1]<0.001) u_value[0] = 573.;
     
#else
// =================================================
  // xp[]=(xp,ypzp) u_value[]=(u,v,w,p)
//   u_value[0] = 0.*xp[0]* (1-xp[0]) *xp[1]* (2-xp[1]);
//   if (xp[2]<0.1)
  u_value[0] = 573.;
#endif
}

// ========================================================
/// This function  defines the boundary conditions for the energy equation:
void MGSolT::bc_read(int bc_gam,int bc_mat, double xp[],int bc_Neum[],int bc_flag[]) {
  // =================================================
  double ILref = 1./_lref;
#if DIMENSION==2
// xp[]=(xp,yp) bc_flag[T]=0-Dirichlet 1-Neumann
  //    boundary conditions box
// if (xp[0]> LXE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  }   // right  
//   if (xp[1]> LYE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  }  // top
//   if (xp[0]< LXB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=2;  }// left
//   if (xp[1]< LYB*ILref+BDRY_TOLL) {    bc_Neum[0]=0;       bc_flag[0]=2;  }// bottom

  // cyl
 
  if (xp[1]> LYE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  }  // top
  if (xp[0]< LXB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  }// left
  if (xp[1]< LYB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  }// bottom
  if (xp[0]> LXE*ILref-BDRY_TOLL) {    bc_Neum[0]=0;       bc_flag[0]=0;  }   // right  

// ext q
 //if (xp[0]> LXE*ILref-BDRY_TOLL)  {     bc_Neum[0]=0;       bc_flag[0]=0;  }   // right  
 // if (xp[1]> LYE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  }  // top
 // if (xp[0]< LXB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  }// left
 // if (xp[1]< LYB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  }// bottom

  // int q
//   if (xp[0]> LXE*ILref-BDRY_TOLL)  {     bc_Neum[0]=0;       bc_flag[0]=0;  }   // right  
//   if (xp[1]> LYE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  }  // top
//   if (xp[0]< LXB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  }// left
//   if (xp[1]< LYB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  }// bottom

  
//   if (xp[0]< BDRY_TOLL && xp[1]< 2.14+BDRY_TOLL) {// left bottom
//     bc_Neum[0]=1;    
//     bc_flag[0]=0;
//   }
  
//   if (xp[0]< LXB*ILref+BDRY_TOLL && xp[1]> 1.059+BDRY_TOLL) {bc_flag[0]=2; bc_Neum[0]=1;} // left top
//   if (xp[0]< BDRY_TOLL && xp[1]< 1.+BDRY_TOLL) { bc_flag[0]=0; bc_Neum[0]=1;} //left bottom
//   if (xp[1] > 1. - BDRY_TOLL && xp[1] < 1.059 + BDRY_TOLL && xp[0]< LXB*ILref+BDRY_TOLL) {bc_flag[0]=2;bc_Neum[0]=0;} //gradino
  
  
 
#endif

#if DIMENSION==3
// =================================================
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
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
  
  // square bundle
     double p=0.00615; //P/D 1.5
//      double p=0.00533;  //P/D 1.3
//   double p=0.005;   //P/D 1.2
  bc_flag[0]=0;bc_Neum[0]=1;
   if (xp[1]*xp[1]+(xp[0]-p)*(xp[0]-p)<  (0.0041+0.5e-6)*(0.0041+0.5e-6)){
      bc_flag[0]=2;bc_Neum[0]=1;
         }
   if(  xp[2] < 1.e-6 ) { bc_flag[0]=0;bc_Neum[0]=0;}
//   

//   //  tri bundle  
//   bc_flag[0]=0;bc_Neum[0]=1;
// //      if (xp[1]*xp[1]+(xp[0]-0.00533)*(xp[0]-0.00533)< 0.0041*0.0041+BDRY_TOLL){                           //rod arc   P/D=1.3
// //      if (xp[1]*xp[1]+(xp[0]-0.00574)*(xp[0]-0.00574)< 0.0041*0.0041+BDRY_TOLL){                           //rod arc   P/D=1.4
//         if (xp[1]*xp[1]+(xp[0]-0.00615)*(xp[0]-0.00615)< 0.0041*0.0041+BDRY_TOLL){                           //rod arc    P/D=1.5
//             bc_flag[0]=2;bc_Neum[0]=1;}
//    if(  xp[2] < 1.e-6 ) { bc_flag[0]=0;bc_Neum[0]=0;} //inlet    
// //     if(  xp[1] < 1.e-8 && xp[0] < 0.00615-0.0041-1.e-6 ) { bc_flag[0]=0;bc_Neum[0]=1;} 
// //     if(   xp[1] >0.003551+ 1.e-8) { bc_flag[0]=0;bc_Neum[0]=1;}
  
#endif

  return;
} // end boundary conditions ==========================


#endif

