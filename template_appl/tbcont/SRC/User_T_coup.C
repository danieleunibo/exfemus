#include "Equations_conf.h"

// ============================================
#ifdef T_COUP_EQUATIONS // 3D-2D Energy adjoint equation
#if T_COUP_EQUATIONS == 1  // standard quadratic elements
// ============================================
#include <sstream>
#include "MGGeomEl.h"
// configuration files -----------
#include "MGSclass_conf.h"
#include "MGFE_conf.h"
// #include "MGSTBKconf.h"
// #include "MGSTTBKconf.h"   // turbulent energy
#include "Printinfo_conf.h"
#include "MGEquationsSystem.h"

// class local configuration -------
#include "MGSolverT_COUP.h"

// local include -------------------
#include "MGMesh.h"
#include "MGSystem.h"
#include "numeric_vectorM.h"
#include "dense_vectorM.h"
#include "sparse_matrixM.h"
#include "dense_matrixM.h"
#include "linear_solverM.h"

#include "parallelM.h"




// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
// =================================================
void MGSolT_COUP::ic_read(int bc_gam,int bc_mat, double xp[],double u_value[]) {
// =================================================
#if DIMENSION==2
// xp[]=(xp,yp) u_value[]=(T,T_ad)
 
  u_value[0] =  1.;
  u_value[1] =  0.;

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
  
//   //    boundary conditions box
//   if (xp[0]> LXE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;       bc_flag[0]=0;  bc_Neum[1]=1;       bc_flag[1]=0; }   // right 
//   if (xp[0]< LXB*ILref+BDRY_TOLL) {    bc_Neum[0]=0;       bc_flag[0]=0;  bc_Neum[1]=0;       bc_flag[1]=0; }   // left
//   if (xp[1]> LYE*ILref-BDRY_TOLL) {    bc_Neum[0]=0;       bc_flag[0]=1;  bc_Neum[1]=0;       bc_flag[1]=0; }   // top
//   if (xp[1]< LYB*ILref+BDRY_TOLL) {    bc_Neum[0]=0;       bc_flag[0]=1;  bc_Neum[1]=0;       bc_flag[1]=0; }   // bottom
  
  // boundary layer box
  if (xp[0]> LXE*ILref-BDRY_TOLL) {bc_Neum[0]=1;  bc_flag[0]=0;    bc_Neum[1]=1;  bc_flag[1]=0; }   // right 
  if (xp[1]> LYE*ILref-BDRY_TOLL) {bc_Neum[0]=1;  bc_flag[0]=0;    bc_Neum[1]=1;  bc_flag[1]=0; }   // top
  if (xp[0]< LXB*ILref+BDRY_TOLL) {bc_Neum[0]=0;  bc_flag[0]=2;    bc_Neum[1]=0;  bc_flag[1]=0; }   // left
  if (xp[1]< LYB*ILref+BDRY_TOLL) {bc_Neum[0]=0;  bc_flag[0]=1;    bc_Neum[1]=0;  bc_flag[1]=0; }   // bottom
  
  
  // channel with injection
//   if (xp[0]> LXE*ILref-BDRY_TOLL) {bc_Neum[0]=1;  bc_flag[0]=0;    bc_Neum[1]=1;  bc_flag[1]=0; }   // right 
//   if (xp[1]> LYE*ILref-BDRY_TOLL) {bc_Neum[0]=0;  bc_flag[0]=4;    bc_Neum[1]=0;  bc_flag[1]=0; }   // top
//   if (xp[0]< LXB*ILref+BDRY_TOLL) {bc_Neum[0]=0;  bc_flag[0]=2;    bc_Neum[1]=0;  bc_flag[1]=0; }   // left
//   if (xp[1]< LYB*ILref+BDRY_TOLL) {bc_Neum[0]=0;  bc_flag[0]=4;    bc_Neum[1]=0;  bc_flag[1]=0; }   // bottom
//   if (xp[1]< LYB*ILref+BDRY_TOLL && xp[0] > 1.-BDRY_TOLL && xp[0] <2.+ BDRY_TOLL) {bc_Neum[0]=0;  bc_flag[0]=1;    bc_Neum[1]=0;  bc_flag[1]=0; }   // bottom
//   
//   

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
