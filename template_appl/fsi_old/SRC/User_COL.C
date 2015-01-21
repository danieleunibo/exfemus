#include "Equations_conf.h"



// ============================================
#ifdef COLOR_EQUATIONS // 3D-2D Energy equation
// ============================================

#include "MGSolverCOL.h"
#include "MGGeomEl.h"

// configuration files -----------
#include "MGFE_conf.h"
#include "Printinfo_conf.h"
#include "MGEquationsSystem.h"

// class local configuration -------


// local include -------------------
// #include "MGMesh.h"
 #include "MGSystem.h"
// #include "numeric_vectorM.h"
// #include "dense_vectorM.h"
// #include "sparse_matrixM.h"
// #include "dense_matrixM.h"
// #include "linear_solverM.h"
// #include "parallelM.h"





// ========================================
/// This function  defines the boundary conditions for the NS system:
void MGSolCOL::bc_intern_read(
  int bc_gam,
  int mat_gam,
  double xp[],   // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], // normal
  int bc_flag[]  // boundary condition flag
) {// ===================================
//     mat_gam=2;
//   if(xp[0]>0.2-0.0001  && xp[0]<0.29+0.0001 && xp[1]>0.3-0.0001  && xp[1]< 0.7+0.0001&& xp[2]< 0.4+0.0001 ) {mat_gam=4;}
// for(int ivar=0;ivar<_n_vars;ivar++) {
// //     if (xp[0]<1-1.e-6 ||xp[0]> 1.01 +1.e-6 || (xp[1]<0.51 -1.e-6 && xp[1]>0.49 +1.e-6) ) { // fluid
//   if(mat_gam==2){
//     bc_flag[ivar]=1;     bc_Neum[ivar]=1 ;
//   } else {
// //     if (xp[0]<1+1.e-6 ||xp[0]> 1.01 -1.e-6 || (xp[1]<0.51 +1.e-6 && xp[1]>0.49 -1.e-6) )  {  // interface boundary
//       if(bc_gam==1000){
//       bc_flag[ivar]=1;        bc_Neum[ivar]=5;
//     } else { //solid
//       bc_flag[ivar]=1;         bc_Neum[ivar]=3 ;
//     }
//   }
// }
  
     if (mat_gam==2){ bc_flag[0]=1;bc_Neum[0]=1; }
   if (mat_gam==4){ bc_flag[0]=1;bc_Neum[0]=3 ;}
   if (bc_gam==1000)   { bc_flag[0]=1;bc_Neum[0]=5;}
  return;
}




void MGSolCOL::ic_read(
  int /*bc_gam*/,
  int /*mat_gam*/,
  double /*xp*/[],
  double u_value[]
) {// =======================================
  for(int ivar=0;ivar<_n_vars;ivar++)  u_value[ivar] = 0.;
  return;
}

// ========================================
/// This function  defines the boundary conditions for the NS system:
void MGSolCOL::bc_read(
  int bc_gam,
  int mat_gam,
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], // normal
  int bc_flag[]         // boundary condition flag
) {// ===================================
  //     0 ->  single component
  //     +4 ->  nonhomogeneous
  //     +2 -> tg
  //     +1 ->  normal
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
//     mat_gam=2;
//   if(xp[0]>0.2-0.0001  && xp[0]<0.29+0.0001 && xp[1]>0.3-0.0001  && xp[1]< 0.7+0.0001&& xp[2]< 0.4+0.0001 ) {mat_gam=4;}
// for(int ivar=0;ivar<_n_vars;ivar++) {
// //     if (xp[0]<1-1.e-6 ||xp[0]> 1.01 +1.e-6 || (xp[1]<0.51 -1.e-6 && xp[1]>0.49 +1.e-6) ) { // fluid
//   if(mat_gam==2){
//     bc_flag[ivar]=1;     bc_Neum[ivar]=1 ;
//   } else {
// //     if (xp[0]<1+1.e-6 ||xp[0]> 1.01 -1.e-6 || (xp[1]<0.51 +1.e-6 && xp[1]>0.49 -1.e-6) )  {  // interface boundary
//       if(bc_gam==1000){
//       bc_flag[ivar]=1;        bc_Neum[ivar]=5;
//     } else { //solid
//       bc_flag[ivar]=1;         bc_Neum[ivar]=3 ;
//     }
//   }
// }

   if (mat_gam==2){ bc_flag[0]=1;bc_Neum[0]=1; }
   if (mat_gam==4){ bc_flag[0]=1;bc_Neum[0]=3 ;}
   if (bc_gam==1000)   { bc_flag[0]=1;bc_Neum[0]=5;}
  return;
 }
#endif
