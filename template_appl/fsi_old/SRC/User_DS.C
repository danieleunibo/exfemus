#include "Equations_conf.h"



// ============================================
#ifdef DS_EQUATIONS // 3D-2D Energy equation
// ============================================

#include "MGSolverDS.h"
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
void MGSolDS::bc_intern_read(
  int bc_gam,
  int mat_gam,
  double xp[],   // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], // normal
  int bc_flag[]  // boundary condition flag
) {// ===================================
// if (xp[0]<1.-1e-6) mat_gam=2;
// if (xp[0]>1.-1e-6) mat_gam=4;
//default as before
// for(int ivar=0;ivar<_n_vars;ivar++) {
// //    bc_flag[ivar]=1;bc_Neum[ivar]=1;
// //    if (mat_gam==4){
// //      bc_flag[ivar]=1;bc_Neum[ivar]=3;
// //   }
//     //   if (bc_gam==100)   interface
//    if (xp[0]<1-1.e-6){ bc_flag[ivar]=1;bc_Neum[ivar]=1; }
//    else{
//    if (xp[0]>1+1.e-6){ bc_flag[ivar]=1;bc_Neum[ivar]=3 ;}
//    else   { bc_flag[ivar]=1;bc_Neum[ivar]=5;}
//    }
#if DIMENSION == 2
// // // // // pressure box
//    if (xp[0]<1-1.e-6){ bc_flag[0]=1;bc_Neum[0]=1; }
//    if (xp[0]>1+1.e-6){ bc_flag[0]=1;bc_Neum[0]=3 ;}
//    if (fabs(xp[0]-1.)<1.e-6)   { bc_flag[0]=1;bc_Neum[0]=5;}

// // // // bar.neu box
   if (mat_gam==2){ bc_flag[0]=1;bc_Neum[0]=1; }
   if (mat_gam==4){ bc_flag[0]=1;bc_Neum[0]=3 ;}
   if (bc_gam==1000)   { bc_flag[0]=1;bc_Neum[0]=5;}
// // // valv_vc_mat
//   if (xp[0]<1-1.e-6 ||xp[0]> 1.01 +1.e-6 || (xp[1]<0.51 -1.e-6 && xp[1]>0.49 +1.e-6) ) { // fluid
  
  
  
//    bc_flag[0]=0;     bc_Neum[0]=1 ;
  
  
// //turek
// if(mat_gam==2){
//     bc_flag[0]=0;     bc_Neum[0]=1 ;
// //     bc_flag[1]=1;     bc_Neum[1]=1 ;
// //     bc_Neum[2]=1;    // pressure flag
//   } else {
// //solid
//       bc_flag[0]=1;         bc_Neum[0]=3 ;
// //       bc_flag[1]=1;      bc_Neum[1]=3 ;
// //       bc_Neum[2]=4;  // pressure flag
//     }
//       if(bc_gam==1000){
//       bc_flag[0]=1;        bc_Neum[0]=5;
// //       bc_flag[1]=1;      bc_Neum[1]=5;
// //       bc_Neum[2]=5;  // pressure flag
//     }
// //turek
#endif
#if DIMENSION==3
// // // // bar.neu box
  mat_gam=2;
  if(xp[0]>0.2-0.0001  && xp[0]<0.29+0.0001 && xp[1]>0.3-0.0001  && xp[1]< 0.7+0.0001&& xp[2]< 0.4+0.0001 ) {mat_gam=4;}
   if (mat_gam==2){ bc_flag[0]=1;bc_Neum[0]=1; }
   if (mat_gam==4){ bc_flag[0]=1;bc_Neum[0]=3 ;}
   if (bc_gam==1000)   { bc_flag[0]=1;bc_Neum[0]=5;}
#endif
  return;
}

void MGSolDS::ic_read(
  int /*bc_gam*/,
  int /*mat_gam*/,
  double /*xp*/[],
  double u_value[]
) {// =======================================
#if DIMENSION==2 // --------------  2D //--------------------------
   u_value[0] = 0.;
//    u_value[1] = 0.;
#else  // //-----------------------//------------------------------------
  // xp[]=(xp,yp,zp) u_value[]=(u,v,w,p)
  u_value[0]= ((_dir==0)? 0.:0.);
//   u_value[0]=0.; 
//   u_value[1] =1.;
//   u_value[2] = 2.;
#endif // //-------------------------------------------
  return;
}

// ========================================
/// This function  defines the boundary conditions for the NS system:
void MGSolDS::bc_read(
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
// ===================================
//   double ILref = 1./_lref;

#if DIMENSION==2  // ----- 2D boundary conditions ----------
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
// bc_flag[0]=-1;bc_Neum[0]=-1;
//     bc_flag[0]=0;     bc_Neum[0]=0 ;
  
//     //box 
//    if (xp[0]<1-1.e-6){
//      if (_dir==0) {bc_flag[0]=0;bc_Neum[0]=1; }
//      if (_dir==1) {bc_flag[0]=0;bc_Neum[0]=0; }
//      if (xp[0]<1.e-6){   bc_flag[0]=0;bc_Neum[0]=0; }
//   }
//    if (xp[0]>1+1.e-6){ bc_flag[0]=0;bc_Neum[0]=3 ;}
//    if (fabs(xp[0]-1.)<1.e-6)   { bc_flag[0]=0;bc_Neum[0]=5;}
// //endbox
    //box 
   if (mat_gam==2){
 {bc_flag[0]=0;bc_Neum[0]=0; }
  }
   if (mat_gam==4){ bc_flag[0]=0;bc_Neum[0]=2 ;}
   if (bc_gam==1000)   { bc_flag[0]=0;bc_Neum[0]=4;}
//endbox

//   //turek 
//     if (mat_gam==2){
//  if(bc_gam==100 || bc_gam==200  ) {bc_flag[0]=0; bc_Neum[0]=0 ;}  
//     } else {
//       if(bc_gam==100 || bc_gam==200  ) {bc_flag[0]=0; bc_Neum[0]=2 ;}  
//     }
// 
//    if(bc_gam==100){  bc_flag[0]=0;     bc_Neum[0]=0 ;}
//     if(bc_gam==120){  bc_flag[0]=0;     bc_Neum[0]=0 ;}
//     if(bc_gam==140){  bc_flag[0]=0;     bc_Neum[0]=0 ;
// //       std::cout<<"heyyyyyyyy";
//     }
//     if(bc_gam==1000){
//   bc_flag[0]=0; bc_Neum[0]=4 ;
// }
// //turek
  #endif  // //----------------------//---------------------------

 
 #if DIMENSION==3 
   mat_gam=2;
  if(xp[0]>0.2-0.0001  && xp[0]<0.29+0.0001 && xp[1]>0.3-0.0001  && xp[1]< 0.7+0.0001&& xp[2]< 0.4+0.0001 ) {mat_gam=4;}
   if (mat_gam==2){
 {bc_flag[0]=0;bc_Neum[0]=0; }
  }
   if (mat_gam==4){ bc_flag[0]=0;bc_Neum[0]=2 ;}
   if (bc_gam==1000)   { bc_flag[0]=0;bc_Neum[0]=4;}
 #endif

  return;
 }
#endif
