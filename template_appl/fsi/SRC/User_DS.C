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
#if DIMENSION == 2
     if (mat_gam==2){ bc_flag[0]=1;bc_Neum[0]=1; }
   if (mat_gam==4){ bc_flag[0]=1;bc_Neum[0]=3 ;}
   if (bc_gam==1000)   { bc_flag[0]=1;bc_Neum[0]=5;}
#endif
#if DIMENSION==3
// // // // bar.neu box
   if (mat_gam==2){ bc_flag[0]=1;bc_Neum[0]=1; }
   if (mat_gam==4){ bc_flag[0]=1;bc_Neum[0]=3 ;}
   if (bc_gam==1000)   { bc_flag[0]=1;bc_Neum[0]=5;}
   
   
//    if (xp[0]<1-1.e-6){ bc_flag[0]=1;bc_Neum[0]=1; }
//    if (xp[0]>1+1.e-6){ bc_flag[0]=1;bc_Neum[0]=3 ;}
//    if (fabs(xp[0]-1.)<1.e-6)   { bc_flag[0]=1;bc_Neum[0]=5;}
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
     if (bc_gam< 20){
     if (_dir==0) {bc_flag[0]=0;bc_Neum[0]=0; }
     if (_dir==1) {bc_flag[0]=0;bc_Neum[0]=0; }
     if (_dir==2) {bc_flag[0]=0;bc_Neum[0]=0; }
  }
    else if (bc_gam<30) { bc_flag[0]=0;bc_Neum[0]=3 ;}
   else if (bc_gam>900)   { bc_flag[0]=0;bc_Neum[0]=5;}
  #endif  // //----------------------//---------------------------

 
 #if DIMENSION==3 
     //box 
   if (bc_gam< 20){
     if (_dir==0) {bc_flag[0]=0;bc_Neum[0]=0; }
     if (_dir==1) {bc_flag[0]=0;bc_Neum[0]=0; }
     if (_dir==2) {bc_flag[0]=0;bc_Neum[0]=0; }
  }
    else if (bc_gam<30) { bc_flag[0]=0;bc_Neum[0]=3 ;}
   else if (bc_gam>900)   { bc_flag[0]=0;bc_Neum[0]=5;}
//endbox
   
 #endif

  return;
 }
#endif
