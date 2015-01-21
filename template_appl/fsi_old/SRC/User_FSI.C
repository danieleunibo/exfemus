#include "Equations_conf.h"

// ===============================
#ifdef FSI_EQUATIONS
// ==============================

#include "MGSclass_conf.h"

// config file -------------
#include "MGGeomEl.h"
#include "MGFE_conf.h"
#include "Printinfo_conf.h"
// local class include -----------
#include "MGSolverFSI.h"
// #include "MGSolverT.h"

// local include
#include "MGSystem.h"
#include "MGEquationsSystem.h"  // Equation map class
#include "MGMesh.h"

// local alg lib
#include "dense_matrixM.h"
#include "sparse_matrixM.h"
#include "dense_vectorM.h"
#include "numeric_vectorM.h"
#include "linear_solverM.h"
#include <string.h>

// ========================================
/// This function  defines the boundary conditions for the NS system:
void MGSolFSI::bc_intern_read(
  int bc_gam,
  int mat_gam,
  double xp[],   // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], // normal
  int bc_flag[]  // boundary condition flag
)  // ===================================
{


#if DIMENSION==2
// // // Pressure box
//   if (xp[0]<1-1.e-6) { // fluid
//     bc_flag[0]=1;     bc_Neum[0]=1 ;
//     bc_flag[1]=1;     bc_Neum[1]=1 ;
//     bc_Neum[2]=1;    // pressure flag
//   } else {
//     if (fabs(xp[0]-1.)<1.e-6)   {  // interface boundary
//       bc_flag[0]=1;        bc_Neum[0]=5;
//       bc_flag[1]=1;      bc_Neum[1]=5;
//       bc_Neum[2]=5;  // pressure flag
//     } else { //solid
//       bc_flag[0]=1;         bc_Neum[0]=3 ;
//       bc_flag[1]=1;      bc_Neum[1]=3 ;
//       bc_Neum[2]=4;  // pressure flag
//     }
//   }
//    if (fabs(xp[0]-1.)<1.e-6)   {  // interface boundary
//       bc_flag[0]=1;        bc_Neum[0]=5;
//       bc_flag[1]=1;      bc_Neum[1]=5;
//       bc_Neum[2]=5;  // pressure flag
//     }
  
  
  // bar.neu
  if (mat_gam==2) { // fluid
    bc_flag[0]=1;     bc_Neum[0]=1 ;
    bc_flag[1]=1;     bc_Neum[1]=1 ;
    bc_Neum[2]=1;    // pressure flag
  } else {
    if (bc_gam==1000)   {  // interface boundary
      bc_flag[0]=1;        bc_Neum[0]=5;
      bc_flag[1]=1;      bc_Neum[1]=5;
      bc_Neum[2]=5;  // pressure flag
    } else { //solid
      bc_flag[0]=1;         bc_Neum[0]=3 ;
      bc_flag[1]=1;      bc_Neum[1]=3 ;
      bc_Neum[2]=2;  // pressure flag
    }
  }
   if (bc_gam==1000)   {  // interface boundary
      bc_flag[0]=1;        bc_Neum[0]=5;
      bc_flag[1]=1;      bc_Neum[1]=5;
      bc_Neum[2]=5;  // pressure flag
    }
  
// // // // begin vlav
// //   if (xp[0]<1-1.e-6 ||xp[0]> 1.01 +1.e-6 || (xp[1]<0.51 -1.e-6 && xp[1]>0.49 +1.e-6) ) { // fluid
//   if(mat_gam==2){
//     bc_flag[0]=1;     bc_Neum[0]=1 ;
//     bc_flag[1]=1;     bc_Neum[1]=1 ;
//     bc_Neum[2]=1;    // pressure flag
// //      std::cout<<"fluid";
//   } else {
// //     if (xp[0]<1+1.e-6 ||xp[0]> 1.01 -1.e-6 || (xp[1]<0.51 +1.e-6 && xp[1]>0.49 -1.e-6) )  {  // interface boundary
//       if(bc_gam==1000){
//       bc_flag[0]=1;        bc_Neum[0]=5;
//       bc_flag[1]=1;      bc_Neum[1]=5;
//       bc_Neum[2]=5;  // pressure flag
// //       std::cout<<"iterface";
//     } else { //solid
//       bc_flag[0]=1;         bc_Neum[0]=3 ;
//       bc_flag[1]=1;      bc_Neum[1]=3 ;
//       bc_Neum[2]=4;  // pressure flag
// //             std::cout<<"solid";
//     }
//   }
//   // // // end vlav
  
  
  
//   // // // begin turek_1
//   if(mat_gam==2){
//     bc_flag[0]=1;     bc_Neum[0]=1 ;
//     bc_flag[1]=1;     bc_Neum[1]=1 ;
//     bc_Neum[2]=1;    // pressure flag
// //      std::cout<<"fluid";
//   } else {
//  //solid
//       bc_flag[0]=1;         bc_Neum[0]=3 ;
//       bc_flag[1]=1;      bc_Neum[1]=3 ;
//       bc_Neum[2]=2;  // pressure flag
//     }
//    
//   if(bc_gam==1000){
//       bc_flag[0]=1;        bc_Neum[0]=5;
//       bc_flag[1]=1;      bc_Neum[1]=5;
//       bc_Neum[2]=5;  // pressure flag
// //       std::cout<<"iterface";
//     }
//   // // //turek_1
#endif
#if DIMENSION==3
  // box3d.neu
  mat_gam=2;
  if(xp[0]>0.2-0.0001  && xp[0]<0.29+0.0001 && xp[1]>0.3-0.0001  && xp[1]< 0.7+0.0001&& xp[2]< 0.4+0.0001 ) {mat_gam=4;}
  if (mat_gam==2) { // fluid
    bc_flag[0]=1;     bc_Neum[0]=1 ;
    bc_flag[1]=1;     bc_Neum[1]=1 ;
    bc_Neum[2]=1;    // pressure flag
  } else {
    if (bc_gam==1000)   {  // interface boundary
      bc_flag[0]=1;        bc_Neum[0]=5;
      bc_flag[1]=1;      bc_Neum[1]=5;
      bc_Neum[2]=5;  // pressure flag
    } else { //solid
      bc_flag[0]=1;         bc_Neum[0]=3 ;
      bc_flag[1]=1;      bc_Neum[1]=3 ;
      bc_Neum[2]=4;  // pressure flag
    }
  }
   if (bc_gam==1000)   {  // interface boundary
      bc_flag[0]=1;        bc_Neum[0]=5;
      bc_flag[1]=1;      bc_Neum[1]=5;
      bc_Neum[2]=5;  // pressure flag
    }
#endif

  return;
}








// =========================================================
/// This function generates the initial conditions for the NS system:
void MGSolFSI::ic_read(
  int bc_gam,
  int mat_gam,
  double xp[],
  double u_value[]
)  // =======================================
{
//xp[] is the NON-DIMENSIONAL node coordinate
  double ILref = 1./_lref;


// ===================================================================
#if DIMENSION==2 // --------------  2D //--------------------------
  // xp[]=(xp,yp) u_value[]=(u,v,p)

#if FSI_EQUATIONS==2
  if(_dir==0)  { u_value[0] = 0.; }
  if(_dir==1)   {
    u_value[0] = 0.;
//     if(xp[1]<0.01) {
    if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) { u_value[0] =0.*0.777+ 0*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref; }
//     }
  }
#else
//   u_value[0] = -cos(2.*3.14159265359*xp[0])*sin(2.*3.14159265359*xp[1]);
//   u_value[1] =  sin(2.*3.14159265359*xp[0])*cos(2.*3.14159265359*xp[1]);

  u_value[0] = 0.;    u_value[1] = -0.0; u_value[2] = 0.;
 if(xp[0]<0.0000001 )   u_value[0] =0.1;
//     if(mat_gam==4 ) u_value[1] =-0.0001;
   if(bc_gam==100) {
//       u_value[0] =20/0.1681*xp[1]*(0.41-xp[1]);
  }

#endif
#if FSI_EQUATIONS==1 //coupled
  double pref = _refvalue[DIMENSION];
  u_value[DIMENSION] =0.;//0.*(LYE*ILref - xp[1])/pref;
//      u_value[DIMENSION] = -0.25*0.0018*(cos(4.*3.14159265359*xp[0]) + cos(4.*3.14159265359*xp[1]))/pref;
#endif

#endif  // //-----------------------//------------------------------------
#if DIMENSION==3 // --------------  2D //--------------------------
 u_value[0] = 0.;    u_value[1] = -0.0; u_value[2] = 0.; u_value[3] = 0.;
  #endif
  
// ===================================================================
#if DIMENSION==3  // --------------  3D //--------------------------

  u_value[0]= 0.;
      if(xp[0] < 0+0.0001)  u_value[0]= 5.5;
  u_value[1] =0.;
  u_value[2] =0.;
  u_value[3] =0.;
#endif // //-------------------------------------------

  return;
}

// ========================================
/// This function  defines the boundary conditions for the NS system:
void MGSolFSI::bc_read(
  int bc_gam,
  int mat_gam,
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], // normal
  int bc_flag[]         // boundary condition flag
)  // ===================================
{
  //     0 ->  single component
  //     +4 ->  nonhomogeneous
  //     +2 -> tg
  //     +1 ->  normal
// ===================================
//   double ILref = 1./_lref;
#if DIMENSION==2  // ----- 2D boundary conditions ----------
//     bc_flag[0]=0;     bc_Neum[0]=1 ;
//     bc_flag[1]=0;     bc_Neum[1]=1 ;
//     bc_Neum[2]=1;    // pressure flag 
    
// // //     pressure box
//   if (xp[0]<1-1.e-6) { // fluid
//     bc_flag[0]=0;     bc_Neum[0]=1 ;
//     bc_flag[1]=0;     bc_Neum[1]=0 ;
//     bc_Neum[2]=1;    // pressure flag
//     
//       if (xp[0]< 1.e-6){   bc_flag[0]=4;     bc_Neum[0]=1 ;}
//     
//   } else {
// 
//       bc_flag[0]=0;         bc_Neum[0]=2 ;
//       bc_flag[1]=0;      bc_Neum[1]=2 ;
//       bc_Neum[2]=4;  // pressure flag
//         if (xp[1]> 1- 1.e-6||xp[1]< 1.e-6){ 
// 	   bc_flag[0]=0;      bc_Neum[0]=3 ;
// 	  bc_flag[1]=0;      bc_Neum[1]=3 ;
// 	  
// 	}
//     
//   }
//    if (fabs(xp[0]-1.)<1.e-6)   {  // interface boundary
//       bc_flag[0]=0;        bc_Neum[0]=5;
//       bc_flag[1]=0;      bc_Neum[1]=4;
//       bc_Neum[2]=5;  // pressure flag
//     }
  // //     pressure box  
  
  
  
  //bar.neu
    bc_flag[0]=0;     bc_Neum[0]=0 ;
    bc_flag[1]=0;     bc_Neum[1]=0 ;
    bc_Neum[2]=0;    // pressure flag
    
      if (xp[0]< 1.e-6){   bc_flag[0]=4;     bc_Neum[0]=0 ;}
      if (xp[0]> 2.- 1.e-6){   bc_flag[0]=0;     bc_Neum[0]=1 ;}
    
if (xp[0]>1. &&xp[0]<1.05 && xp[1]<0.5){
      bc_flag[0]=0;         bc_Neum[0]=2 ;
      bc_flag[1]=0;      bc_Neum[1]=2 ;
      bc_Neum[2]=2;  // pressure flag
  
}
//         if (xp[1]< 1.e-6){ 
// 	   bc_flag[0]=0;      bc_Neum[0]=3 ;
// 	  bc_flag[1]=0;      bc_Neum[1]=3 ;
	  
    
  
   if (bc_gam==1000)   {  // interface boundary
      bc_flag[0]=0;        bc_Neum[0]=4;
      bc_flag[1]=0;      bc_Neum[1]=4;
      bc_Neum[2]=5;  // pressure flag
    }
  
  
  
  
  
    
// // valv
//   if (mat_gam==2) { // fluid
//     bc_flag[0]=0;     bc_Neum[0]=1 ;
//     bc_flag[1]=0;     bc_Neum[1]=1 ;
//     bc_Neum[2]=1;    // pressure flag 
//     if(bc_gam==300) {  bc_flag[1]=4;     bc_Neum[1]=0 ;   bc_flag[0]=0;     bc_Neum[0]=0 ; bc_Neum[2]=0;} //inlet
//     if(bc_gam==100) {//lateral symmetry axis
//       bc_flag[0]=0;     bc_Neum[0]=0 ;
//       bc_flag[1]=0;     bc_Neum[1]=1 ;
//       bc_Neum[2]=0;    // pressure flag
//     }
//   if(bc_gam==200){ //outlet 
//     bc_flag[0]=0;     bc_Neum[0]=10 ;
//     bc_flag[1]=0;     bc_Neum[1]=1 ;
//   }
//   }
//   else {
//     if (bc_gam=1000)  {  // interface boundary
//       bc_flag[0]=0;        bc_Neum[0]=5;
//       bc_flag[1]=0;     bc_Neum[1]=4;
//       bc_Neum[2]=5;  // pressure flag
//     } else { //solid
//       bc_flag[0]=0;      bc_Neum[0]=3 ;
//       bc_flag[1]=0;      bc_Neum[1]=3 ;
//       bc_Neum[2]=4;  // pressure flag
// 
//       if(bc_gam==300) { //inlet
//         bc_flag[0]=0;        bc_Neum[0]=3;
//         bc_flag[1]=0;     bc_Neum[1]=2;
//         bc_Neum[2]=4;  // pressure flag
//       }
//         if(bc_gam==200) { //outlet
//         bc_flag[0]=0;        bc_Neum[0]=3;
//         bc_flag[1]=0;     bc_Neum[1]=2;
//         bc_Neum[2]=4;  // pressure flag
//       }
//     }
//   }
// // endvalv


// // // turek_1
//   if (mat_gam==2) { // fluid
//     bc_flag[0]=0;     bc_Neum[0]=1 ;
//     bc_flag[1]=0;     bc_Neum[1]=1 ;
//     bc_Neum[2]=1;    // pressure flag 
//     if(bc_gam==100) {  bc_flag[0]=4;     bc_Neum[0]=0 ;   bc_flag[1]=0;     bc_Neum[1]=0 ; bc_Neum[2]=0;} //inlet
//     if(bc_gam==200){ //cyrcle 
//     bc_flag[0]=0;     bc_Neum[0]=0 ;
//     bc_flag[1]=0;     bc_Neum[1]=0 ;
//    bc_Neum[2]=0;
//   }
//     if(bc_gam==140){ //outlet 
//     bc_flag[0]=0;     bc_Neum[0]=1 ;
//     bc_flag[1]=0;     bc_Neum[1]=0 ;
//     bc_Neum[2]=0;
//   }
//     if(bc_gam==120 ||(xp[1]<0.00001 || xp[1]>0.41-0.0001)) {//lateral 
//       bc_flag[0]=0;     bc_Neum[0]=0 ;
//       bc_flag[1]=0;     bc_Neum[1]=0 ;
//       bc_Neum[2]=0;    // pressure flag
//     }
// 
// 
//   }
//   else { 
// 
//       bc_flag[0]=0;      bc_Neum[0]=2 ;
//       bc_flag[1]=0;      bc_Neum[1]=2 ;
//       bc_Neum[2]=2;  // pressure flag
// 
// //        if(bc_gam==200) { //tube
// // 	 std::cout<<"inner bc  solid ";
// //         bc_flag[0]=0;        bc_Neum[0]=2;
// //         bc_flag[1]=0;     bc_Neum[1]=2;
// //         bc_Neum[2]=2;  // pressure flag
// //       }
//     
//   }
//          if(bc_gam==1000) { //tube
// 	 std::cout<<"inner bc  solid ";
//         bc_flag[0]=0;        bc_Neum[0]=4;
//         bc_flag[1]=0;     bc_Neum[1]=4;
//         bc_Neum[2]=4;  // pressure flag
//       }
// // // endturek1
#endif


#if DIMENSION==3

  mat_gam=2;
  if(xp[0]>0.2-0.0001  && xp[0]<0.29+0.0001 && xp[1]>0.3-0.0001  && xp[1]< 0.7+0.0001&& xp[2]< 0.4+0.0001 ) {mat_gam=4;}
  //bar.neu
    if (mat_gam==2) { // fluid
    bc_flag[0]=0;     bc_Neum[0]=0 ;
    bc_flag[1]=0;     bc_Neum[1]=0 ;
    bc_flag[2]=0;     bc_Neum[2]=0 ;
    bc_Neum[3]=1;    // pressure flag
    
      if (xp[0]< 1.e-6){   bc_flag[0]=4;     bc_Neum[0]=0 ;}
    
  } else {

      bc_flag[0]=0;         bc_Neum[0]=2 ;
      bc_flag[1]=0;      bc_Neum[1]=2 ;
          bc_flag[2]=0;      bc_Neum[2]=2 ;
      bc_Neum[3]=4;  // pressure flag
//         if (xp[1]< 1.e-6){ 
// 	   bc_flag[0]=0;      bc_Neum[0]=3 ;
// 	  bc_flag[1]=0;      bc_Neum[1]=3 ;
	  
	}
    
  
   if (bc_gam==1000)   {  // interface boundary
      bc_flag[0]=0;        bc_Neum[0]=4;
      bc_flag[1]=0;      bc_Neum[1]=4;
      bc_flag[2]=0;      bc_Neum[2]=4;
      bc_Neum[3]=5;  // pressure flag
    }



#endif


}


#ifdef FSIP_EQUATIONS

// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
// =================================================
void MGSolFSIP::ic_read(double xp[],double u_value[])
{
// xp[]=(xp,yp) u_value[]=(u,v,p)
  u_value[0] =0.;
}

// ========================================================
/// This function  defines the boundary conditions for the energy equation:
void MGSolFSIP::bc_read(double xp[],int bc_Neum[],int bc_flag[])
{
  // =================================================
  const double Lref = _mgphys.get_par("Lref");
  double ILref = 1./Lref;
#if DIMENSION==2
// xp[]=(xp,yp) bc_flag[T]=0-Dirichlet 1-Neumann
  //    boundary conditions box

//   if(xp[0]< LXB*ILref+BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }   // left
//
//   if(xp[1]< LYB*ILref+BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }   // bottom
//   if(xp[0]> LXE*ILref-BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=0;  }   // right
//   if(xp[1]> LYE*ILref-BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=0;  }   // top

  bc_flag[0]=0; bc_Neum[0]=1;
  if(xp[0] > 2.5 -BDRY_TOLL ) {bc_flag[0]=0; bc_Neum[0]=0;} //right fluid
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

#endif

#endif // ENDIF FSI_EQUATIONS==0





