// ======================================================================================
// --------------   NAVIER-STOKES system [NS_F] -----------------------------------------
// ======================================================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
// ======================================================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ======================================================================================

// class files --------------------------------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverNS.h"       // Navier-Stokes class header file


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




// ======================================================================================
/// This function generates the initial conditions for the NS system:
void MGSolNS::ic_read(
  int bc_gam,
  int bc_mat,
  double xp[],
  double u_value[]
) {// ===================================================================================
//xp[] is the NON-DIMENSIONAL node coordinate
  double ILref = 1./_lref;

#if DIMENSION==2 // --------------  2D //--------------------------
  
  #if NS_EQUATIONS==2
    if (_dir==0)  { u_value[0] = 0.; }
    if (_dir==1)   {
      u_value[0] = 0.; 
      if (/*xp[0] >LXB*ILref+BDRY_TOLL &&*/ xp[0] <LXE*ILref-BDRY_TOLL) { u_value[0] =1.*0.25+ 0*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref; }
    }  
  #else
  u_value[0] = 0.;
  u_value[1] = 0.5648;
  // cyl
//     if (xp[0] <LXE*ILref-BDRY_TOLL) u_value[1] =-0.25e+3*(xp[0]-LXE*ILref)*(LXE*ILref+xp[0])/_uref;//0.0165; 
    // annulus
//     if (xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) u_value[1] =1.; 
  #endif
  
  #if NS_EQUATIONS==1 //coupled
    double pref = _refvalue[DIMENSION];
    u_value[DIMENSION] =0.*(LYE*ILref - xp[1])/pref;
  #endif

#endif  // //-----------------------//------------------------------------

// ===================================================================

#if DIMENSION==3  // --------------  3D //--------------------------

//   #if NS_EQUATIONS==2
//     if (_dir==0)  { u_value[0] = 0.; }
//     if (_dir==1)   {
//       u_value[0] = 0.;
//       if (xp[1]<0.01) {
// 	if (xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) { u_value[0] =1.*0.777+ 0*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref; }
//       }
//     }
//     if (_dir==2)  { u_value[0] = 0.; }
//   #else
//   // xp[]=(xp,yp,zp) u_value[]=(u,v,w,p)
//   u_value[0]= 0.;
//   u_value[1] =0.;
// //   if (xp[1]<0.01) {
// //     if (xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) { u_value[0] =1.*0.777+ 0*(xp[0] -LXB*ILref)*(LXE*ILref-xp[0])/_uref; }
// //   }
//   u_value[2] = 0.1832;
//   #endif
  
//   #if NS_EQUATIONS==1 //coupled
//   double pref = _refvalue[DIMENSION];
//   u_value[DIMENSION] =0.*(LYE*ILref - xp[1])/pref;
//   #endif

  
  // ====================================================================================
  // ================================= return _map_med;===================================================
  // Test1
  // ====================================================================================
  // ====================================================================================
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
  // Initial conditions
  u_value[0]= 0.;
  u_value[1] =0.;
  u_value[2] = 0.;
  
  #if NS_EQUATIONS==1 //coupled
  double pref = _refvalue[DIMENSION]; u_value[DIMENSION] =0.*(LYE*ILref - xp[1])/pref;
  #endif
  // ====================================================================================
  // ==================================================================================== 
#endif
  return;
}

// ========================================
// This function  defines the boundary conditions for the NS system:
void MGSolNS::bc_read(
  int bc_gam,
  int bc_mat,
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[], 	// normal
  int bc_flag[]         // boundary condition flag
) {// ===================================
  //     0  ->  single component
  //     +4 ->  nonhomogeneous
  //     +2 ->  tg
  //     +1 ->  normal
// ===================================
  double ILref = 1./_lref;

#if DIMENSION==2  // ----- 2D boundary conditions ----------
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann

#if NS_EQUATIONS==2
  if (_dir==0) { //    boundary conditions box
    if (xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } // top
    if (xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } // bottom
    if (xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } // left
    if (xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } // right
  }
  if (_dir==1) { //    boundary conditions box
    if (xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; } // top
    if (xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=4; bc_Neum[0]=0; } // bottom
    if (xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; } // left
    if (xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; } // right
  }

#else


//    boundary conditions box
#if NS_EQUATIONS==1
  bc_flag[2]=0;
#endif  
  // channel 3m long
//   if (xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=0; bc_Neum[1]=1;} //  top
//   if (xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=4; bc_Neum[0]=0; bc_Neum[1]=0;} //  bottom
//   if (xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=2; bc_Neum[0]=0; bc_Neum[1]=1;} //  left
//   if (xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=2; bc_Neum[0]=0; bc_Neum[1]=1;} //  right

// cyl
  if (xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=0; bc_Neum[1]=1;} //  top
  if (xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=0; bc_Neum[1]=1;} //  bottom
  if (xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=0; bc_Neum[1]=1;} //  left
  if (xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=2; bc_Neum[0]=0; bc_Neum[1]=1;} //  right
// annulus
//    if (xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=0; bc_Neum[0]=0; bc_Neum[1]=1;} //  top
//    if (xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=4; bc_Neum[0]=0; bc_Neum[1]=1;} //  bottom
//    if (xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=2; bc_Neum[0]=0; bc_Neum[1]=1;} //  left
//    if (xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_flag[1]=2; bc_Neum[0]=0; bc_Neum[1]=1;} //  right

#endif // //---------------------------------------------

#endif // //---------------------------------------------

#if DIMENSION==3 // -------- 3D boundary conditions ---------
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
  //    boundary conditions box
#if NS_EQUATIONS==2
  if (_dir==0) { //    boundary conditions box
    if (xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; } //  left
    if (xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; } //  right
    if (xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } // top
    if (xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } // bottom
    if (xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } //  left
    if (xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } //  right
  }
  if (_dir==1) { //    boundary conditions box
    if (xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; } //  left
    if (xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; } //  right
    if (xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=1; } // top
    if (xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=4; bc_Neum[0]=0; } // bottom
    if (xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } //  left
    if (xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } //  right
  }
  if (_dir==2) { //    boundary conditions box
    if (xp[1] > LYE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } // top
    if (xp[1] < LYB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } // bottom
    if (xp[0] < LXB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } //  left
    if (xp[0] > LXE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } //  right
    if (xp[2] < LZB*ILref + BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } //  left
    if (xp[2] > LZE*ILref - BDRY_TOLL) {bc_flag[0]=0; bc_Neum[0]=0; } //  right
  }

#else

//   bc_flag[3]=0;
//   if (xp[1]< LYB*ILref + BDRY_TOLL){bc_flag[0]=0;bc_flag[1]=0;bc_flag[2]=0; bc_Neum[0]=1;bc_Neum[1]=1;bc_Neum[2]=1;} //  INLET
//   if (xp[1]> LYE*ILref - BDRY_TOLL){bc_flag[0]=0;bc_flag[1]=0;bc_flag[2]=0; bc_Neum[0]=1;bc_Neum[1]=1;bc_Neum[2]=1;} //  OUTLET
//   if (xp[0]< LXB*ILref + BDRY_TOLL){bc_flag[0]=0;bc_flag[1]=0;bc_flag[2]=0; bc_Neum[0]=1;bc_Neum[1]=1;bc_Neum[2]=1;}
//   if (xp[0]> LXE*ILref - BDRY_TOLL){bc_flag[0]=0;bc_flag[1]=0;bc_flag[2]=0; bc_Neum[0]=1;bc_Neum[1]=1;bc_Neum[2]=1;}
//   if (xp[2]< LZB*ILref + BDRY_TOLL){bc_flag[0]=0;bc_flag[1]=2;bc_flag[2]=0; bc_Neum[0]=0;bc_Neum[1]=1;bc_Neum[2]=0;} //current
//   if (xp[2]> LZE*ILref - BDRY_TOLL){bc_flag[0]=0;bc_flag[1]=2;bc_flag[2]=0; bc_Neum[0]=0;bc_Neum[1]=1;bc_Neum[2]=0;} //current


// // square bundle
     double p=0.00615;  // 1.5
//    double p=0.00533; // 1.3
//   double p=0.005;   // 1.2
  bc_flag[3]=0;
  bc_flag[0]=0;bc_flag[1]=0;bc_flag[2]=0; bc_Neum[0]=1;bc_Neum[1]=1;bc_Neum[2]=1;
  
  
   if (xp[1]*xp[1]+(xp[0]-p)*(xp[0]-p)< (0.0041+.5e-6)*(0.0041+.5e-6)){
// //      bc_flag[0]=0;bc_flag[1]=0;bc_flag[2]=2; bc_Neum[0]=0;bc_Neum[1]=0;bc_Neum[2]=1;}
     bc_flag[0]=0;bc_Neum[0]=0;bc_flag[1]=0;bc_Neum[1]=0;bc_flag[2]=2;bc_Neum[2]=1;}

   if(xp[1]< BDRY_TOLL)  {bc_flag[1]=0;bc_Neum[1]=0;}
   if(xp[0]< BDRY_TOLL)  {bc_flag[0]=0;bc_Neum[0]=0;}
//      if(xp[1]> 0.00533-xp[0]  -BDRY_TOLL)  {bc_flag[0]=2;bc_Neum[0]=1;bc_flag[1]=0;bc_Neum[1]=0;}
//      
     
     
 //  tri bundle  
//    bc_flag[3]=0;
//    bc_flag[0]=0;bc_flag[1]=0;bc_flag[2]=0; bc_Neum[0]=1;bc_Neum[1]=1;bc_Neum[2]=1;
// //  
//    if(xp[1]< BDRY_TOLL)  {bc_flag[1]=0;bc_Neum[1]=1;}                                                   //down side
//    if(xp[0]< BDRY_TOLL)  {bc_flag[0]=0;bc_Neum[0]=1;}						        //left side
//    if (xp[1]*xp[1]+(xp[0]-0.00533)*(xp[0]-0.00533)< 0.0041*0.0041+BDRY_TOLL){                           //rod arc   P/D=1.3
//      bc_flag[0]=0;bc_Neum[0]=0;bc_flag[1]=0;bc_Neum[1]=0;bc_flag[2]=2;bc_Neum[2]=1;}
//      if (xp[1]*xp[1]+(xp[0]-0.00574)*(xp[0]-0.00574)< 0.0041*0.0041+BDRY_TOLL){                           //rod arc    P/D=1.4
//      bc_flag[0]=0;bc_Neum[0]=0;bc_flag[1]=0;bc_Neum[1]=0;bc_flag[2]=2;bc_Neum[2]=1;}
//      if (xp[1]*xp[1]+(xp[0]-0.00615)*(xp[0]-0.00615)< 0.0041*0.0041+BDRY_TOLL){                           //rod arc    P/D=1.5
//      bc_flag[0]=0;bc_Neum[0]=0;bc_flag[1]=0;bc_Neum[1]=0;bc_flag[2]=2;bc_Neum[2]=1;}

#endif



  // =================================================
  // =================================================
  // Test1 3D
  // =================================================
  // =================================================
  // =================================================
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
  //    boundary conditions box Test1
   int imesh=atoi(_mgutils.get_file("MESHNUMBER").c_str());
   bc_Neum[0]=1; bc_Neum[1]=1;  bc_Neum[2]=1; bc_Neum[3]=1;
   bc_flag[0]=0; bc_flag[1]=0;  bc_flag[2]=0; bc_flag[2]=0; 
   
   if(imesh==1){
   // mesh1 simple_hex27.med
   if( bc_gam==100){
     bc_Neum[0]=0; bc_Neum[1]=0;  bc_Neum[2]=0;    
     bc_flag[0]=0; bc_flag[1]=0;  bc_flag[2]=0; 
  }   // lateral
   if( bc_gam==50){     
     bc_Neum[0]=0; bc_Neum[1]=0;  bc_Neum[2]=0;    
     bc_flag[0]=0; bc_flag[1]=0;  bc_flag[2]=4; 
  }   // bottom 
    if( bc_gam==11){     
     bc_Neum[0]=0; bc_Neum[1]=0;  bc_Neum[2]=0;    
     bc_flag[0]=0; bc_flag[1]=0;  bc_flag[2]=0; 
  }   // lateral
   if( bc_gam==211){ 
     bc_Neum[0]=0; bc_Neum[1]=0;  bc_Neum[2]=1;    
     bc_flag[0]=0; bc_flag[1]=0;  bc_flag[2]=4; 
  }   // outlet 
  
       if (xp[0] < 0. +BDRY_TOLL || xp[0] > 1.-BDRY_TOLL ||
         xp[1] < 0. +BDRY_TOLL || xp[1] > 1.-BDRY_TOLL) {
       bc_Neum[0]=0; bc_Neum[1]=0;  bc_Neum[2]=0;    
       bc_flag[0]=0; bc_flag[1]=0;  bc_flag[2]=0;
     }
   }  
 
   if(imesh==2){
   // mesh 2  simple_hex27_2.med
     if( bc_gam==211){
         bc_Neum[0]=0; bc_Neum[1]=0;  bc_Neum[2]=0;    
         bc_flag[0]=0; bc_flag[1]=0;  bc_flag[2]=4;
    } // inlet Dirichtlet coupled  
     if( bc_gam==200){ 
         bc_Neum[0]=0; bc_Neum[1]=0;  bc_Neum[2]=0;    
         bc_flag[0]=0; bc_flag[1]=0;  bc_flag[2]=0;
    }  // lateral D 
     if( bc_gam==150){     
     bc_Neum[0]=0; bc_Neum[1]=0;  bc_Neum[2]=1;    
     bc_flag[0]=0; bc_flag[1]=0;  bc_flag[2]=0; 
  }   // outlet 
   }
   // =================================================
   // =================================================
  
#endif  //  end DIMENSION ==3 ---------------------------

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

  if (xp[0]< LXB*ILref+BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }  // left
  if (xp[0]> LXE*ILref-BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }  // right
  if (xp[1]< LYB*ILref+BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }  // bottom
  if (xp[1]> LYE*ILref-BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=0;  }  // top
#endif

#if DIMENSION==3
// =================================================
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
  //    boundary conditions box
  if (xp[0] < LXB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
  if (xp[0] > LXE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
  if (xp[1] < LYB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
  if (xp[2] < LZB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
  if (xp[2] > LZE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
  if (xp[1] > LYE*ILref-BDRY_TOLL) {    bc_Neum[0]=0;    bc_flag[0]=0;  }
#endif
  return;
} // end boundary conditions ==========================





#endif // ENDIF NS_EQUATIONS==0
#endif  //ENDIF NS_EQUATIONS

