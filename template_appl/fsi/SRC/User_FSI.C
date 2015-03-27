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
    if (xp[0]<(LXE-LXB)/2-BDRY_TOLL) { // fluid
    bc_flag[0]=1;     bc_Neum[0]=1 ;
    bc_flag[1]=1;     bc_Neum[1]=1 ;
                      bc_Neum[2]=1;   
  }  else { //solid
      bc_flag[0]=1;      bc_Neum[0]=3 ;
      bc_flag[1]=1;      bc_Neum[1]=3 ;
                         bc_Neum[2]=3;  // pressure flag 
    }
   if (fabs(xp[0]-(LXE-LXB)/2 )<BDRY_TOLL)   {  // interface boundary
      bc_flag[0]=0;   bc_Neum[0]=5;
      bc_flag[1]=0;   bc_Neum[1]=5;
                      bc_Neum[2]=5;
    }
#endif
#if DIMENSION==3
  if (mat_gam==2) { // fluid
    bc_flag[0]=1;     bc_Neum[0]=1 ;
    bc_flag[1]=1;     bc_Neum[1]=1 ;
    bc_flag[2]=1;    bc_Neum[2]=1;   
   bc_flag[3]=0;  bc_Neum[3]=1;    // pressure flag
  }  else { //solid
      bc_flag[0]=1;         bc_Neum[0]=3 ;
      bc_flag[1]=1;      bc_Neum[1]=3 ;
      bc_flag[2]=1;      bc_Neum[2]=3;  // pressure flag 
   bc_flag[3]=0;  bc_Neum[3]=3;    // pressure flag
    }
   if (bc_gam==1000)   {  // interface boundary
      bc_flag[0]=0;        bc_Neum[0]=4;
      bc_flag[1]=0;      bc_Neum[1]=5;
      bc_flag[2]=0;      bc_Neum[2]=4;
     bc_flag[3]=0;  bc_Neum[3]=5;    // pressure flag
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
    u_value[0]=0.; u_value[1]=0.; u_value[2]=0.;
    
    if (bc_gam==10) { u_value[0]=10.; }
#if DIMENSION==2
 if (xp[0]<LXB+BDRY_TOLL) u_value[0]=0.8;
#endif
#if DIMENSION==3
    u_value[3]=0;
#endif
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
    if (xp[0]<LXB+BDRY_TOLL) { //inlet
    bc_flag[0]=4;     bc_Neum[0]=0 ;
    bc_flag[1]=0;     bc_Neum[1]=0 ;
                      bc_Neum[2]=1 ;
    }
    if (xp[0]<(LXE-LXB)/2-BDRY_TOLL){//liquid border
      if (xp[1]<LYB+BDRY_TOLL || xp[1]>LYE-BDRY_TOLL ) { 
	bc_flag[0]=0;     bc_Neum[0]=0 ;
	bc_flag[1]=0;     bc_Neum[1]=0 ;
			  bc_Neum[2]=1 ;
	}
    }else if (xp[0]>(LXE-LXB)/2+BDRY_TOLL){
      	bc_flag[0]=0;     bc_Neum[0]=3 ;//solid border
	bc_flag[1]=0;     bc_Neum[1]=3 ;
			  bc_Neum[2]=3 ;
      if (xp[0]>LXE-BDRY_TOLL) { //solid end
    bc_flag[0]=0;     bc_Neum[0]=2 ;
    bc_flag[1]=0;     bc_Neum[1]=2 ;
                      bc_Neum[2]=3 ;
    }
    }

    if(fabs(xp[0]-(LXE-LXB)/2 )<BDRY_TOLL)  { //interface
    bc_flag[0]=0;     bc_Neum[0]=5 ;
    bc_flag[1]=0;     bc_Neum[1]=4 ;
                      bc_Neum[2]=5 ;
    }
#endif


#if DIMENSION==3


// pressbox
    if (bc_gam==11 || bc_gam==12 ||bc_gam==13 || bc_gam==14  ) { 
    bc_flag[0]=0;     bc_Neum[0]=0 ;
    bc_flag[1]=0;     bc_Neum[1]=1 ;
    bc_flag[2]=0;     bc_Neum[2]=0 ;
      bc_Neum[3]=0 ; bc_Neum[3]=1 ;
    }
     if (bc_gam==15) { 
    bc_flag[0]=0;     bc_Neum[0]=0 ;
    bc_flag[1]=4;     bc_Neum[1]=1 ;
    bc_flag[2]=0;     bc_Neum[2]=0 ;
     bc_Neum[3]=0 ; bc_Neum[3]=1 ;
    }
    if (bc_gam==21 ) { 
    bc_flag[0]=0;     bc_Neum[0]=3 ;
    bc_flag[1]=0;     bc_Neum[1]=3 ;
    bc_flag[2]=0;     bc_Neum[2]=3 ;
        bc_Neum[3]=0 ; bc_Neum[3]=3 ;
    } 
    if (bc_gam==23) { 
    bc_flag[0]=0;     bc_Neum[0]=3 ;
    bc_flag[1]=0;     bc_Neum[1]=3 ;
    bc_flag[2]=0;     bc_Neum[2]=3 ;
        bc_Neum[3]=0 ; bc_Neum[3]=3 ;
    } 
    if (bc_gam==22 ) { 
    bc_flag[0]=0;     bc_Neum[0]=3 ;
    bc_flag[1]=0;     bc_Neum[1]=3 ;
    bc_flag[2]=0;     bc_Neum[2]=3 ;
      bc_Neum[3]=0 ; bc_Neum[3]=3 ;
    } 
    if (bc_gam==24 ) { 
    bc_flag[0]=0;     bc_Neum[0]=3 ;
    bc_flag[1]=0;     bc_Neum[1]=3 ;
    bc_flag[2]=0;     bc_Neum[2]=3 ;
      bc_Neum[3]=0 ; bc_Neum[3]=3 ;
    } 
    if (bc_gam==25) { 
    bc_flag[0]=0;     bc_Neum[0]=2 ;
    bc_flag[1]=0;     bc_Neum[1]=2 ;
    bc_flag[2]=0;     bc_Neum[2]=2 ;
      bc_Neum[3]=0 ; bc_Neum[3]=3 ;
    } 
        if (bc_gam==1000) {
    bc_flag[0]=0;     bc_Neum[0]=4 ;
    bc_flag[1]=0;     bc_Neum[1]=5 ;
    bc_flag[2]=0;     bc_Neum[2]=4 ; 
       bc_Neum[3]=0 ; bc_Neum[3]=5 ;
	} 
// end pressbox


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





