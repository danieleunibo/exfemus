// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef NS_EQUATIONS
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

// class files --------------------------------------------------
#include "MGSolverNS.h"       // Navier-Stokes class header file


// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation   
#include "Printinfo_conf.h"  // Print options

// local Femus class include -----------------------------------
#include "MGMesh.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsMap.h"  // Equation map class

// standard lib -----------------------------------------------
#include <string.h>          // string library

// local alg lib -----------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ==============================================================




// =========================================================
/// This function generates the initial conditions for the NS system:
void MGSolNSC::ic_read(
  int bc_gam,
  int mat_gam, 
  double xp[],
  double u_value[]
) {// =======================================
//xp[] is the NON-DIMENSIONAL node coordinate
  double ILref = 1./_lref;
  
  u_value[0]=0.; u_value[1]=0.;
  
  
  return;
}

// ========================================
// This function  defines the boundary conditions for the NS system:
void MGSolNSC::bc_read(
  int bc_gam,
  int mat_gam, 
  double xp[],          // xp[] is the NON-DIMENSIONAL node coordinates
  int bc_Neum[],        // normal
  int bc_flag[]         // boundary condition flag
){


  const double Lref = _mgphys.get_par("Lref");
  double ILref = 1./Lref;

  bc_flag[0]=0; bc_flag[1]=0;
  bc_Neum[0]=1; bc_Neum[1]=1; 


if (xp[2]>13.485-BDRY_TOLL){  // top of the tube
  
    bc_flag[1]=0; bc_flag[0]=0;
    bc_Neum[1]=0; bc_Neum[0]=1;
    
}

if (xp[2]<0.+BDRY_TOLL){  // bottom of the tube
  
    bc_flag[1]=1; bc_flag[0]=0;
    bc_Neum[1]=0; bc_Neum[0]=1;
    
}

  return;
}


#endif  //ENDIF NS_EQUATIONS

