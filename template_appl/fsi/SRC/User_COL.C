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
  
   if (bc_gam<20){ bc_flag[0]=1;bc_Neum[0]=1; }
   else if (bc_gam<30){ bc_flag[0]=1;bc_Neum[0]=3 ;}
   if (bc_gam==1000)   { bc_flag[0]=1;bc_Neum[0]=5;}
   
//         if (xp[0]<1-1.e-6){ bc_flag[0]=1;bc_Neum[0]=1; }
//    if (xp[0]>1+1.e-6){ bc_flag[0]=1;bc_Neum[0]=3 ;}
//    if (fabs(xp[0]-1.)<1.e-6)   { bc_flag[0]=1;bc_Neum[0]=5;}
  return;
 }
#endif
