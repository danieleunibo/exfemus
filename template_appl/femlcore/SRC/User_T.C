 
#include "ReactData.h"
#include "Equations_conf.h"

// ============================================
#ifdef T_EQUATIONS // 3D-2D Energy equation
// ============================================
#include <sstream>

// // class local configuration -------
 #include "MGSolverT.h"
// 
// // local include -------------------
 #include "MGSystem.h"



/// This function generates the initial conditions for the energy equation:
// =================================================
void MGSolT::ic_read(int bc_gam, int mat_gam, double xp[],double u_value[]) {
// =================================================
  
    u_value[0] = 673.15;
    
}

////bc gambit
// ========================================================
/// This function  defines the boundary conditions for the energy equation:
void MGSolT::bc_read(int bc_gam, int mat_gam, double xp[],int bc_Neum[],int bc_flag[]) {
// =================================================
  
  
  std::string msh_name=_mgutils.get_file("MESHNAME");
  
  if (msh_name=="elsy") { 
    bc_Neum[0]=1;   bc_flag[0]=0;
 /*
    if ((bc_gam==3100)||(bc_gam==3200)||(bc_gam==3300)||(bc_gam==3400)||(bc_gam==3500)){
       bc_Neum[0]=0; bc_flag[0]=2;}*/

 
    if (bc_gam==4000){
       bc_Neum[0]=0; bc_flag[0]=2;}
 
 
  }
  
  if (msh_name=="box") {
    
    const double Lref = 1.; //_mgphys.get_par("Lref");
    double ILref = 1./Lref;
    
    if (xp[2]< LZB*ILref+BDRY_TOLL )
    {bc_Neum[0]=0; bc_flag[0]=2;}

  }
  
  if (msh_name=="elsy_core") { 
    
    if (xp[2]< BDRY_TOLL )
    {bc_Neum[0]=0; bc_flag[0]=2;}

  }
  

  return;
  
} // end boundary conditions ==========================



/// This function generates the initial conditions for the energy equation:
// =================================================
void MGSolTC::ic_read(int bc_gam, int mat_gam, double xp[],double u_value[]) {
// =================================================
  
    u_value[0] = 608.15;
    
    }

////bc gambit
// ========================================================
/// This function  defines the boundary conditions for the energy equation:
void MGSolTC::bc_read(int bc_gam, int mat_gam, double xp[],int bc_Neum[],int bc_flag[]) {
// ========================================================

    const double Lref = 1.; //_mgphys.get_par("Lref");
    double ILref = 1./Lref;

   bc_Neum[0]=1; bc_flag[0]=0;
  
   if (xp[2]> 13.485-BDRY_TOLL )  { // top of the tube
     bc_flag[0]=2;    
     bc_Neum[0]=0;
   }
   
//    if (xp[2]< 12.785)  { // fixed temperature under the gv
//      bc_flag[0]=2;    
//      bc_Neum[0]=0;
//    }

  return;
  
} // end boundary conditions ==========================


#endif
