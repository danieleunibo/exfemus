// libc+++ include
#include <iostream>
#include <cstdlib>
#include <sstream>

// configuration files -------------------------
#include   "Printinfo_conf.h"


// solver library -------------------------------------
#include  "Solverlib_conf.h"  // Solver library options 
// Petsc
#ifdef HAVE_PETSCM
#include "petsc.h" // for Petsc solver
#endif
// Mpi
#ifdef HAVE_MPI
#include <mpi.h>   //For MPI_COMM_WORLD
#endif

// include local class
#include "MGFemusInit.h"
#include "MGUtils.h"
#include "MGSystem.h"
#include "MGGeomEl.h"
#include "MGMesh.h"
#include "MGFEMap.h"
#include "MGFE.h"
#include "MGEquationsSystem.h"
#include "MGTimeLoop.h"
#include "FEMUS.h"

#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#endif
/// Set up
// =======================================
// Main program
// =======================================

int main(int argc, char** argv) {

  argc = argc ; argv=argv;  // no unused warning
 
  std::cout<<" ============ MGUtils ===================================== \n";
  std::cout<<" =========================================================== \n";
  // setting MGUtils -> parameters and file name ------------------------
  std::vector<MGUtils*> mgutils;
  std::string mesh_nameP[NUM_MESH];
  std::ostringstream filenameP[2];  std::ostringstream osfilenameP[2];

  for(int i_mesh=0; i_mesh< NUM_MESH; i_mesh++) {
    // MGUtils constructor ----------------------------------------------------
    mgutils.push_back(new MGUtils(i_mesh+1));
    // mesh names -------------------------------------------------------------
    mesh_nameP[i_mesh]= mgutils[i_mesh]->get_file("F_MESH_READ"); // name mesh
    int posP = mesh_nameP[i_mesh].find(".");  // position of "live" in str
    filenameP[i_mesh] <<   mesh_nameP[i_mesh].substr(0,posP)  << "_gen.med" ;
    osfilenameP[i_mesh]<< mgutils[i_mesh]->_mesh_dir <<filenameP[i_mesh].str();
    std::cout<<" \n P mesh file "<< i_mesh+1 << "= "
                                  << osfilenameP[i_mesh].str().c_str() <<"\n "; 
  }
  std::cout<<" ============ end loop mesh ================================ \n";
  std::cout<<" =========================================================== \n";
  
// FEM class -----------------------------------------------------------
  MGFEMap *mgfemap; mgfemap=new MGFEMap();
  // MGFEMap mg_femap;
  MGFE *dfe_q;    dfe_q=new MGFE(2,ELTYPE); dfe_q->init_qua();
  mgfemap->set_FE(dfe_q); //// initialize quadratic fem
  MGFE *dfe_l;  dfe_l=new MGFE(1,ELTYPE); dfe_l->init_lin();
  mgfemap->set_FE(dfe_l); //initialize linear fem
  MGFE *dfe_k; dfe_k=new MGFE(0,ELTYPE);  dfe_k->init_pie();
  mgfemap->set_FE(dfe_k); //  initialize piecewise fem

  // MGGeomEl ----------------------------------------------------------
  MGGeomEl *mggeomel;  mggeomel=new  MGGeomEl();
  
  // system problem =========================================================
  std::vector<FIELDS> myproblemP; myproblemP.resize(4);
  // Problem to solve for each mesh
  myproblemP[0]=NS_F; 
  myproblemP[1]=FS_F;  // calling NS adjoint equation
  myproblemP[2]=K_F;
  myproblemP[3]=KTT_F; // calling K adjoint equation
  
  // system 1
  // MGFemusInit --------------------------------------------------------------
  FEMUS P;                                        // constructor
  P.init_param(*mgutils[0]);                      // init parameter
  P.init_fem(*mggeomel,*mgfemap);                 // init fem      
  // setting mesh -------------------------------------------------------------
  P.setMesh();                                    // set MGmesh   
  // setting system -----------------------------------------------------------
  P.setSystem(myproblemP);                         // set system

  // solving
  int    n_steps = mgutils[0]->get_par("nsteps");
  double      dt = mgutils[0]->get_par("dt");
  int print_step = mgutils[0]->get_par("printstep");
  int    itime_0  = mgutils[0]->get_par("itime");
  double time    = 0.;
  
  P.solve_setup(itime_0,time);                    // initial time loop (t=0)

  for(int itime=itime_0; itime< itime_0 + n_steps; itime++) {

    P.solve_onestep(itime_0,itime,print_step,time,dt);    // solving P

  }   // end time loop

  // end ======================================================================
  // --------------------------------------------------------------------------
  P.terminate(); 

  // clean --------------------------------------------------------------------
  mgutils.clear();
  delete dfe_q;  delete dfe_l;   delete dfe_k;  // delete   fem
  delete mggeomel; delete mgfemap;
  return 0;
}