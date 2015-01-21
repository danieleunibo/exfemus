#include <iostream>
#include <cstdlib>
#include <sstream>

// configuration files -------------------------
#include   "Printinfo_conf.h"

// LibMesh library included ------------------------------
#ifdef LM_INIT
#include "libmesh.h" // for Libmesh library 
#endif

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
#include "MGEquationsMap.h"
#include "MGTimeLoop.h"


#include "ReactData.h"

/// Set up
// =======================================
// Main program
// =======================================
// double pass_common[10];
int main(int argc, char** argv) {

  // ================ MGFemusInit ================================
#ifdef LM_INIT
  LibMeshInit init(argc, argv);
//   std::cout << libMesh::processor_id() << std::endl;
#else
#ifndef HAVE_LASPACKM
  MGFemusInit start(argc,argv);
#endif
#endif

  // ============== MGUtils =======================================
  std::vector<MGUtils*> mg_utils;
  for(int mesh=0; mesh<NUM_MESH; mesh++) mg_utils.push_back(new MGUtils(mesh+1)); // MGUtils class: read and print parameters 
  
  // ============== MGSystem ======================================
  int np_data=0;    int ncell_data=2;
#ifdef TBK_EQUATIONS
  np_data =1;
#endif
  std::vector<MGSystem*> mg_phys;
  for(int mesh=0; mesh<NUM_MESH; mesh++) {
    mg_phys.push_back(new MGSystem(*mg_utils[mesh],NULL,np_data,ncell_data));
    mg_phys[mesh]->read_par();
#ifdef PRINT_INFO  // ---- info ---------------
    mg_phys[mesh]->print_par();       // print parameters
#endif   
  }
  
// ======  MGMesh ================================================
  std::vector<MGMesh*> mg_mesh;
  MGGeomEl mg_geomel;  // set element type
  
  for(int mesh=0; mesh<NUM_MESH; mesh++) {
  const double Lref  = mg_phys[mesh]->get_par("Lref");      // reference L
  const int NoLevels= mg_utils[mesh]->get_par("nolevels");  // numb of Level
  mg_mesh.push_back(new MGMesh(start.comm(), *mg_utils[mesh],mg_geomel,Lref));
    if (NoLevels != mg_mesh[mesh]->_NoLevels) {
    std::cout << "Inconsistent Number of Levels between Mesh and SolBase"
              << std::endl;
    abort();
    }
    mg_mesh[mesh]->print(NoLevels-1,0);
    mg_phys[mesh]->set_mesh(mg_mesh[mesh]);
    mg_phys[mesh]->init_data(0);
  }

// ===== MGFEMap =========================================
  MGFEMap mg_femap;   // fem map

  MGFE dfe_q(2,ELTYPE); dfe_q.init_qua();      // initialize quadratic fem
  mg_femap.set_FE(&dfe_q); // 
  MGFE dfe_l(1,ELTYPE); dfe_l.init_lin();      // initialize linear fem
  mg_femap.set_FE(&dfe_l); // 
  MGFE dfe_k(0,ELTYPE); dfe_k.init_pie();      // initialize piecewise fem
  mg_femap.set_FE(&dfe_k); // 
  
  
// ===== MGEEquationsMap =========================================
  std::vector<MGEquationsMap*> mg_equations_map;
  std::vector<ReactData*> react_data;
//   MGEquationsMap *mg_equations_map[NUM_MESH];
  for(int mesh=0; mesh<NUM_MESH; mesh++) {
  mg_equations_map.push_back(new MGEquationsMap(*mg_utils[mesh], *mg_phys[mesh] ,*mg_mesh[mesh],mg_femap));  // MGEquationsMap class
  mg_equations_map[mesh]->init();                                    // adds the equations to the map
  mg_equations_map[mesh]->setDofBcOpIc();                           // set operators
  react_data.push_back(new ReactData());
  react_data[mesh]->datagen(*mg_mesh[mesh],*mg_phys[mesh],*mg_utils[mesh]);
  }
  
//   // ======================  ReactData =====================================================
//   ReactData react_data;
//   react_data.datagen(mgmesh,mg_phys,mg_utils);                         // read fuel and pressure drop
//   react_data.datagen(*mg_mesh[0],mg_phys1,mg_utils1);      
//    if (num_mesh>1) react_data.datagen(*mg_mesh[1],mg_phys2,mg_utils2);
   
   
  // =================== MGTimeLoop ========================================================
  MGTimeLoop mg_time_loop(mg_utils,mg_equations_map);         // MGTimeLoop class
  int      t_in=0;  double    time=0.;                        //  initial time
  const int restart      = mg_utils[0]->get_par("restart");   // restart or not

  mg_time_loop.transient_setup(restart,t_in,time);             //  MGTimeLoop: setup
  mg_time_loop.transient_loop(t_in,time);                      //  MGTimeLoop: Solver & Output

  
  // ==================  clean =========================================================
#ifdef HAVE_PETSCM
  std::string petsc_femus_log = mg_utils[0]->get_file("PETSC_FEMUS_LOG");
  std::ostringstream petsc_log;
  petsc_log <<"./" << petsc_femus_log;
//   PetscLogPrintSummary(MPI_COMM_WORLD,petsc_log.str().c_str());
#endif
  
  return 0;
}