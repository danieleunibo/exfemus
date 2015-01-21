// libc++ include
#include <iostream>
#include <cstdlib>
#include <sstream>

// configuration files -------------------------
#include   "Printinfo_conf.h"

// LibMesh library included ------------------------------
// #ifdef LM_INIT
// #include "libmesh.h" // for Libmesh library
// #endif

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

#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#endif

#include "FEMUS.h"



// struct sBC {
//   std::string name;
//   std::string type;
//   bool isAnalytic;
//   std::string equation;
//   ParaMEDMEM::MEDCouplingFieldDouble * field;
// };
//
// void getInput(FEMUS &, std::vector<sBC> & boundary);


/// Set up
// =======================================
// Main program
// =======================================
// double pass_common[10];
int main(int argc, char** argv) {

  argc = argc ; argv=argv;  // no unused warning
  // setting MGUtils -> parameters and file name ------------------------
  std::vector<MGUtils*> mgutils;
  for(int mesh=0; mesh<2; mesh++) mgutils.push_back(new MGUtils(mesh+1));

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


  // system 1
  // MGFemusInit --------------------------------------------------------------

  FEMUS P;
  P.init_param(*mgutils[0]);
  P.init_fem(*mggeomel,*mgfemap);
  // setting mesh -------------------------------------------------------------
  std::string mesh_dir= mgutils[0]->_mesh_dir;     // name mesh
  P.setMedMesh(mesh_dir+"/test1quad9_group_mat_gen.med");// set med-mesh
  P.setMesh();
  // setting system -----------------------------------------------------------
  std::vector<NS_FIELDS> myproblem; myproblem.resize(1);
  myproblem[0]=DA_F;
  P.setSystem(myproblem);
  P.init_interfaces("/test1quad9_group_mat_gen.med");
  P.set_param_BCzone();

  // system 2
  // MGFemusInit --------------------------------------------------------------
  FEMUS P1;
  P1.init_param(*mgutils[1]);
  P1.init_fem(*mggeomel,*mgfemap);
  // setting mesh ----------------------------------------
  P1.setMedMesh(mesh_dir+"/test1quad9_group_mat_gen.med"); // set femus-mesh (P._mg_mesh[])
  P1.setMesh();
  // setting system ----------------------------------------
  P1.setSystem(myproblem);
  P1.init_interfaces("/test1quad9_group_mat_gen.med");
  P1.set_param_BCzone(P);


  // solving ==================================================================
  P.set_interface_userbc();
  P1.set_interface_userbc();
  P.solve();
  P1.solve();


  // end ======================================================================
  // --------------------------------------------------------------------------
  P.terminate();
  P1.terminate();

  // clean --------------------------------------------------------------------
  mgutils.clear();
  delete dfe_q;  delete dfe_l;   delete dfe_k;  // delete   fem
  delete mggeomel; delete mgfemap;
  return 0;
}


// // =========================================================
// // get input boundary condition
// void getInput(FEMUS &P, std::vector<sBC> & boundary) {
//
//   std::vector<std::string> bcNames = P.get_interface_Names();
//   int nBC = bcNames.size();
// //    ParaMEDMEM::MEDCouplingFieldDouble *bdy=NULL;
//   int face_id=P.get_interface_byname(1000);
//   ParaMEDMEM::MEDCouplingFieldDouble* bdy=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES,ParaMEDMEM::NO_TIME);
//   bdy->setName("MyScalarFieldOnNode");
//   bdy->setMesh(P.get_interface_mesh(face_id));
//   //ParaMEDMEM::DataArrayDouble *array=ParaMEDMEM::DataArrayDouble::New();
//   //array->alloc(bdy->getMesh()->getNumberOfNodes(),1);//Implicitely fieldOnNodes will be a 1 component field.
//   //array->fillWithValue(17.);
//   //bdy->setArray(array);
//   //array->decrRef();
//   bdy->fillFromAnalytic(1,"17.*x*x");
//
//   std::vector<const ParaMEDMEM::MEDCouplingFieldDouble *> v_fs(1,bdy);
//   ParaMEDMEM::MEDCouplingFieldDouble::WriteVTK("bdy.vtk",v_fs);
//
// //    bdy = P.getValuesOnBoundary("1000", "u");
//
// //   // Field
// // //   std::string inFieldFile = "MESH/in_square2_2.med";
// // //    std::string inFieldFile = "RESU/test2_bdy_square2_2_pass.med";
// //
// // //   std::vector<std::string> InFieldNames = MEDLoader::GetAllFieldNames(inFieldFile.c_str());
// // //   std::vector<std::string> InMeshNames  = MEDLoader::GetMeshNames(inFieldFile.c_str());
// //
// // //   MEDCouplingFieldDouble *infield =
// // //      MEDLoader::ReadField(ON_NODES, inFieldFile.c_str(), InMeshNames[0].c_str(),
// // //                           0, InFieldNames[0].c_str(), -1, -1);
// //
// //   // boundary condition construction
// //   std::vector<std::string> bcNames = P.getBoundaryNames();
// //   int nBC = bcNames.size();
// //   boundary.resize(nBC);
// //
// //   for (int i=0; i<nBC; i++) {
// //     boundary[i].name = bcNames[i];
// //     boundary[i].type = "Neumann";
// //     boundary[i].isAnalytic = true;
// //
// //
// //   if( bcNames[i]=="inlet"){
// // //   boundary[i].type = "Dirichlet";
// // //   boundary[i].isAnalytic = true;
// // //   boundary[i].field = NULL;
// // //   boundary[i].equation = "1";
// //   }
// //
// //   if( bcNames[i]=="lateral"){
// //   boundary[i].type = "Dirichlet";
// //   boundary[i].isAnalytic = true;
// // //   boundary[i].field = infield ;
// // //   boundary[i].field = NULL;
// //   boundary[i].equation = "0";
// //  }
// //   if( bcNames[i]=="share"){
// //  //   boundary[i].type = "Dirichlet";
// // //   boundary[i].isAnalytic = true;
// // //   boundary[i].equation = "2.0";
// // //   boundary[i].field = NULL;
// //     }
// //   }
// //
// //
//
//   /// 3D  ////////////
//
//   // Field
// //  std::string inFieldFile = "MESH/in_square2_2.med";
// //    std::string inFieldFile = "RESU/test2_bdy_square2_2_pass.med";
//
// //   std::vector<std::string> InFieldNames = MEDLoader::GetAllFieldNames(inFieldFile.c_str());
// //   std::vector<std::string> InMeshNames  = MEDLoader::GetMeshNames(inFieldFile.c_str());
//
// //   MEDCouplingFieldDouble *infield =
// //      MEDLoader::ReadField(ON_NODES, inFieldFile.c_str(), InMeshNames[0].c_str(),
// //                           0, InFieldNames[0].c_str(), -1, -1);
//
//   // boundary condition construction:
//   // *****************************************
// //   struct sBC (boundary[i]) has to define:
// //        std::string name;
// //        std::string type;
// //        bool isAnalytic;
// //        std::string equation;
// //        ParaMEDMEM::MEDCouplingFieldDouble * field;
//
//
//
//   boundary.resize(nBC);
//
//   for (int i=0; i<nBC; i++) {
//     boundary[i].name = bcNames[i];// sBC struct name
//     boundary[i].type = "UnMark"; // sBC struct type
//     boundary[i].isAnalytic = true;// sBC struct flag
//     boundary[i].equation = "";    // sBC struct expression
//     boundary[i].field = NULL;     // sBC struct num field
//
//     if (bcNames[i]=="1000") {
//       boundary[i].type = "Mark";
//       boundary[i].isAnalytic = false;
// //       boundary[i].field = bdy;
//       boundary[i].field = NULL;
// //           boundary[i].equation = "IVec*1.0";
//     }
//     if (bcNames[i]=="36") {
//       boundary[i].type = "Mark";
//       boundary[i].isAnalytic = true;
// //       boundary[i].field = bdy;
//       boundary[i].field = NULL;
//       boundary[i].equation = "IVec*0";
//     }
// //     if (bcNames[i]=="bottom") {
// //       boundary[i].type = "Dirichlet";
// //       boundary[i].isAnalytic = true;
// //       boundary[i].field = NULL;
// //       boundary[i].equation = "IVec*0.0 + JVec*0.0+ KVec*1.0";
// //     }
// //     if(bcNames[i]=="Group_3") {
// //       boundary[i].type = "Dirichlet";
// //       boundary[i].isAnalytic = true;
// //       boundary[i].field = NULL;
// //       boundary[i].equation = "IVec*0.0 + JVec*0.0+ KVec*0.0";
// //     }
//   }
//
//
//
//
//   return;
// }
// =========================================================
