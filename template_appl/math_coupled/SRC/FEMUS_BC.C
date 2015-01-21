#include <iostream>
#include <cstdlib>
#include <sstream>

#define FDEBUG (0)

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
#include "MGFEMap.h"
#include "MGFE.h"
#include "MGEquationsSystem.h"
#include "MGTimeLoop.h"


// class include
// #include "conf.hxx"
#include "FEMUS.h"

// #ifdef HAVE_MED
// MED includes
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "InterpKernelExprParser.hxx"
// #endif

// LibMeshCpp includes
#include "ParallelMeshExtended.h"
#include "EquationSystemsExtended.h"
// #include "_LibMeshProblem.h"
// #include "Debug.hxx"


// ============================================================================
// These functions initialize interface properties for coupling
// Note:
// _boundary[i].name      = UsersBC struct name 
//         (automatic)
// _boundary[i].type      = UsersBC struct type: 
//          example: Marked or UnMarked    -> (user in init_BCzone())
// _boundary[i].isAnalytic= UsersBC struct flag
//          example: true or false (user in init_BCzone())
// _boundary[i].equation  = UsersBC struct expression
//         example: "IVec*(z+2) + JVec*(y+x)+ KVec*(x+1.0)"
// _boundary[i].field     = UsersBC struct num field  (NULL or pointer)
//         example: NULL or ParaMEDMEM::MEDCouplingFieldDouble bdybdy pointer 
//=============================================================================
// ============================================================================
// This function initializes the interface properties for coupling
void  FEMUS::init_BCzone(std::vector<std::string> &bcNames
) {// =========================================================================
 
//   std::vector<std::string> bcNames = get_interface_Names();
  int nBC = bcNames.size();
  _boundary.resize(nBC);

  for (int i=0; i<nBC; i++) {
    _boundary[i].name = bcNames[i];  // UsersBC struct name
    _boundary[i].type = "UnMark";    // UsersBC struct type
    _boundary[i].from_cmp = 0;       // UsersBC initial cmp
    _boundary[i].n_cmp = 1;          // UsersBC n of cmp
    _boundary[i].order_cmp = 2;      // UsersBC order cmp
    _boundary[i].isAnalytic = false; // UsersBC struct flag
    _boundary[i].equation = "";      // UsersBC struct expression
    _boundary[i].field = NULL;       // UsersBC struct num field
    _boundary[i].on_nodes = true;    // UsersBC node field (true)

    // ============  USER AREA ================================================
    // please mark the interface used in coupling
    // for mat only elements (_boundary[i].on_nodes = false) (automatic)
    if (bcNames[i]=="1"){
     _boundary[i].type = "Mark";
     _boundary[i].on_nodes = false; // (default)
    }
    // for bc only nodes (_boundary[i].on_nodes = true) (at the moment)
    if (bcNames[i]=="1000"){
      _boundary[i].type = "Mark";
      _boundary[i].on_nodes = true; // (default)
    }
    if (bcNames[i]=="36") {
      _boundary[i].type = "Mark";
      _boundary[i].on_nodes = true;// (default)
    } 
    
    // ============  USER AREA ================================================
    
  }

  return;
}


// ============================================================================
// This function sets the interface properties for coupling
// =========================================================
void  FEMUS::set_param_BCzone(){

  std::vector<std::string> bcNames = get_interface_Names();
  int nBC = bcNames.size();
//    ParaMEDMEM::MEDCouplingFieldDouble *bdy=NULL;
//   int face_id=search_idxBC_byid(1000);
//   ParaMEDMEM::MEDCouplingFieldDouble* bdy=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES,ParaMEDMEM::NO_TIME);
//   bdy->setName("MyScalarFieldOnNode");
//   bdy->setMesh(get_interface_mesh(face_id));
//   //ParaMEDMEM::DataArrayDouble *array=ParaMEDMEM::DataArrayDouble::New();
//   //array->alloc(bdy->getMesh()->getNumberOfNodes(),1);//Implicitely fieldOnNodes will be a 1 component field.
//   //array->fillWithValue(17.);
//   //bdy->setArray(array);
//   //array->decrRef();
//   bdy->fillFromAnalytic(2,"JVec*x+IVec*1");

//   std::vector<const ParaMEDMEM::MEDCouplingFieldDouble *> v_fs(1,bdy);
//   ParaMEDMEM::MEDCouplingFieldDouble::WriteVTK("bdy.vtk",v_fs);

//    bdy = P.getValuesOnBoundary("1000", "u");

//   // Field
// //   std::string inFieldFile = "MESH/in_square2_2.med";
// //    std::string inFieldFile = "RESU/test2_bdy_square2_2_pass.med";
//
// //   std::vector<std::string> InFieldNames = MEDLoader::GetAllFieldNames(inFieldFile.c_str());
// //   std::vector<std::string> InMeshNames  = MEDLoader::GetMeshNames(inFieldFile.c_str());
//
// //   MEDCouplingFieldDouble *infield =
// //      MEDLoader::ReadField(ON_NODES, inFieldFile.c_str(), InMeshNames[0].c_str(),
// //                           0, InFieldNames[0].c_str(), -1, -1);
//
//   // boundary condition construction
//   std::vector<std::string> bcNames = P.getBoundaryNames();
//   int nBC = bcNames.size();
//   boundary.resize(nBC);
//
//   for (int i=0; i<nBC; i++) {
//     _boundary[i].name = bcNames[i];
//     _boundary[i].type = "Neumann";
//     _boundary[i].isAnalytic = true;
//
//
//   if( bcNames[i]=="inlet"){
// //   _boundary[i].type = "Dirichlet";
// //   _boundary[i].isAnalytic = true;
// //   _boundary[i].field = NULL;
// //   _boundary[i].equation = "1";
//   }
//
//   if( bcNames[i]=="lateral"){
//   _boundary[i].type = "Dirichlet";
//   _boundary[i].isAnalytic = true;
// //   _boundary[i].field = infield ;
// //   _boundary[i].field = NULL;
//   _boundary[i].equation = "0";
//  }
//   if( bcNames[i]=="share"){
//  //   _boundary[i].type = "Dirichlet";
// //   _boundary[i].isAnalytic = true;
// //   _boundary[i].equation = "2.0";
// //   _boundary[i].field = NULL;
//     }
//   }
//
//

  /// 3D  ////////////

  // Field
//  std::string inFieldFile = "MESH/in_square2_2.med";
//    std::string inFieldFile = "RESU/test2_bdy_square2_2_pass.med";

//   std::vector<std::string> InFieldNames = MEDLoader::GetAllFieldNames(inFieldFile.c_str());
//   std::vector<std::string> InMeshNames  = MEDLoader::GetMeshNames(inFieldFile.c_str());

//   MEDCouplingFieldDouble *infield =
//      MEDLoader::ReadField(ON_NODES, inFieldFile.c_str(), InMeshNames[0].c_str(),
//                           0, InFieldNames[0].c_str(), -1, -1);

  // boundary condition construction:
  // *****************************************
//   struct sBC (_boundary[i]) has to define:
//        std::string name;
//        std::string type;
//        bool isAnalytic;
//        std::string equation;
//        ParaMEDMEM::MEDCouplingFieldDouble * field;



  _boundary.resize(nBC);

  for (int i=0; i<nBC; i++) if(_boundary[i].type =="Mark"){
    
    // default 
    _boundary[i].from_cmp = 0;       // UsersBC initial cmp
    _boundary[i].n_cmp = 1;          // UsersBC n of cmp
    _boundary[i].order_cmp = 2;      // UsersBC order cmp
    _boundary[i].isAnalytic = false; // sBC struct flag
    _boundary[i].equation = "";   // sBC struct expression
    _boundary[i].n_cmp = 1;
    _boundary[i].field = NULL;      // sBC struct num field
    _boundary[i].on_nodes = true;   // UsersBC node field (true)
 
    // set to change the default
    if (bcNames[i]=="1000") {
      _boundary[i].isAnalytic = true;
//       _boundary[i].field = bdy;
      _boundary[i].n_cmp = 1;
     _boundary[i].field = NULL;
          _boundary[i].equation = "IVec*2";
    }
//       if (bcNames[i]=="1") {
//       _boundary[i].isAnalytic = true;
//       _boundary[i].field = NULL;
// //       _boundary[i].n_cmp = 1;
// //    _boundary[i].field = NULL;
//       _boundary[i].equation = "50";
//       _boundary[i].on_nodes = false;   // UsersBC node field (true)
//     }
    
    
    if (bcNames[i]=="36") {
      _boundary[i].isAnalytic = true;
//       _boundary[i].field = bdy;
      _boundary[i].field = NULL;
//        _boundary[i].n_cmp = 1;
      _boundary[i].equation = "IVec*12";
    }
//     if (bcNames[i]=="bottom") {
//       _boundary[i].type = "Dirichlet";
//       _boundary[i].isAnalytic = true;
//       _boundary[i].field = NULL;
//       _boundary[i].equation = "IVec*0.0 + JVec*0.0+ KVec*1.0";
//     }
//     if(bcNames[i]=="Group_3") {
//       _boundary[i].type = "Dirichlet";
//       _boundary[i].isAnalytic = true;
//       _boundary[i].field = NULL;
//       _boundary[i].equation = "IVec*0.0 + JVec*0.0+ KVec*0.0";
//     }
  }




  return;
}

// ============================================================================
// This function sets the interface properties for coupling
// =========================================================
void  FEMUS::set_param_BCzone(FEMUS & P_OLD){
  std::vector<std::string> bcNames = get_interface_Names();
  int nBC = bcNames.size();
  
  ParaMEDMEM::MEDCouplingFieldDouble *bdy=NULL;  
  bdy = P_OLD.getValuesOnBoundary("1000", "DA",0);
  
//   ParaMEDMEM::MEDCouplingFieldDouble *bdy=NULL;
//   int face_id=search_idxBC_byid(1000);
//   ParaMEDMEM::MEDCouplingFieldDouble* bdy=ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES,ParaMEDMEM::NO_TIME);
//   bdy->setName("MyScalarFieldOnNode");
//   bdy->setMesh(get_interface_mesh(face_id));
//   //ParaMEDMEM::DataArrayDouble *array=ParaMEDMEM::DataArrayDouble::New();
//   //array->alloc(bdy->getMesh()->getNumberOfNodes(),1);//Implicitely fieldOnNodes will be a 1 component field.
//   //array->fillWithValue(17.);
//   //bdy->setArray(array);
//   //array->decrRef();
//   bdy->fillFromAnalytic(2,"JVec*x+IVec*1");

//   std::vector<const ParaMEDMEM::MEDCouplingFieldDouble *> v_fs(1,bdy);
//   ParaMEDMEM::MEDCouplingFieldDouble::WriteVTK("bdy.vtk",v_fs);

//    bdy = P.getValuesOnBoundary("1000", "u");

//   // Field
// //   std::string inFieldFile = "MESH/in_square2_2.med";
// //    std::string inFieldFile = "RESU/test2_bdy_square2_2_pass.med";
//
// //   std::vector<std::string> InFieldNames = MEDLoader::GetAllFieldNames(inFieldFile.c_str());
// //   std::vector<std::string> InMeshNames  = MEDLoader::GetMeshNames(inFieldFile.c_str());
//
// //   MEDCouplingFieldDouble *infield =
// //      MEDLoader::ReadField(ON_NODES, inFieldFile.c_str(), InMeshNames[0].c_str(),
// //                           0, InFieldNames[0].c_str(), -1, -1);
//
//   // boundary condition construction
//   std::vector<std::string> bcNames = P.getBoundaryNames();
//   int nBC = bcNames.size();
//   boundary.resize(nBC);
//
//   for (int i=0; i<nBC; i++) {
//     _boundary[i].name = bcNames[i];
//     _boundary[i].type = "Neumann";
//     _boundary[i].isAnalytic = true;
//
//
//   if( bcNames[i]=="inlet"){
// //   _boundary[i].type = "Dirichlet";
// //   _boundary[i].isAnalytic = true;
// //   _boundary[i].field = NULL;
// //   _boundary[i].equation = "1";
//   }
//
//   if( bcNames[i]=="lateral"){
//   _boundary[i].type = "Dirichlet";
//   _boundary[i].isAnalytic = true;
// //   _boundary[i].field = infield ;
// //   _boundary[i].field = NULL;
//   _boundary[i].equation = "0";
//  }
//   if( bcNames[i]=="share"){
//  //   _boundary[i].type = "Dirichlet";
// //   _boundary[i].isAnalytic = true;
// //   _boundary[i].equation = "2.0";
// //   _boundary[i].field = NULL;
//     }
//   }
//
//

  /// 3D  ////////////

  // Field
//  std::string inFieldFile = "MESH/in_square2_2.med";
//    std::string inFieldFile = "RESU/test2_bdy_square2_2_pass.med";

//   std::vector<std::string> InFieldNames = MEDLoader::GetAllFieldNames(inFieldFile.c_str());
//   std::vector<std::string> InMeshNames  = MEDLoader::GetMeshNames(inFieldFile.c_str());

//   MEDCouplingFieldDouble *infield =
//      MEDLoader::ReadField(ON_NODES, inFieldFile.c_str(), InMeshNames[0].c_str(),
//                           0, InFieldNames[0].c_str(), -1, -1);

  // boundary condition construction:
  // *****************************************
//   struct sBC (_boundary[i]) has to define:
//        std::string name;
//        std::string type;
//        bool isAnalytic;
//        std::string equation;
//        ParaMEDMEM::MEDCouplingFieldDouble * field;



  _boundary.resize(nBC);

  for (int i=0; i<nBC; i++) if(_boundary[i].type =="Mark"){
    
        // default 
    _boundary[i].from_cmp = 0;       // UsersBC initial cmp
    _boundary[i].n_cmp = 1;          // UsersBC n of cmp
    _boundary[i].order_cmp = 2;      // UsersBC order cmp
    _boundary[i].isAnalytic = false; // sBC struct flag
    _boundary[i].equation = "";   // sBC struct expression
    _boundary[i].n_cmp = 1;
    _boundary[i].field = NULL;      // sBC struct num field
    _boundary[i].on_nodes = true;   // UsersBC node field (true)
    

    if (bcNames[i]=="1000") {
      _boundary[i].isAnalytic = false;
//       _boundary[i].field = bdy;
//       _boundary[i].n_cmp = 1;
     _boundary[i].field = bdy;
//           _boundary[i].equation = "IVec*2";
     _boundary[i].equation = "";
    }
//       if (bcNames[i]=="1") {
//       _boundary[i].isAnalytic = true;
//       _boundary[i].field = NULL;
//       _boundary[i].n_cmp = 1;
// //      _boundary[i].field = NULL;
//          _boundary[i].equation = "50";
//     }
    
    
    if (bcNames[i]=="36") {
      _boundary[i].isAnalytic = true;
//       _boundary[i].field = bdy;
      _boundary[i].field = NULL;
//        _boundary[i].n_cmp = 1;
      _boundary[i].equation = "IVec*12";
    }
//     if (bcNames[i]=="bottom") {
//       _boundary[i].type = "Dirichlet";
//       _boundary[i].isAnalytic = true;
//       _boundary[i].field = NULL;
//       _boundary[i].equation = "IVec*0.0 + JVec*0.0+ KVec*1.0";
//     }
//     if(bcNames[i]=="Group_3") {
//       _boundary[i].type = "Dirichlet";
//       _boundary[i].isAnalytic = true;
//       _boundary[i].field = NULL;
//       _boundary[i].equation = "IVec*0.0 + JVec*0.0+ KVec*0.0";
//     }
  }




  return;
}

// =========================================================
