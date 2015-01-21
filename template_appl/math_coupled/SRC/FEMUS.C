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

// Libmesh test case includes
// #include "LaplaceProblem.hxx"
// #include "NS_Problem.hxx"
// #include "NS_Problem3D.hxx"
// #include "LinearElasticProblem.hxx"
// #include "Con_DiffProblem.hxx"
// #include "MonodimProblem.hxx"

// ****************************************************************************
// ****************  Constructor Destructor ***********************************

// ============================================================================
// Basic constructor
FEMUS::FEMUS()  :
  _comm(MPI_COMM_WORLD)    // communicator
//   _mg_utils(NULL),          // parameters and files
//   _mg_geomel(NULL),         // fem geom class element
//   _mg_femap(NULL)        // fem
//    _num_mgmesh(NUM_MESH_MAIN)           // n of femus_meshes
//   _mg_mesh(NULL),           // MG meshes
//   _med_mesh(NULL)         // top med mesh
  
  
//   _mg_equations_map(NULL),  // system equation
//   _mg_time_loop(NULL)
  {     // transient system
    
    
 
    
  //=========================================================================

  // Init MPI flag --------------------------------
  int flag=0;  MPI_Initialized(&flag);
  if (flag) {_local_MPI_Init = false; }
  else {_local_MPI_Init = true; }
  // femus init -----------------------------------
  int argc = 0;    char ** argv = NULL;
  _start=new  MGFemusInit(argc,argv,_comm);

  return;
}

// ============================================================================
// This function is a constructor with  communicator
FEMUS::FEMUS( 
  MPI_Comm comm
):
  _comm(comm)    // communicator
//   _mg_utils(NULL),          // parameters and files
//   _mg_geomel(NULL),         // fem geom class element
//   _mg_femap(NULL),          // fem
//   _num_mgmesh(NUM_MESH_MAIN)          // n of femus_meshes
//   _mg_mesh(NULL),           // MG meshes
//   _med_mesh(NULL)          // top med mesh
//   _mg_equations_map(NULL),  // system equation
//   _mg_time_loop(NULL) 
  {     // transient system



  // Init MPI flag
  int flag=0;  MPI_Initialized(&flag);
  if (flag) {_local_MPI_Init = false; }
  else {_local_MPI_Init = true; }

  // femus init
  int argc = 0;    char ** argv = NULL;
  _start=new  MGFemusInit(argc,argv);

  return;
}
// =======================================================================
void FEMUS::init_param(
  MGUtils &   mgutils
) { // ====================================================================
  _mg_utils=&mgutils;
  return;
}
// =======================================================================
void FEMUS::init_fem(
  MGGeomEl & mggeomel,
  MGFEMap & mgfemap
) { // ====================================================================
  // setting MGUtils
//  if(mgfemap[0] == NULL) {
//      std::cout<< "FEMUS::init: no _mg_utils or _mg_femap";
//      abort();
//    }
//   // setting MGFEMap (fem)
  _mg_geomel=&mggeomel;
  _mg_femap=&mgfemap;
  return;
}
// ============================================================================
// This function is the destructor
FEMUS::~FEMUS() {
  // ==========================================================================
  clean();
//   if (_problem) delete _problem;     // problem
  if (_med_mesh) _med_mesh->decrRef();       // med-mesh
  _interface_mesh_vect.clear();
  _boundary.clear(); 
  delete _start;
}

void FEMUS::clean() {
  // ==========================================================================
//   if (_problem) delete _problem;     // problem

  for (int i=0; i<(int) _interface_mesh_vect.size(); i++)
    if (_interface_mesh_vect[i].support) _interface_mesh_vect[i].support->decrRef();

  for (int i=0; i<(int) _boundary.size(); i++)
    if (_boundary[i].field) _boundary[i].field->decrRef();
}
// ============================================================================
// This function is the problem destructor
void FEMUS::terminate(
) {// =========================================================================
//    if (_problem) _problem->terminate();
#ifdef HAVE_PETSCM
//   std::string petsc_femus_log = _mg_utils[0]->get_file("PETSC_FEMUS_LOG");
//   std::ostringstream petsc_log;
//   petsc_log <<"./" << petsc_femus_log;
//   PetscLogPrintSummary(MPI_COMM_WORLD,petsc_log.str().c_str());
#endif
}

// // // ****************************************************************************
// // // ****************    end Constructor Destructor *****************************
// //
// // // ****************************************************************************
// // // ****************    Set    *************************************************


// // // =============================================================================
// // // This function sets the type of problem
// // // Poisson
// // // Elasticity
// // // convdiff     monodim
// // // navierstokes navierstokes3D
void FEMUS::setSystem(
  const std::vector<NS_FIELDS> & pbName
) {// ==========================================================================

  // MgSystem
//   int np_data=0;  int ncell_data=2;
//   for(int mesh=0; mesh<_num_mesh; mesh++) {
// //     _mg_phys.push_back(new MGSystem(*_mg_utils[mesh],NULL,np_data,ncell_data));
//     _mg_phys.push_back(new EquationSystemsExtended(*_mg_utils[mesh],*_mg_mesh[mesh],*_mg_femap,1));
//     _mg_phys[mesh]->read_par();
//
// #ifdef PRINT_INFO  // ---- info ---------------
//     _mg_phys[mesh]->print_par();       // print parameters
// #endif
// //     _mg_phys[mesh]->set_mesh(*_mg_mesh[mesh]);
//     _mg_phys[mesh]->init_data(0);
//   }


  //   MGEquationsMap *mg_equations_map[NUM_MESH];
//   for (int imesh=0; imesh<_num_mgmesh; imesh++) {
    _mg_equations_map.push_back(new EquationSystemsExtended(
                                  *_mg_utils,*_mg_mesh[0],*_mg_femap,1,0,0));  // MGEquationsMap class
    _mg_equations_map[0]->read_par();
#ifdef PRINT_INFO  // ---- info ---------------
    _mg_equations_map[0]->print_par();       // print parameters
#endif
    _mg_equations_map[0]->init_data(0);
    _mg_equations_map[0]->init(pbName);                                    // adds the equations to the map
    _mg_equations_map[0]->setDofBcOpIc();                            // set operators
    _mg_equations_map[0]->set_mesh_femus(*_mg_mesh[0]);
    _mg_equations_map[0]->set_mesh_med(*_med_mesh);
//   }
  assert(_mg_equations_map[0]!=NULL);

//   // Reading boundary conditions ------------------------------
//   std::string mesh_dir=_mg_utils[0]->_mesh_dir;
//   std::string localFile=mesh_dir+"/test1quad9_group_mat.med";
//     std::vector<std::string> meshNames = MEDLoader::GetMeshNames(localFile.c_str());
//   if(meshNames.size() < 1) {
//     std::cout<<  " FEMUS::setMesh : no meshes in the file'";
//   }
//    std::string localMeshName = meshNames[0];
//   std::vector<std::string> GroupNames =
//     MEDLoader::GetMeshGroupsNames(localFile.c_str(), localMeshName.c_str());
//
//   int nGroups = GroupNames.size();
//   std::vector<std::string> vG(1);
//    _interface_mesh_vect.resize(nGroups);
//   for(int i=0; i<nGroups; i++) {
//     vG[0] = GroupNames[i];
//     _interface_mesh_vect[i].name = GroupNames[i];
//     _interface_mesh_vect[i].support = MEDLoader::ReadUMeshFromGroups(localFile.c_str(), meshNames[0].c_str(), -1, vG);
//     if(_interface_mesh_vect[i].support) {
//      _interface_mesh_vect[i].support->zipCoords();
//      _interface_mesh_vect[i].id = defineBoundary(_interface_mesh_vect[i].support);
//     } else
//       _interface_mesh_vect[i].id = -1;
//   }


  //time loop
  _mg_time_loop=new  MGTimeLoop(*_mg_utils,_mg_equations_map);
  assert(_mg_time_loop!=NULL);

// //
// //   std::string sPbName = pbName;
// //
// // //   if (sPbName == "Poisson")                        // Laplace 1D 2D 3D
// // //     _problem = new LaplaceProblem(_comm);
// // //   else if (sPbName == "Elasticity")                // Elaticity
// // //     _problem = new LinearElasticProblem(_comm);
// // //   else if (sPbName == "convdiff")                  // Heat equation
// // //     _problem = new Con_DiffProblem(_comm);
// // //   else if (sPbName == "navierstokes")              // Navier-Stokes 2D
// // //     _problem = new NS_Problem(_comm);
// // //   else if (sPbName == "navierstokes3D")            // Navier-Stokes 3D
// // //     _problem = new NS_Problem3D(_comm);
// // //   else if (sPbName == "monodim")                   // monodimensional
// // //     _problem = new Monodim_Problem(_comm);
// // //   else
// // //     _problem = NULL;
  return;
}

// ========================================================================
// This function gets the Group mesh from the med file (medfile_name)
// and sets the id and the mesh to the interfaces
// vector (_interface_mesh_vect[i])
// through the index of the interface-functions
// (from EquationSystemsExtended)
void FEMUS::init_interfaces(
  const std::string & medfile_name, // medfile name    (in)
  const int index_mgmesh,           // mgmesh index    (in)
  const int index_medmesh                       // med-mesh index  (in)
) {// =========================================================================


  // Reading mesh Names from med file  ------------------------------
  std::string mesh_dir=_mg_utils->_mesh_dir;
  std::string localFile=mesh_dir+medfile_name;
  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(localFile.c_str());
  // Reading mesh from the index_medmesh med-mesh name
  if (meshNames.size() < 1) {
    std::cout<<  " FEMUS::setMesh : no meshes in the file'";
  }
  std::string localMeshName = meshNames[ index_medmesh];
  // Reading group names
  std::vector<std::string> GroupNames =
    MEDLoader::GetMeshGroupsNames(localFile.c_str(), localMeshName.c_str());

  init_BCzone(GroupNames);

  int nGroups0 = GroupNames.size();
//   int nGroups=0;
//   for (int i=0; i< nGroups0; i++)
//     if (atoi(GroupNames[i].c_str())>=10) nGroups++;

  // From goup names to interface_mesh (id,name,support)

  std::vector<std::string> vG(1);
  _interface_mesh_vect.resize(nGroups0);
//   int iGroups=0;
  for (int i=0; i<nGroups0; i++) {
//     if ( _boundary[i].type == "Mark") {

    vG[0] = GroupNames[i];
    _interface_mesh_vect[i].name = GroupNames[i]; // interface_mesh name
    int id_name=atoi(GroupNames[i].c_str());            // id group name (int)
    _interface_mesh_vect[i].id =id_name;
    int id_level=-1; bool on_nodes=true; 
    if (id_name<10) {id_level=0; on_nodes=false;}
    _interface_mesh_vect[i].support = NULL;  
    if (_boundary[i].type == "Mark") {
      _interface_mesh_vect[i].support =             // interface_mesh support
        MEDLoader::ReadUMeshFromGroups(
          localFile.c_str(), meshNames[0].c_str(), id_level, vG);
      if (_interface_mesh_vect[i].support) {             }// interface_mesh identity
      _interface_mesh_vect[i].support->zipCoords();
//      _interface_mesh_vect[iGroups].id = defineBoundary(_interface_mesh_vect[i].support);
     
      std::cout << "FEMUS::setInterfaces: support  set to boundary with name "<<id_name  <<"\n";
      _mg_equations_map[0]->addfunction_BC(
                id_name, _interface_mesh_vect[i].support,
		_boundary[i].from_cmp,_boundary[i].n_cmp,_boundary[i].order_cmp,
		on_nodes);

//       iGroups++;
    } else {
      _interface_mesh_vect[i].id = -1;
//        _interface_mesh_vect[iGroups].name =
      std::cout << "FEMUS::setInterfaces: support not set to boundary "<<id_name  <<"\n";
    }


  } // for (int i



  return;
}
// //
// //
// //
// =============================================================================
// This function sets the mesh from med-mesh (m) to libmesh
void FEMUS::setMesh(
//   const ParaMEDMEM::MEDCouplingUMesh * m
) {// ==========================================================================
//   if (_problem) _problem->setMesh(m);


  //MGMesh  (std::vector<MGMesh*> mg_mesh)
//   for (int mesh=0; mesh<_num_mgmesh; mesh++) {
//     const double Lref  = _mg_phys[mesh]->get_par("Lref");      // reference L
    const int NoLevels= _mg_utils->get_par("nolevels");  // numb of Level
    _mg_mesh.push_back(new ParallelMeshExtended(*_med_mesh,_start->comm(), *_mg_utils,*_mg_geomel,NULL));
//     _mg_mesh.push_back(new MGMesh(_start->comm(), *_mg_utils[mesh],*_mg_geomel,1.));
    if (NoLevels != _mg_mesh[0]->_NoLevels) {
      std::cout << "Inconsistent Number of Levels between Mesh and SolBase"
                << std::endl;
      abort();
    }
    // prind mesh at level NoLevels-1 (linear connectivity)
    _mg_mesh[0]->print(NoLevels-1,0);

    // prind mesh at level NoLevels-1 (med format)
    std::string mesh_name= _mg_utils->get_file("F_MESH_READ");//://= _mg_utils.get_file("F_MESH_READ");
    unsigned pos = mesh_name.find(".");         // position of "live" in str
    std::ostringstream name;
    name << _mg_utils->_mesh_dir <<  mesh_name.substr(0,pos) << "_fine.med" ;
   // _mg_mesh[mesh]->print_med(NoLevels-1,name.str().c_str());
//   }

  return;
}

// =============================================================================
//This function reads and sets the med-mesh
// (also the libmesh calling the other setMesh function)
void FEMUS::setMedMesh(const std::string & dataFile) {

  // processor
//    int iProc=0; MPI_Comm_rank(_start->comm(), &iProc);
  // set subdomain number
//   int i=0;
// #ifdef HAVE_PETSCM
//   MPI_Comm_rank(MPI_COMM_WORLD, &i);
//   std::cout <<  "FEMUS::setMesh:  proc: " << i << std::endl;
// #endif



  // Reading MED mesh ----------------------------------------------------
  // Mesh  filename
  std::string localFile;
  int l = dataFile.size();
  if (dataFile.substr(l-4) == ".med")    {localFile = dataFile;}
  else {
    std::cout<<"FEMUS::setMesh:"<<  dataFile <<"does not exist!";
  }
  std::cout << " FEMUS::setMesh: MED file "  << ": "<< localFile << std::endl;

  // Mesh names (inside MEDfile one can have different meshes
  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(localFile.c_str());
//   std::vector<std::string> meshNames = MEDLoader::GetMeshNames("/homesd/msandro/software/femus/USER_APPL/MESH/test1quad9_group_mat_gen.med");
  if (meshNames.size() < 1) {
    std::cout<<  " FEMUS::setMesh : no meshes in the file'";
  }
  // The first name is the good one
  std::string localMeshName = meshNames[0];
  _med_mesh = MEDLoader::ReadUMeshFromFile(localFile.c_str(), localMeshName.c_str(), 0);
  if (_med_mesh == NULL) {
    std::cout<<  " FEMUS::setMesh : unable to read the med-mesh'";
  }




  return;
}

// *******************************************************************
// ********************* set and get  boundary conditions ************

// ===================================================================
// This function gets the name list of interfaces
std::vector<std::string> FEMUS::get_interface_Names(
) { //================================================================
  std::vector<std::string> names(_interface_mesh_vect.size());
  for (int i = 0; i<(int) _interface_mesh_vect.size(); i++)
    names[i] =  _interface_mesh_vect[i].name;
  return names;
}
// ===================================================================
// This function gets the name of interface i
std::string FEMUS::get_interface_Name(
  int i  //  vector index
) { // ===============================================================
  assert(i>=0 && i<(int)_interface_mesh_vect.size());
  return _interface_mesh_vect[i].name;
}

// // ============================================================================
// // This functions set the boundary conditions to  EquationSystemsExtended
// void FEMUS::setBCType(
//   const std::string & nameBC,   // boundary condition name (in)
//   const std::string & typeBC    // boundary condition type (in)
// ) { // ========================================================================
//   // get the identity
//   int id =search_idBC_byname(nameBC);
//
// //   _mg_equations_map[0]->setBoundaryConditionType(id, typeBC.c_str());
// //    _mg_equations_map[0]->setBC(id, "0.0");
// //   if (_problem) {
// //     _problem->setBoundaryConditionType(id, typeBC.c_str());
// //     _problem->setAnalyticBoundaryValues(id, "0.0");
// //  }
// }
// ===================================================================
void FEMUS::setAnalyticBoundaryValues(
  const std::string & bcName,           // boundary name
   int n_cmp,
  const std::string & /*bcType*/,           // boundary type
  const std::string & bcExpression      // boundary symbolic expr
) {
//   if (FDEBUG) fDebugPos << bcName << " " << bcType << std::endl;
//   setBCType(bcName, bcType);
  setAnalyticBCValues(bcName,n_cmp, bcExpression);
  return;
}


void FEMUS::setFieldBoundaryValues(
  const std::string & bcName,
 int n_cmp,
 const std::string & /*bcType*/,
 const ParaMEDMEM::MEDCouplingFieldDouble * bcField) {
//   if (FDEBUG) fDebugPos << bcName << " " << bcType << std::endl;
//   setBCType(bcName, bcType);
  if (bcField ==NULL) return;
  setFieldBCValues(bcName,n_cmp, bcField);
  return;
}
// ============================================================================
void FEMUS::setAnalyticBCValues(
  const std::string & nameBC,
  int n_cmp,
  const std::string & f
) {
  int id = search_idBC_byname(nameBC); 
  _mg_equations_map[0]->setBC(id,n_cmp , f.c_str());
//   if (_problem) _problem->setAnalyticBoundaryValues(id, f.c_str());
}
// =================================================
void FEMUS::setFieldBCValues(
  const std::string & nameBC,
  int n_cmp,
  const ParaMEDMEM::MEDCouplingFieldDouble *f) {
  int id = search_idBC_byname(nameBC);
  _mg_equations_map[0]->setFieldBoundaryValues(id,n_cmp, f);

  // set field to equation system
//   if (f->getTypeOfField() == ParaMEDMEM::ONCELLS)
//     _mg_equations_map[0]->setBC(id, f);
//    else
//   _mg_equations_map[0]->setNodeBC(id, n_cmp, f);
//   if (_problem) _problem->setFieldBoundaryValues(id, f);
  return;
}

// int FEMUS::defineBoundary(const ParaMEDMEM::MEDCouplingUMesh * mesh_bc) {
// //    return (_problem) ? _problem->defineBoundary(mesh_bc) : -1;
//  int j=_mg_equations_map[0]->defineIdInterface(mesh_bc);
//
//   return j;
// //   return NULL;
// }




// ============================================================================
/// This function gets the values of the variable with names in
/// std::vector<char *> vName on boundary with identity id
ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::getValuesOnBoundary(// field (out) 
  const std::string & nameBC,                    // boundary name (char*) (in)
  const std::string & systemName                      // system name      (in)
)  {
  int id = search_idBC_byname(nameBC);
  return _mg_equations_map[0]->getValuesOnBoundary_nodes(id, systemName.c_str());
}


// ============================================================================
/// This function gets all the values on boundary with identity id
ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::getValuesOnBoundary(  // field (out) 
  const std::string & nameBC,                     // boundary name (char*) (in)
 const std::string & systemName,                  // system name           (in)
  int iCmp                                        // component             (in)
)  {
  int id =search_idBC_byname(nameBC);
  return _mg_equations_map[0]->getValuesOnBoundary_nodes(id,systemName.c_str(),iCmp);
}

// ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::getValuesOnBoundary_nodes
// (const std::string & nameBC,  std::vector<char *> vName) const {
//   int id = searchBoundaryCondition(nameBC);
//   return _problem ? _problem->getValuesOnBoundary_nodes(id, vName) : NULL;
// }






// *******************************************************************
// **************** source *******************************************

// ===================================================================
// This function sets an analytic source
void FEMUS::setAnalyticSource(
  const std::string & f                       // analytic source
) { // ================================================================
//     if (_problem) _problem->setSource(_mesh,f.c_str());

  _mg_equations_map[0]->setSource(_med_mesh,f.c_str());
}
// ===================================================================
// This function sets an analytic source
// void FEMUS::setAnalyticSource(
// //   const  ParaMEDMEM::MEDCouplingUMesh * medmesh,
//   const std::string & f                       // analytic source
// ) { // ================================================================
//    if (_problem) _problem->setSource(medmesh,f.c_str());
// //   if (_mg_equations_map[0]) _mg_equations_map[0]->setSource(medmesh,f.c_str());
// }
// ===================================================================
// This function sets a numerical source (coupling)
void FEMUS::setSource(
  const ParaMEDMEM::MEDCouplingFieldDouble *f  // numerical sources
) { // ================================================================
  _mg_equations_map[0]->setSource(_med_mesh,f);
//  if (_problem) _problem->setSource(f);
}

// ===================================================================
// This function sets a numerical source (coupling)
void FEMUS::setAverageSources(
  int /*n*/,              //  number of average sources <-
  double /*val*/[]        //  average source            <-
) {
//   if (_problem) _problem->setAverageSources(n,val);
}

// // // *******************************************************************
// // // *******************************************************************
// // // **************** Solve  *******************************************
//====================================================================
// This function solves the problem
void FEMUS::solve() {
//    if (_problem) _problem->solve();

  int      t_in=0;
  double    time=0.;                        //  initial time
  const int restart      = _mg_utils->get_par("restart");   // restart or not

  _mg_time_loop->transient_setup(restart,t_in,time);             //  MGTimeLoop: setup
  _mg_time_loop->transient_loop(t_in,time);


}
// //
// // // *******************************************************************
// // // **************** Field solution************************************
// //
// // //====================================================================
// // //====================================================================
// // // This function gets the solution in MED format
// // ParaMEDMEM::MEDCouplingFieldDouble *FEMUS::getOutputField(
// //   const std::string & vName         // field name <-
// // ) const { //==========================================================
// // //   return (_problem) ? _problem->getOutputField(vName.c_str()) : NULL;
// // }
// //
// // //====================================================================
// // // This function prints the solution in MED format
// // void FEMUS::saveOutputField(
// //   const std::string & vName,                // field name <-
// //   const std::string & prefix                // file name  <-
// // ) const { // ==========================================================
// // //   ParaMEDMEM::MEDCouplingFieldDouble * f = getOutputField(vName);
// // //   std::string t = replaceEnv(prefix);
// //
// //   int iProc;  MPI_Comm_rank(_comm, &iProc);
// // //   std::ostringstream s;
// // //   s << prefix << "_" << vName << "_" << iProc+1 << ".med";
// // //   MEDLoader::WriteField(s.str().c_str(), f, true);
// //   MPI_Barrier(_comm);
// // }
// //
// //
// //
// //
// // //====================================================================
// // // This function print the solution in general format
// // ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::getInputFieldTemplate(const std::string &name) {
// // //    if (name == "")
// // //       return getOutputField("???");
// // //    else
// // //       return getValuesOnBoundary(name, "???");
// // }
// //
// //
// //
// // // *******************************************************************
// // // **************** Debug *******************************************
// //
// // // ===================================================================
// // // This function sets the file for logging
// // void FEMUS::setDebugFile(
// //   const std::string & s
// // ) { // ===============================================================
// // //   if (FDEBUG) fDebug.open(s.c_str(), _comm);
// // }
// //
// // // ===================================================================
// // // This function ends the  logging
// // void FEMUS::endDebug(
// // ) {// ===============================================================
// // //   fDebug.close();
// // }


// ============================================================================
// ==========================  Private ========================================
// ============================================================================

// ============================================================================
// This function gives the id in the _interface_mesh_vect from
// the name boundary
int FEMUS::search_idBC_byname(        // id name (return)
  const std::string & nameBC          //  boundary name  (in)
) const { // ====================================================

  for (std::vector<interface_mesh>::const_iterator
       i = _interface_mesh_vect.begin(); i != _interface_mesh_vect.end(); i++)
    if (i->name == nameBC) return i->id;

  std::cout<< " Non existent boundary region named in that way: "<< nameBC;
  abort();

}
// ============================================================================
// This function gives the id in the _interface_mesh_vect from
// the name boundary
int FEMUS::search_idxBC_byid( // vector index              (return)
  const int  inameBC          //  id boundary name (int)   (in)
) const { // ==================================================================
  int icount=0;
  for (std::vector<interface_mesh>::const_iterator
       i = _interface_mesh_vect.begin(); i != _interface_mesh_vect.end(); i++) {
    if (i->id == inameBC) return icount;
    icount++;
  }
  std::cout<< " Non existent boundary region with id in that way: "<< inameBC;
  abort();

}

// ============================================================================
// This function gives the id in the _interface_mesh_vect from
// the name boundary
void FEMUS::set_interface_userbc() {// vector index
  // Boundary Probelm 1 (mesh1)
  for (int i=0; i<(int)_boundary.size(); i++) {
    UsersBC & b = _boundary[i];
    if (b.type=="Mark") {
      if (b.isAnalytic)
        setAnalyticBoundaryValues(b.name,b.n_cmp, b.type, b.equation);
      else
        setFieldBoundaryValues(b.name,b.n_cmp, b.type, b.field);
    }
  }
  return;
}
