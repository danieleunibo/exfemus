#include <iostream>
#include <cstdlib>
#include <sstream>


// configuration files -------------------------
#include   "Printinfo_conf.h"

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
#include "FEMUS.h"

#ifdef HAVE_MED
// MED includes
#include "InterfaceFunctionM.h"
#include "MEDLoader.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#endif

#include "MeshExtended.h"
#include "EquationSystemsExtendedM.h"


// ****************************************************************************
// ****************  Constructor Destructor ***********************************

// ============================================================================
// Basic constructor
FEMUS::FEMUS()  :
  _comm(MPI_COMM_WORLD)    // communicator
{
  // Init MPI flag --------------------------------
  int flag=0;  MPI_Initialized(&flag);
  if(flag) {_local_MPI_Init = false; }
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
  _comm(comm)   // communicator
//   _num_mgmesh(0)
  {        // n of femus_meshes
  // transient system
  // Init MPI flag
  int flag=0;  MPI_Initialized(&flag);
  if(flag) {_local_MPI_Init = false; }
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
// ============================================================================
void FEMUS::init_fem(
  MGGeomEl & mggeomel,
  MGFEMap & mgfemap
) { // ========================================================================
// A) setting MGGeomEl
  _mg_geomel=&mggeomel;  // ***************************************************
  if(_mg_geomel == NULL) {
    std::cout<< "FEMUS::init_fem: no _mg_geomel"; abort();
  }
  /// B) setting MGFEMap (fem)
  _mg_femap=&mgfemap;  // *****************************************************
  if(_mg_femap == NULL) {
    std::cout<< "FEMUS::init_fem: no _mg_femap"; abort();
  }

  return;
}
// ============================================================================
// This function is the destructor
FEMUS::~FEMUS() {
  // ==========================================================================
  delete _start;
  #ifdef HAVE_MED
//   if(_med_mesh) _med_mesh->decrRef();        // med-mesh
#endif

}

// // ============================================================================
// // This function is the problem destructor
// void FEMUS::terminate(
// ) {// =========================================================================
// 
// }

// // // ****************************************************************************
// // // ****************    end Constructor Destructor *****************************
// //
// // // ****************************************************************************
// // // ****************    Set    *************************************************


// // // =============================================================================
// // // This function sets the type of problem
void FEMUS::setSystem(
  const std::vector<FIELDS> & pbName,
  int n_data_points,
  int n_data_cell
) {// ==========================================================================


  _mg_equations_map=new EquationSystemsExtendedM(*_mg_utils,*_mg_mesh,*_mg_femap,n_data_points,n_data_cell);  // MGEquationsMap class
  _mg_equations_map->read_par();
#ifdef PRINT_INFO  // ---- info ---------------
  _mg_equations_map->print_par();       // print parameters
#endif
  _mg_equations_map->init_data(0);
  _mg_equations_map->init(pbName);                              // adds the equations to the map
  _mg_equations_map->setDofBcOpIc();                            // set operators
  _mg_equations_map->set_mesh_mg(*_mg_mesh);
#ifdef HAVE_MED
  _mg_equations_map->set_mesh_med(*_med_mesh);
#endif
//   }
  if(_mg_geomel == NULL) {
    std::cout<< "FEMUS::setSystem: no _mg_equations_map"; abort();
  }

  //time loop
  _mg_time_loop=new  MGTimeLoop(*_mg_utils,*_mg_equations_map);
  if(_mg_time_loop == NULL) {
    std::cout<< "FEMUS::setSystem: no _mg_time_loop"; abort();
  }

  return;
}


// =============================================================================
// This function sets the mesh from med-mesh (m) to libmesh
void FEMUS::setMesh(
) {// ==========================================================================

  const int NoLevels= _mg_utils->get_par("nolevels");  // numb of Level
  _mg_mesh=new MeshExtended(_start->comm(), *_mg_utils,*_mg_geomel);
  // check insanity
  if(_mg_mesh == NULL) {
    std::cout<< "FEMUS::setMesh: no _mg_mesh"; abort();
  }
  if(NoLevels != _mg_mesh->_NoLevels) {
    std::cout << "Inconsistent Number of Levels between Mesh and SolBase"
              << std::endl;
    abort();
  }
  // prind mesh at level NoLevels-1 (linear connectivity)
  _mg_mesh->print(NoLevels-1,0);

#ifdef HAVE_MED
  // prind mesh at level NoLevels-1 (med format)
  std::string mesh_name= _mg_utils->get_file("F_MESH_READ");//://= _mg_utils.get_file("F_MESH_READ");
  unsigned pos = mesh_name.find(".");         // position of "live" in str
  std::ostringstream name;
  if(_mg_mesh->_iproc==0) {
    name << _mg_utils->_mesh_dir <<  mesh_name.substr(0,pos) << "_fine.med" ;
    _mg_mesh->print_med(NoLevels-1,name.str().c_str());
  }
#endif

  return;
}

// // // *******************************************************************
// // // **************** Solve  *******************************************
// // // *******************************************************************
/// This function sets up the intial set
void FEMUS::solve_setup(
  int        & t_in,                 ///< initial time iteration
  double     &  time                 ///< actual time
) {
  const int restart      = _mg_utils->get_par("restart");   // restart or not
  _mg_time_loop->transient_setup(restart,t_in,time);     //  MGTimeLoop: setup
  return;
}

//=============================================================================
// This function solves one step  for transient problems
void FEMUS::solve_onestep(
  const int  & t_in,                 ///< initial time iteration
  const int  & t_step,               ///< actual time iteration
  const int  & print_step,            ///< print every
  double     &  time,                ///< actual time
  double     &  dt                   ///< step time
) { // ========================================================================
  _mg_time_loop->transient_onestep(t_in,t_step,print_step,time,dt);    ///< step time
  return;
}

#ifdef HAVE_MED
// ========================================================================
// This function gets the Group mesh from the med file (medfile_name)
// and sets the id and the mesh to the interfaces
// vector (_interface_mesh_vect[i])
// through the index of the interface-functions
// (from EquationSystemsExtendedM)
void FEMUS::init_interface(
  const int interface_name,
  int interface_id,
  int order_cmp,
  const std::string & medfile_name, // medfile name    (in)
  bool on_nodes,
  const int index_medmesh                       // med-mesh index  (in)
) {// =========================================================================

  ostringstream name_id; name_id <<interface_id;
  std::vector<std::string> vG(1);
  std::string id_name=name_id.str().c_str();
  vG[0] =name_id.str().c_str();

  // Reading mesh Names from med file  ------------------------------

  std::string mesh_dir=_mg_utils->_mesh_dir;
  std::string localFile=mesh_dir+medfile_name;
  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(localFile.c_str());
  // Reading mesh from the index_medmesh med-mesh name
  if(meshNames.size() < 1) {
    std::cout<<  " FEMUS::setMesh : no meshes in the file'";
  }
  std::string localMeshName = meshNames[index_medmesh];
  // Reading group names
  std::vector<std::string> GroupNames =
    MEDLoader::GetMeshGroupsNames(localFile.c_str(), localMeshName.c_str());


  int nGroups0 = GroupNames.size();
  ParaMEDMEM::MEDCouplingUMesh * support;
  int id_level=-1;  if(interface_id<10) {id_level=0;}
  support = MEDLoader::ReadUMeshFromGroups(localFile.c_str(), meshNames[0].c_str(), id_level,vG);
  support->zipCoords();

  std::cout << "FEMUS::setInterfaces: support  set to boundary with name "<<
            interface_id  << "\n";
  _mg_equations_map->add_interface_fun(interface_name, interface_id, support,on_nodes,order_cmp);


  return;
}


// ========================================================================
// This function gets the Group mesh from the med file (medfile_name)
// and sets the id and the mesh to the interfaces
// vector (_interface_mesh_vect[i])
// through the index of the interface-functions
// (from EquationSystemsExtendedM)
void FEMUS::init_interface(
  const int interface_name,
  int interface_id,
  int order_cmp,
  const std::string & medfile_name, // medfile name    (in)
  const std::string & medfile_name_old, // medfile name    (in)
  const FEMUS & P_old,
  const int index_mgmesh,           // mgmesh index    (in)
  const int index_medmesh                      // med-mesh index  (in)

) {// =========================================================================


  ostringstream name_id; name_id <<interface_id;
  std::vector<std::string> vG(1);
  std::string id_name=name_id.str().c_str();
  vG[0] =name_id.str().c_str();

  // Reading mesh Names from med file  ------------------------------
  std::string mesh_dir=_mg_utils->_mesh_dir;
  std::string localFile=mesh_dir+medfile_name;
  std::string localFile_old=mesh_dir+medfile_name_old;
  std::vector<std::string> meshNames = MEDLoader::GetMeshNames(localFile.c_str());
  // Reading mesh from the index_medmesh med-mesh name
  if(meshNames.size() < 1) {
    std::cout<<  " FEMUS::setMesh : no meshes in the file'";
  }
  std::string localMeshName = meshNames[index_medmesh];
  // Reading group names
  std::vector<std::string> GroupNames =
    MEDLoader::GetMeshGroupsNames(localFile.c_str(), localMeshName.c_str());

  int nGroups0 = GroupNames.size();

  // From group names to interface_mesh (id,name,support)
  ParaMEDMEM::MEDCouplingUMesh * support;
  int id_level=-1; bool on_nodes=true;
  if(interface_id<10) {id_level=0; on_nodes=false;}
  support= MEDLoader::ReadUMeshFromGroups(localFile_old.c_str(), meshNames[index_medmesh].c_str(), id_level, vG);
  support->zipCoords();

  std::cout << "FEMUS::setInterfaces: support  set to boundary with name "<<id_name  <<"\n";
  _mg_equations_map->add_interface_fun(interface_name,interface_id, support,on_nodes,order_cmp);

  return;
}

// =========================================================================
/// This function sets the value from first_cmp to end_cmp of the field
/// on the old solution x_old of the interface id
void FEMUS::write_Boundary_value(
  int id_boundary_name ,      ///< identity interface name (in)
  std::string mgsystem_name, ///< system name          (in)
  int n_cmp,             ///< from variable system (in)
  int first_cmp               ///< to variable system   (in)

) {
  _mg_equations_map->write_Boundary_value(
    id_boundary_name,
    mgsystem_name,n_cmp,first_cmp);
  return;
}


// =============================================================================
//This function reads and sets the med-mesh
// (also the libmesh calling the other setMesh function)
void FEMUS::setMedMesh(const std::string & dataFile) {

  // Reading MED mesh ----------------------------------------------------
  // Mesh  filename
  std::string localFile;
  int l = dataFile.size();
  if(dataFile.substr(l-4) == ".med")    {localFile = dataFile;}
  else {
    std::cout<<"FEMUS::setMesh:"<<  dataFile <<"does not exist!";
  }
  std::cout << " FEMUS::setMesh: MED file "  << ": "<< localFile << std::endl;

  // Mesh names (inside MEDfile one can have different meshes
   std::vector<std::string> meshNames = MEDLoader::GetMeshNames(localFile.c_str());
//   std::vector<std::string> meshNames = MEDLoader::GetMeshNames("/homesd/msandro/software/femus/USER_APPL/MESH/test1quad9_group_mat_gen.med");
   if(meshNames.size() < 1) {
     std::cout<<  " FEMUS::setMesh : no meshes in the file'";abort();
   }
  // The first name is the good one
  std::string localMeshName = meshNames[0];
  _med_mesh = MEDLoader::ReadUMeshFromFile(localFile.c_str(), localMeshName.c_str(), 0);
  if(_med_mesh == NULL) {
    std::cout<<  " FEMUS::setMesh : unable to read the med-mesh'"; abort();
  }

  return;
}



// ===================================================================
void FEMUS::setAnalyticBoundaryValues(
  int interface_name,
  int n_cmp,
  const std::string & bcExpression      // boundary symbolic expr
) {
  _mg_equations_map->setBC(interface_name,n_cmp , bcExpression.c_str());
  return;
}

// ===================================================================
void FEMUS::setAnalyticSource(
//   const std::string & bcName,           // boundary name
  int interface_name,
  int n_cmp,
  const std::string & bcExpression      // boundary symbolic expr
) {
  _mg_equations_map->setBC(interface_name,n_cmp , bcExpression.c_str());
  return;
}

void FEMUS::setFieldBoundaryValues(
  int interface_name,
  int n_cmp,
  const ParaMEDMEM::MEDCouplingFieldDouble * bcField) {

  if(bcField ==NULL) return;

  _mg_equations_map->setBC(interface_name,n_cmp, bcField);

  return;
}

void FEMUS::setFieldSource(
  int interface_name,
  int n_cmp,
  const ParaMEDMEM::MEDCouplingFieldDouble * srcField) {

  if(srcField ==NULL) return;
  _mg_equations_map->setBC(interface_name,n_cmp, srcField);

  return;
}


// ============================================================================
/// This function gets all the values on boundary with identity id
ParaMEDMEM::MEDCouplingFieldDouble * FEMUS::getValuesOnBoundary(  // field (out)
  int  interface_name,                     // boundary name (char*) (in)
  const std::string & systemName,                  // system name           (in)
  int n_cmp,                                        // component             (in)
  int first_cmp                                        // component             (in)
)  {
  return _mg_equations_map->getValuesOnBoundary_nodes(interface_name,systemName.c_str(),n_cmp,first_cmp);
}
#endif