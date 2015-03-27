#include <set>
#include <cstring>

#include "EquationSystemsExtendedM.h"
#include "MeshExtended.h"
#include "numeric_vectorM.h"

#ifdef HAVE_MED
#include "InterfaceFunctionM.h"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"

// using namespace libMesh;
using namespace ParaMEDMEM;
#endif
// ========================================================================
/// Constructor
EquationSystemsExtendedM::EquationSystemsExtendedM(
  MGUtils& mgutils_in,
  MeshExtended& mgmesh_in,
  MGFEMap& mgfemap_in,
  int npoint_data,
  int ncell_data
) :
  MGEquationsSystem(mgutils_in,mgmesh_in,mgfemap_in,npoint_data,ncell_data),
  _mg_mesh(&mgmesh_in) {  // mesh class in

}


// ========================================================================
EquationSystemsExtendedM::~EquationSystemsExtendedM() {
#ifdef HAVE_MED
  std::map<int, InterfaceFunctionM *>::iterator it = _interfaceFunMap.begin();
  for(; it != _interfaceFunMap.end(); it++) {
    if(it->second) delete it->second;
    it->second = NULL;
  }
#endif
return;
}


#ifdef HAVE_MED
// =======================================================================
void EquationSystemsExtendedM::setBC(
  int name,
  int  n_cmp,
  const MEDCouplingFieldDouble *field
) {// =======================================================================

  InterfaceFunctionM * fct = get_interface_fun(name);  // interface-function
  if(fct == NULL) return;
  fct->set_field(field);
#ifdef PRINT_MED
//   fct->printOn(std::cout,name);
#endif

  return;
}

// =========================================================================
void EquationSystemsExtendedM::setBC(
  int name,
  int n_cmp,
  const char *s
) { // ======================================================================

  InterfaceFunctionM * fct = get_interface_fun(name);
  if(fct == NULL) return;
  if(name>9)
    fct->set_analytic_field(s, n_cmp);
  else
    fct->set_analytic_field_elem(s, n_cmp);
  
#ifdef  PRINT_MED
//  fct->printOn(std::cout,name);
#endif

  return;
}


// ============================================================================
// ===============  Interface function routines ===============================

// ============================================================================
// This functon add an interface-function (med-mesh,femus-mesh, field)
//  in the interface-function boundary map (_interfaceFunMap):
//  the field in the interface-function is not assigned here
void EquationSystemsExtendedM::add_interface_fun(       ///<map position(return)
  const int  interface_name,
  const int  interface_id,    ///< boundary id (int)  (in)
  const ParaMEDMEM::MEDCouplingUMesh * b,  ///< med-submesh        (in)
//   const int /*from_cmp*/,               ///< initial id component
//   const int /*n_cmp*/,                  ///< n components
  const bool on_nodes,                     ///< values on nodes (true)
  const int order_cmp                      ///< order component (2=quad;lin=1)
) { // ========================================================================

  // new function
  InterfaceFunctionM *fun = new InterfaceFunctionM;
  // set the mesh interface
  if(on_nodes) {
    fun->set_mesh_interface_nodeID(_mg_mesh,b,interface_id,order_cmp);
  } else { // on elements (mat)
    fun->set_mesh_interface_elemID(_mg_mesh,b,interface_id);
  }
  // setting the fun in the map _interfaceFunMap at bd_name_id
  _interfaceFunMap[interface_name] = fun;  // added to the map

  // print
#ifdef PRINT_MED
  std::cout << " Added  InterfaceFunction on interface "<< interface_id << "\n";
  fun->printOn(std::cout, interface_id);
#endif
  return; 
}

// ============================================================================
// This function gets the  the value of the variable with id number
//  "variable_id" on nodes on the boundary with identity "id" in the
//  system "system_name"
MEDCouplingFieldDouble * EquationSystemsExtendedM::getValuesOnBoundary_nodes(
  int                name,           // int boundary identity   (in)
  const char *system_name,           // system name             (in)
  int               n_cmp,           //  first variable       (in)
  int           first_cmp            // n variables       (in)
)  {// ========================================================================
  // from LibMesh function fct
  InterfaceFunctionM * fct = get_interface_fun(name);
  if(fct == NULL)  return NULL;
  int nNodes = fct->getSupport()->getNumberOfNodes();
//   std::map<int, int> & NodeID = fct->NodeID();
  int n_nodes_mg  = fct->get_n();
  int * map_mg  = fct->get_map_mg();
  int * map_med = fct->get_map_med();

  // from MGMesh
  int Level=_mg_mesh->_NoLevels-1;
  int offset=_mg_mesh->_NoNodes[Level];
  // from  MGSystem* system;
  MGSolBase * mgsyst=get_eqs(system_name);
//   int nComp =n_variable;
//   std::cout << "\n EquationSystemsExtendedM nComp = " << n_cmp << std::endl;

  // A new LibMesh function f
  MEDCouplingFieldDouble * f = MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES);
  f->setMesh(fct->getSupport());
  f->setName(system_name);

  // array function to fill f
  DataArrayDouble *array = DataArrayDouble::New();
  array->alloc(nNodes,n_cmp);
  for(int i_mg= n_nodes_mg; i_mg <nNodes; i_mg++) 
    for(int j=  first_cmp; j< first_cmp+n_cmp; j++) 
       array->setIJ(i_mg,j-first_cmp,0.);
  // filling array
  for(int i_mg=0; i_mg < n_nodes_mg; i_mg++) {
    int node_mg   = map_mg[i_mg];  // mg  node
    int node_med  = map_med[i_mg];  // med node

    for(int j=  first_cmp; j< first_cmp+n_cmp; j++) {
      const int kdof_top = mgsyst->_node_dof[Level][node_mg+j*offset];
      double v = (*(mgsyst->x_old[Level]))(kdof_top);
      array->setIJ(node_med,j-first_cmp,v);
    }
  }
  f->setArray(array);   // set array -> f
  array->decrRef();     // delete array
  f->checkCoherency();  // check f
//  fct->printOn(std::cout,name);   // print fct

  return f;

}

// ============================================================================
// This function gets the  the value of the variable with id number
//  "variable_id" on nodes on the boundary with identity "id" in the
//  system "system_name"
MEDCouplingFieldDouble * EquationSystemsExtendedM::getValuesOnBoundary_elem(
  int                  id,           // int boundary identity   (in)el_conn
  const char *system_name,           // system name             (in)
  int               n_cmp,           //  first variable       (in)
  int           first_cmp            // n variables       (in)
)  {// ========================================================================
  // from LibMesh function fct
  InterfaceFunctionM * fct = get_interface_fun(id);
  if(fct == NULL)  return NULL;
  int nNodes = fct->getSupport()->getNumberOfNodes();
//   std::map<int, int> & NodeID = fct->NodeID();
  int n_nodes_mg  = fct->get_n();
  int * map_mg  = fct->get_map_mg();
  int * map_med = fct->get_map_med();

  // from MGMesh
  int Level=_mg_mesh->_NoLevels-1;
  int offset=_mg_mesh->_NoNodes[Level];
  // from  MGSystem* system;
  MGSolBase * mgsyst=get_eqs(system_name);
//   int nComp =n_variable;
  std::cerr << "\n EquationSystemsExtendedM nComp = " << n_cmp << std::endl;

  // A new LibMesh function f
  MEDCouplingFieldDouble * f = MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES);
  f->setMesh(fct->getSupport());
  f->setName(system_name);

  // array function to fill f
  DataArrayDouble *array = DataArrayDouble::New();
  array->alloc(nNodes,n_cmp);
  // filling array
  for(int i_mg=0; i_mg < n_nodes_mg; i_mg++) {
    int node_mg   = map_mg[i_mg];  // mg  node
    int node_med  = map_med[i_mg];  // med node

    for(int j=  first_cmp; j< first_cmp+n_cmp; j++) {
      const int kdof_top = mgsyst->_node_dof[Level][node_mg+j*offset];
      double v = (*(mgsyst->x_old[Level]))(kdof_top);
      array->setIJ(node_med,j,v);
    }
  }
  f->setArray(array);   // set array -> f
  array->decrRef();     // delete array
  f->checkCoherency();  // check f
//   fct->printOn(std::cout,id);   // print fct

  return f;

}
// =========================================================================
void EquationSystemsExtendedM::write_Boundary_value(
  int id_boundary_name,       ///< identity interface name (in)
  std::string mgsystem_name, ///< system name             (in)
  int n_cmp,                 ///<  variable system    (in)
  int first_cmp              ///< from variable system      (in)

) { // ======================================================================

  InterfaceFunctionM * fct = get_interface_fun(id_boundary_name);  // interface-function
  if(fct == NULL)  return;

  const ParaMEDMEM::MEDCouplingUMesh *support= fct->getSupport();
  int nNodes_med = support->getNumberOfNodes();
  int n_nodes_mg  = fct->get_n();
  int * map_mg  = fct->get_map_mg();
  int * map_med = fct->get_map_med();


  MGSolBase * mgsyst=get_eqs(mgsystem_name.c_str());
  int Level=_mg_mesh->_NoLevels-1;
  int offset= _mg_mesh->_NoNodes[Level];
  NumericVectorM *old_sol_top= mgsyst->x_old[Level];


  MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *>(&_mgmesh);
  std::vector<int> mat=ext_mesh->_mat_id;
  std::vector<int> bc_id=ext_mesh->_bc_id;

//   std::vector<double> /*src_value*/(3);
  double bc_value[3];
  // map
  for(int i_mg=0; i_mg < n_nodes_mg; i_mg++) {
    int gl_node_bd = map_mg[i_mg];  // mg  node
    int node_med  = map_med[i_mg];  // med node

    fct->eval(node_med,n_cmp, bc_value);
    for(int jdim=first_cmp; jdim<first_cmp+n_cmp; jdim++) {
      const int kdof_top = mgsyst->_node_dof[Level][gl_node_bd+jdim*offset];
      old_sol_top->set(kdof_top , bc_value[jdim-first_cmp]);
    }
  }

  return;
}
#endif