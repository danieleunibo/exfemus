#define FDEBUG 0
#include <set>
#include <cstring>
// #include "libmesh/equation_systems.h"
// #include "libmesh/node.h"
// #include "libmesh/boundary_info.h"
#include "MGEquationsSystem.h"
#include "EquationSystemsExtended.h"
#include "Debug.hxx"
#include "LibMeshUtil.h"
#include "LibMeshFunction.h"
#include "numeric_vectorM.h"

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "Femus_BCTypes.h"
#include <MEDCouplingRefCountObject.hxx>
#include "MEDCouplingRemapper.hxx"

// using namespace libMesh;
using namespace ParaMEDMEM;

// ========================================================================
/// Constructor
EquationSystemsExtended::EquationSystemsExtended(
  MGUtils& mgutils_in,
  ParallelMeshExtended& mgmesh_in,
  MGFEMap& mgfemap_in,
  int nComp,
  int npoint_data,
  int ncell_data
) :
  MGEquationsSystem(mgutils_in,mgmesh_in,mgfemap_in,npoint_data,ncell_data),
// //    _mgutils(mgutils_in),  // file name class in
  _mesh(&mgmesh_in),    // mesh class in
//   *mg_utils[mesh],NULL,np_data,ncell_data)
  _nComp(nComp),
//  _source(NULL)
  _source(new LibMeshFunction) {
  int _nAverageSources=1;
  _AverageSources=new double[_nAverageSources];
  for(int i=0; i<_nAverageSources; i++) _AverageSources[i]=0.;
  _nComp=1;               // components

//   _comm(comm),           // communicator
//   problemName(pbName),   // problem name
//   _nComp(nComp),

//   _mesh= mgmesh_in;           // LIBMESH-mesh
//   _initLibMesh(NULL)
  // ==============================================================
//   _problemName=;
  // bc table -------------------------------------
//   _BCTypes["Mark"] = Mark;
//   _BCTypes["mark"] = Mark;
//   _BCTypes["UnMark"] = UnMark;
//   _BCTypes["Unmark"] = UnMark;
//   _BCTypes["unmark"] = UnMark;
//   _BCTypes["unMark"] = UnMark;
}

// ========================================================================
EquationSystemsExtended::~EquationSystemsExtended() {
  if(_source) delete _source;

  std::map<boundary_id_type, LibMeshFunction *>::iterator it = _interfaceFunMap.begin();
  for(; it != _interfaceFunMap.end(); it++) {
    if(it->second) delete it->second;
    it->second = NULL;
  }
}
// ========================================================================
LibMeshFunction * EquationSystemsExtended::getSource() {
  return _source;
}
// ========================================================================
// This function sets the source function (med-mesh,femus-mesh,field)
void EquationSystemsExtended::setSource(
  const MEDCouplingUMesh * medmesh,   // med-mesh
  const char *s                       // analytical expression
) {// ======================================================
//   MGMesh & mesh = get_mesh(); // femus-mesh

  _source->set_mesh_interface(false,&_mgmesh, medmesh); // set mesh in the source function
  _source->set_analytic_field_interface(s, _nComp);            // se
//   if (FDEBUG) fDebugPos << "setSource (string " << s << ")" << std::endl;
//   if (FDEBUG) _source->printOn(fDebug);
  return;
}
// ========================================================================
// This function sets the source function (med-mesh,femus-mesh,field)
void EquationSystemsExtended::setSource(
  const ParaMEDMEM::MEDCouplingUMesh * medmesh,
  const ParaMEDMEM::MEDCouplingFieldDouble * f
) {// ======================================================
//   MGMesh & mesh = get_mesh(); // femus-mesh

  _source->set_mesh_interface(1,&_mgmesh, medmesh); // set mesh in the source function
  _source->set_field_source(f);            // se
//   if (FDEBUG) fDebugPos << "setSource (string " << s << ")" << std::endl;
//   if (FDEBUG) _source->printOn(fDebug);
  return;
}
// ========================================================================
void EquationSystemsExtended::eraseSource() {
  if(_source) {
    delete _source;
    _source = NULL;
  }
}

// // =================================================
// void EquationSystemsExtended::setBoundaryConditionType(
//   int id,
//   const char * typeBC
// ) { // =============================================
//
//   // check
//   std::string s(typeBC);
//   if (_BCTypes.find(s) == _BCTypes.end()) {
//     std::cout<< "setBoundaryConditionType: bc type not found";
//     abort();
//   }
//   // set BCType
//   setBCType(id, _BCTypes[s]);
// }

// // ========================================================================
// void EquationSystemsExtended::setBCType(
//   boundary_id_type id,
//   FemusBCType type
// ) {
// //   if (FDEBUG) fDebugPos << "id " << id << " BCType : " << type << std::endl;
//   _bc_type[id] = type;
// }
// // =======================================================================
// FemusBCType EquationSystemsExtended::getBCType(boundary_id_type id) {
//   return (_bc_type.find(id) == _bc_type.end()) ? Undefined : _bc_type[id];
// }
// =======================================================================
void EquationSystemsExtended::setBC(
  boundary_id_type id,
  int  n_cmp,
  const MEDCouplingFieldDouble *field
) {// =======================================================================

  LibMeshFunction * fct = get_interface_function_BC(id);  // interface-function
  if(fct == NULL) return;
//
//   MGMesh & mesh = get_mesh();
//   std::map<std::pair<int, int>, int>::iterator itF = fct->FaceID().begin();
//   for (; itF != fct->FaceID().end(); itF++) {
//     int elem_id = itF->first.first;
//     int side = itF->first.second;
//     mesh.boundary_info->add_side(elem_id, side, id);
//   }
//
  fct->setField(field);
//   if (FDEBUG) {
//      fDebugPos << "setBC (type " << _bc_type[id]
// 	           << " field " << f->getName() << ")" << std::endl;
//      fct->printOn(fDebug);
//   }
  return;
}

void EquationSystemsExtended::setFieldBoundaryValues(
  int bc_id,
  int  n_cmp,
  const MEDCouplingFieldDouble *f
) {
  if(f->getTypeOfField() == ON_CELLS)
    setBC(bc_id,  n_cmp, f);
  else
    setNodeBC(bc_id,   n_cmp, f);

  return;
}

// ===========================================================
void EquationSystemsExtended::setNodeBC(
  boundary_id_type id,
  int  n_cmp,
  const MEDCouplingFieldDouble *f
) { //==========================================================

  LibMeshFunction * fct = get_interface_function_BC(id);
  if(fct == NULL) return;

//   if (FDEBUG) fDebug.array(f->getArray(), "ici f 3");
//
//   MEDCouplingFieldDouble * ff=MEDCouplingFieldDouble::New(ON_CELLS);
//   const MEDCouplingUMesh * mesh = fct->getSupport();
//   ff->setMesh(mesh);
//   int nCells = mesh->getNumberOfCells();
//   int nComp = f->getNumberOfComponents();
//   DataArrayDouble * d = DataArrayDouble::New();
//   const DataArrayDouble * df = f->getArray();
//   d->alloc(nCells, nComp);
//
//   const DataArrayInt * c = mesh->getNodalConnectivity();
//   const DataArrayInt * cI = mesh->getNodalConnectivityIndex();
//
//   for (int iCell=0; iCell < nCells; iCell++) {
//
//     INTERP_KERNEL::NormalizedCellType t
//       = mesh->getTypeOfCell(iCell);
//     int n_nodes_elem;
//
//     switch (t) {
//     case INTERP_KERNEL::NORM_SEG2:
//       n_nodes_elem = 2;
//       break;
//     case INTERP_KERNEL::NORM_SEG3:
//       n_nodes_elem = 3;
//       break;
//     case INTERP_KERNEL::NORM_TRI3:
//       n_nodes_elem = 3;
//       break;
//     case INTERP_KERNEL::NORM_TETRA4:
//       n_nodes_elem = 4;
//       break;
//     case INTERP_KERNEL::NORM_QUAD4:
//       n_nodes_elem = 4;
//       break;
//     case INTERP_KERNEL::NORM_QUAD8:
//       n_nodes_elem = 8;
//       std::cout << "\n I am quad8  \n";
//       break;
//     case INTERP_KERNEL::NORM_QUAD9:
//       n_nodes_elem = 9;
//       std::cout << "\n I am quad9  \n";
//       break;
//     case INTERP_KERNEL::NORM_TETRA10:
//       n_nodes_elem = 10;
//       std::cout << "\n I am tetra10  \n";
//       break;
//     case INTERP_KERNEL::NORM_HEXA8:
//       n_nodes_elem = 8;
//       std::cout << "\n I am hexa8  \n";
//       break;
//     case INTERP_KERNEL::NORM_HEXA20:
//       n_nodes_elem = 20;
//       std::cout << "\n I am hexa20  \n";
//       break;
//     case INTERP_KERNEL::NORM_HEXA27:
//       n_nodes_elem = 27;
//       std::cout << "\n I am hexa27  \n";
//       break;
//     default:
//       n_nodes_elem = 0;
//     }
//
//     int cstart = cI->getIJ(iCell, 0);
//     //   int cend = cI->getIJ(iCell+1, 0);
//
//     double v = 0.0;
//     for (int i=0; i<n_nodes_elem; i++) {
//       double x = df->getIJ(c->getIJ(cstart+1+i, 0), 0);
//       v += x;
//     }
//     v /= n_nodes_elem;
//     d->setIJ(iCell, 0, v);
//   }
//
//   ff->setArray(d);
//   ff->checkCoherency();

//   if (FDEBUG) fDebugPos << "setNodeBC (type " << _bc_type[id]
// 	    << " field " << ff->getName() << ")" << std::endl;
  setBC(id, n_cmp,f);
}

// =========================================================================
void EquationSystemsExtended::setBC(
  boundary_id_type id,
  int n_cmp,
  const char *s
) { // ======================================================================

  LibMeshFunction * fct = get_interface_function_BC(id);
  if(fct == NULL) return;
  fct->set_analytic_field_interface(s, n_cmp);
//
//   MeshBase & mesh = get_mesh();
//   std::map<std::pair<int, int>, int>::iterator itF = fct->FaceID().begin();
//   for (; itF != fct->FaceID().end(); itF++) {
//     int elem_id = itF->first.first;
//     int side = itF->first.second;
//     mesh.boundary_info->add_side(elem_id, side, id);
//   }
//   if (FDEBUG) fDebugPos << "setBC (type " << _bc_type[id]
// 	    << " string " << s <<")" << std::endl;
//   if (FDEBUG) fct->printOn(fDebug);
}


// =======================================================================
LibMeshFunction * EquationSystemsExtended::get_interface_function_BC(
  boundary_id_type id
) { // ===================================================================
  return _interfaceFunMap.find(id) == _interfaceFunMap.end() ? NULL : _interfaceFunMap[id];
}

// =======================================================================
void EquationSystemsExtended::erase_interface_function_BC(
  boundary_id_type id
) {// =========================================================================
  _interfaceFunMap.erase(id);
}
// ============================================================================
// This functon add an interface-function (med-mesh,femus-mesh, field)
//  in the interface-function boundary map (_interfaceFunMap):
//  the field in the interface-function is not assigned here
 int EquationSystemsExtended::addfunction_BC(       ///<map position(return)
    const boundary_id_type  bd_name_id,    ///< boundary id (int)  (in)
    const ParaMEDMEM::MEDCouplingUMesh * b, ///< med-submesh        (in)
    const int from_cmp,                     ///< initial id component
    const int n_cmp,                        ///< n components
    const int order_cmp,                    ///< order component (2=quad;lin=1)
    const bool on_nodes                     ///< values on nodes (true)                     ///< order component (2=quad;lin=1)
  ){// ========================================================================

  // the index in the map  _interfaceFunMap[i] (not used) 
//   int i=0;
//   for(; i<_interfaceFunMap.size() + 1; i++) {
//     if(_interfaceFunMap.find(bd_name_id) == _interfaceFunMap.end()) break;
//   }
  
  // new function
  LibMeshFunction *fun = new LibMeshFunction; 
  
  // set the mesh interface
  if(bd_name_id>9){// on bc nodes (bc_flag)
     if(on_nodes){
       fun->set_mesh_interface_nodeID(_mesh,b,order_cmp);
     }
     else{std::cout<< "EquationSystemsExtended::addfunction_BC: node_search_only=0 ? \n" ;abort();}
  }
 else{  // on elements (mat)
    fun->set_mesh_interface(false,_mesh, b);  // set the mesh and med-submesh
 }
 fun->printOn(std::cout);  // print
 // setting the fun in the map _interfaceFunMap at bd_name_id
 _interfaceFunMap[bd_name_id] = fun;                               // add to the map

  return bd_name_id; // local map positions
}

// ===========  LIBmesh ============================

// // =====================================================
// /// This function gives the identity to interfaces-meshes
// /// by adding a new interface-function
// boundary_id_type EquationSystemsExtended::defineIdInterface(
//   const boundary_id_type  bd_name_int, ///< global name id (int)
//   const MEDCouplingUMesh * b           ///< bc mesh
// 
// ) { // ========================================================
// 
//   /// A) call addfunction_BC(bd_name_int,b) // ********************************
//   boundary_id_type bc_id = addfunction_BC(bd_name_int,b);
//   //   mesh->boundary_info->build_node_list_from_side_list();
//   return bc_id;
// }
// =======================================================
void printField(const MEDCouplingFieldDouble *f) {
//   if (FDEBUG) {
//     fDebugPos << "field name " << f->getName() << std::endl;
//     printMesh(f->getMesh());
//     const DataArrayDouble * arrayf = f->getArray();
//     fDebug.array(arrayf, "values");
//   }
  return;
}
// ============================================================================
// This function gets the all the values on nodes on the boundary
// with identity "id" in the system "system_name"
MEDCouplingFieldDouble * EquationSystemsExtended::getValuesOnBoundary_nodes(
  int id,                           ///< int boundary identity   (in)
  const char *system_name           ///< system name             (in)
)  {// ========================================================================
  // from LibMesh function fct
  LibMeshFunction * fct = get_interface_function_BC(id);
  if(fct == NULL)  return NULL;
  int nNodes = fct->getSupport()->getNumberOfNodes();
  std::map<int, int> & NodeID = fct->NodeID();

  // from MGMesh
  int Level=_mesh->_NoLevels-1;  int offset=_mesh->_NoNodes[Level];
  // from  MGSystem* system;
  MGSolBase * eq_system=get_eqs(system_name);
  int nComp =eq_system->_n_vars;
  std::cerr << "\n EquationSystemsExtended nComp = " << nComp << std::endl;

  // A new LibMesh function f
  MEDCouplingFieldDouble * f = MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES);
  f->setMesh(fct->getSupport());
  f->setName(system_name);
  //   f->setNature(ParaMEDMEM::ConservativeVolumic);

  // array function to fill f
  DataArrayDouble *array = DataArrayDouble::New();
  array->alloc(nNodes,nComp);
  // filling array
  for(std::map<int, int>::iterator iN = NodeID.begin(); iN != NodeID.end(); iN++) {
    for(int iComp = 0; iComp < nComp; iComp++) {
      int node_med = iN->second;
      int node_mg = iN->first;
      double v = (*(eq_system->x_old[Level]))(node_mg+iComp*offset);
      array->setIJ(node_med,iComp,v);
    }
  }

  f->setArray(array);   // set array -> f
  array->decrRef();     // delete array
  f->checkCoherency();  // check f
//   printField(f);        // print f
//   std::ofstream fout("map.dat");
//   fct->printOn(fout);   // print fct
  return f;

}

// ============================================================================
// This function gets the  the value of the variable with id number
//  "variable_id" on nodes on the boundary with identity "id" in the
//  system "system_name"
MEDCouplingFieldDouble * EquationSystemsExtended::getValuesOnBoundary_nodes(
  int id,                            // int boundary identity   (in)
  const char *system_name,           // system name             (in)
  int         variable_id            // variable identity       (in)
)  {// ========================================================================
  // from LibMesh function fct
  LibMeshFunction * fct = get_interface_function_BC(id);
  if(fct == NULL)  return NULL;
  int nNodes = fct->getSupport()->getNumberOfNodes();
  std::map<int, int> & NodeID = fct->NodeID();

  // from MGMesh
  int Level=_mesh->_NoLevels-1;  int offset=_mesh->_NoNodes[Level];
  // from  MGSystem* system;
  MGSolBase * eq_system=get_eqs(system_name);
  int nComp =1;
  std::cerr << "\n EquationSystemsExtended nComp = " << nComp << std::endl;

  // A new LibMesh function f
  MEDCouplingFieldDouble * f = MEDCouplingFieldDouble::New(ParaMEDMEM::ON_NODES);
  f->setMesh(fct->getSupport());
  f->setName(system_name);
  //   f->setNature(ParaMEDMEM::ConservativeVolumic);

  // array function to fill f
  DataArrayDouble *array = DataArrayDouble::New();
  array->alloc(nNodes,nComp);
  // filling array
  for(std::map<int, int>::iterator iN = NodeID.begin(); iN != NodeID.end(); iN++) {
    int node_med = iN->second;
    int node_mg = iN->first;
    double v = (*(eq_system->x_old[Level]))(node_mg+variable_id *offset);
    array->setIJ(node_med,0,v);
  }
  f->setArray(array);   // set array -> f
  array->decrRef();     // delete array
  f->checkCoherency();  // check f
  
//     printField(f);        // print f
//    std::ofstream fout("map.dat");
//    fct->printOn(fout);   // print fct
  
  
  return f;

}

MEDCouplingFieldDouble * EquationSystemsExtended::getValuesOnBoundary(
  int id,
  std::vector<char *> vName
) {
//   LibMeshFunction * fct = equation_systems->getBC(id);
//   if (fct == NULL)  return NULL;
//   System* system; std::vector<int> vars_id;
//
//
//   if (!(system = _getVarIds(vName,vars_id)))  return NULL;
//
// //  system= &(equation_systems->get_system(0));
// //   vars_id.push_back(0);vars_id.push_back(1);
//
//   int nComp = vars_id.size();
//   std::cerr << "nComp = " << nComp << std::endl;
//
//   if (FDEBUG) fDebugPos << " _LibMeshProblem::getValuesOnBoundary "
//                           << std::endl;
//
//    MEDCouplingFieldDouble * f     = MEDCouplingFieldDouble::New(ON_CELLS);
//   f->setMesh(fct->getSupport());
//   f->setName(vName[0]);
//   f->setNature(ParaMEDMEM::ConservativeVolumic);
//
//   long nCells = fct->getSupport()->getNumberOfCells();
//   DataArrayDouble *array = DataArrayDouble::New();
//   array->alloc(nCells,nComp);
//
//   std::map<std::pair<int, int>, int> & faceID = fct->FaceID();
//
//   for (int iComp = 0; iComp < nComp; iComp++) {
//     int iVar = vars_id[iComp];
//     for (std::map<std::pair<int, int>, int>::iterator iF = faceID.begin(); iF != faceID.end(); iF++) {
//       int elemID = iF->first.first;
//       int side = iF->first.second;
//       int idF = iF->second;
//
//       double v = 0.0;
//       int nv = 0;
//
//       const Elem * E = mesh->elem(elemID);
//       for (int j=0; j<E->n_nodes(); j++)
//         if (E->is_node_on_side(j, side)) {
//           nv ++;
//           v += system->current_solution(E->get_node(j)->dof_number(0,iVar,0));
//         }
//       if (nv > 0) v /= nv;
//       array->setIJ(idF,iComp,v);
//     }
//   }
//
//   f->setArray(array);
//   f->setName(vName[0]);
//   array->decrRef();
//   f->checkCoherency();
//   printField(f);
//   return f;
  return NULL;
}
// // =================================================================
// MEDCouplingFieldDouble * EquationSystemsExtended::getValuesOnBoundary_nodes(
//   int id,
//   std::vector<char *> vName
// ) const {
//
//   LibMeshFunction * fct = equation_systems->getBC(id);
//   if (fct == NULL)  return NULL;
//   System* system; std::vector<int> vars_id;
//
//
//   if (!(system = _getVarIds(vName,vars_id)))  return NULL;
//
// //  system= &(equation_systems->get_system(0));
// //   vars_id.push_back(0);vars_id.push_back(1);
//
//   int nComp = vars_id.size();
//   std::cerr << "nComp = " << nComp << std::endl;
//
//   if (FDEBUG) fDebugPos << " _LibMeshProblem::getValuesOnBoundary "
//                           << std::endl;
//
//   MEDCouplingFieldDouble * f  = MEDCouplingFieldDouble::New(ON_NODES);
//   f->setMesh(fct->getSupport());
//   f->setName(vName[0]);
// //   f->setNature(ParaMEDMEM::ConservativeVolumic);
//
//   long nNodes = (fct->getSupport())->getNumberOfNodes();
//   DataArrayDouble *array = DataArrayDouble::New();
//   array->alloc(nNodes,nComp);
//
//
//
//   const std::vector<int> & id_share = mesh->getNodeID();
//
//   mesh->boundary_info->print_info(std::cerr);
//
//
//   MeshBase::const_element_iterator el = mesh->active_local_elements_begin();
//   const MeshBase::const_element_iterator end_el = mesh->active_local_elements_end();
//
//   for (; el != end_el; ++el) {
//     Elem const* elem = *el;
//     for (int iComp = 0; iComp < nComp; iComp++) {
//       int iVar = vars_id[iComp];
//
//       if (elem->on_boundary()) {
//         for (unsigned int s=0; s<elem->n_sides(); s++)
//           if (elem->neighbor(s) == NULL)
//             if (mesh->boundary_info->has_boundary_id(elem, s, id))
//               for (int j=0; j<elem->n_nodes(); j++) {
//                 if (elem->is_node_on_side(j, s)) {
//                   const unsigned int id_nr = elem->get_node(j)->id();
//                   const unsigned int dof_nr = elem->get_node(j)->dof_number(0,iVar,0);
//                   array->setIJ(fct->NodeID()[id_nr],iComp,system->current_solution(dof_nr));
//                   { std::cerr << "current " << fct->NodeID()[id_nr] << " " << dof_nr << "\n"; }
//                 }
//               }
//       }
//
// //       unsigned int k = 0;
// //       for(k=0; (k<nNodes) && (k < id_share.size());++k) {
// // //         if (FDEBUG) { fDebugPos << "k " << k << " " << std::endl; }
// //         int iN = id_share[k];
// //         if (FDEBUG) { fDebugPos << "iN " << iN << " " << std::endl; }
// //         const Node & N = mesh->node(iN);
// //         if( mesh->boundary_info->has_boundary_id(&N, id))
// //         {
// //           const unsigned int dof_nr = N.dof_number(0,iVar,0);
// //           if (FDEBUG) { fDebugPos << "current " << k << " "
// //                                 << system->current_solution(dof_nr) << "\n"; }
// //           array->setIJ(fct->NodeID()[dof_nr],iComp,system->current_solution(dof_nr));
// //         }
// //       }
// //       for (; k<nNodes; ++k) array->setIJ(k,iComp,0.0);
//     }
//   }
//
//
//
//
// //   LibMeshFunction * fct = equation_systems->getBC(id);
// //   if (fct == NULL)  return NULL;
// //   System* system; std::vector<int> vars_id;
// //
// //
// //   if (!(system = _getVarIds(vName,vars_id)))  return NULL;
// //
// // //  system= &(equation_systems->get_system(0));
// // //   vars_id.push_back(0);vars_id.push_back(1);
// //
// //   int nComp = vars_id.size();
// //   std::cerr << "nComp = " << nComp << std::endl;
// //
// //   if (FDEBUG) fDebugPos << " _LibMeshProblem::getValuesOnBoundary "
// //                         << std::endl;
// //
// //   MEDCouplingFieldDouble * f  = MEDCouplingFieldDouble::New(ON_NODES);
// //   f->setMesh(fct->getSupport());
// //   f->setName(vName[0]);
// // //   f->setNature(ParaMEDMEM::ConservativeVolumic);
// //
// //   long nNodes = (fct->getSupport())->getNumberOfNodes();
// //   DataArrayDouble *array = DataArrayDouble::New();
// //   array->alloc(nNodes,nComp);
// //
// //
// //   std::map<std::pair<int, int>, int> & faceID = fct->FaceID();
// //   std::map<int, int> & nodeID = fct->NodeID();
// //
// //   for (int iComp = 0; iComp < nComp; iComp++) {
// //     int iVar = vars_id[iComp];
// //     for(std::map<std::pair<int, int>, int>::iterator iF = faceID.begin();iF != faceID.end(); iF++) {
// //       int elemID = iF->first.first;
// //       int side = iF->first.second;
// //       int idF = iF->second;
// //
// //       double v = 0.;
// //
// //       const Elem * E = mesh->elem(elemID);
// //       for (int j=0; j<E->n_nodes(); j++)
// //         if (E->is_node_on_side(j, side)) {
// //           v = 100.*( system->current_solution(E->get_node(j)->dof_number(0,iVar,0)));
// //           int id_node = (E->get_node(j))->id();
// //           array->setIJ(nodeID[id_node],iComp,v);
// //         }
// //     }
// //   }
//
//
//
//
//   f->setArray(array);
//   f->setName(vName[0]);
//   array->decrRef();
//   f->checkCoherency();
//   printField(f);
//   return f;
// }

// void _initSystem(){
//  set_mesh_femus(ParallelMeshExtended & mg_mesh_femus_in);
//   set_mesh_med(MEDCouplingUMesh  & mg_mesh_med_in);
//
// };


// // =====================================================================
void EquationSystemsExtended::_getVarIds(
  const char *name,
  std::vector<int> & varIds
) const { // ============================================================


  varIds.resize(1);
  varIds[0]=0;
//   libMesh::System * system = NULL;
//   int i_system, n_systems = equation_systems->n_systems();
//
//
//    std::string datavar(name);
//   std::string dataFile(name);
//   std::string prefix;
//    size_t p = dataFile.find_last_of("*");
//    size_t q = dataFile.find_last_of("V");
//    p = (p == std::string::npos) ? 0 : p+1;
//   prefix = dataFile.substr(p, q-p);
//
//
//
//    int lv = std::strlen(prefix.c_str());
//  int lv = std::strlen(name);
//   varIds.resize(0);

//   for (i_system=0; i_system<n_systems; i_system++) {
//     libMesh::System * s = &(equation_systems->get_system(i_system));
//     int j, n_vars = s->n_vars();
//     for (j=0; j<n_vars; j++) {
//       int lv_lib = strlen((s->variable_name(j)).c_str());
// //       if (s->variable_name(j).compare(0, lv, name) == 0)
//       std::string lib_str=s->variable_name(j);
//       std::size_t found = datavar.find(lib_str);
//       if (found<lv)
//         varIds.push_back(j);
//     }
//     if (varIds.size() > 0) {
//       system = s;
//       break;
//     }
//   }
//
//   return system;
}
//
// // =====================================================================
void EquationSystemsExtended::_getVarIds(
  std::vector< char *>name,
  std::vector<int> & varIds
) const { // ============================================================
//   MGEquationsSystem * system = NULL;
//   int i_system, n_systems = equation_systems->n_systems();
//
//   varIds.resize(0);
//
//   for (i_system=0; i_system<n_systems; i_system++) {
//     libMesh::System * s = &(equation_systems->get_system(i_system));
//     int j, n_vars = s->n_vars();
//     for (j=0; j<n_vars; j++) {
//       for (int k=0; k<name.size(); k++) {
//         int lv = strlen(name[k]);
//         if (s->variable_name(j).compare(0, lv, name[k]) == 0)
//           varIds.push_back(j);
//       }
//     }
//     if (varIds.size() > 0) {
//       system = s;
//       break;
//     }
//   }
//
//    return system;
  varIds.resize(1);
  varIds[0]=0;
  return;
}
