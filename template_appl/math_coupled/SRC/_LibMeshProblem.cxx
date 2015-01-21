#define FDEBUG 0

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"

#include "libmesh/fe.h"

#include "libmesh/quadrature_gauss.h"

#include "libmesh/dof_map.h"
#include "libmesh/point.h"


#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

#include "libmesh/perf_log.h"

#include "libmesh/elem.h"

#include "libmesh/string_to_enum.h"

#include "ParallelMeshExtended.hxx"
#include "_LibMeshProblem.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "InterpKernelExprParser.hxx"
#include "LibMeshFunction.hxx"
#include "EquationSystemsExtended.hxx"
#include "MEDCouplingNatureOfField.hxx"

#include "Debug.hxx"
#include "libmesh/boundary_info.h"

using namespace ParaMEDMEM;
using namespace libMesh;

// =============================================================================
_LibMeshProblem::_LibMeshProblem(
  const char *pbName,   // name
  int nComp,            // number of component
  MPI_Comm comm         // communicator
) :
  _comm(comm),           // communicator 
  problemName(pbName),   // problem name
  _nComp(nComp),         // components  
  med_mesh(NULL),        // MED-mesh
  mesh(NULL),            // LIBMESH-mesh
  equation_systems(NULL),// equation system
  _initLibMesh(NULL) {
  // ==============================================================

  int argc = 0;
  char ** argv = new char *[1];
  argv[0] = strdup("LibMesh_Salome_Component");
//   _initLibMesh = new LibMeshInit(argc, argv);
//     ,      (comm != MPI_COMM_NULL ? comm : MPI_COMM_WORLD));

//   free(argv[0]);  delete [] argv;
  // bc table -------------------------------------
  _BCTypes["Dirichlet"] = Dirichlet;
  _BCTypes["dirichlet"] = Dirichlet;
  _BCTypes["Neumann"] = Neumann;
  _BCTypes["neumann"] = Neumann;
  _BCTypes["Displacement"] = Displacement;
  _BCTypes["displacement"] = Displacement;
  _BCTypes["Force"] = Force;
  _BCTypes["force"] = Force;
  // ------------------------------------------------
  if (FDEBUG) { fDebugPos << "ok ici  _BCTypes " << _BCTypes.size() << std::endl; }
  return;
}

_LibMeshProblem::~_LibMeshProblem() {
  delete _initLibMesh;
}

// // ========================================================
// void _LibMeshProblem::setMesh(
//   const MEDCouplingUMesh * m           // MED-mesh
// ) { // =====================================================
//   med_mesh = m;
//   mesh = new ParallelMeshExtended(*m, _comm);
//   equation_systems = new EquationSystemsExtended(*mesh, _nComp);
//   _initSystem();
// }
// 
// // ===============================================================
// void _LibMeshProblem::solve() {
//   // =============================================================
//   LinearImplicitSystem& system= equation_systems->get_system<LinearImplicitSystem>(problemName.c_str());
//   system.solve();
// 
//   std::vector <void (*)(libMesh::EquationSystems&)>::iterator it;
// 
//   for (it = _post_processing_functions.begin(); it != _post_processing_functions.end(); it++) {
//     (*it)(*equation_systems);
//   }
// }
// 
// // =====================================================================
// libMesh::System * _LibMeshProblem::_getVarIds(
//   const char *name,
//   std::vector<int> & varIds
// ) const { // ============================================================
//   libMesh::System * system = NULL;
//   int i_system, n_systems = equation_systems->n_systems();
// 
// 
//   std::string datavar(name);
// //   std::string dataFile(name);
// //   std::string prefix;
// //   size_t p = dataFile.find_last_of("*");
// //   size_t q = dataFile.find_last_of("V");
// //   p = (p == std::string::npos) ? 0 : p+1;
// //   prefix = dataFile.substr(p, q-p);
// //
// //
// //
// //   int lv = strlen(prefix.c_str());
//   int lv = strlen(name);
//   varIds.resize(0);
// 
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
// }
// 
// // =====================================================================
// libMesh::System * _LibMeshProblem::_getVarIds(
//   std::vector< char *>name,
//   std::vector<int> & varIds
// ) const { // ============================================================
//   libMesh::System * system = NULL;
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
//   return system;
// }
// 
// 
// // ============================================================
// MEDCouplingFieldDouble * _LibMeshProblem::_getOutputField(
//   const char *vName,
//   TypeOfField typeF
// ) const {
//   System* system;
//   std::vector<int> vars_id;
// 
//   if (!(system = _getVarIds(vName,vars_id))) return NULL;
// 
//   int nComp = vars_id.size();
// 
//   MEDCouplingFieldDouble * f   = MEDCouplingFieldDouble::New(typeF);
//   f->setMesh(med_mesh);
//   f->setName(vName);
//   DataArrayDouble *array = DataArrayDouble::New();
// 
//   if (typeF == ON_NODES) {
//     long n = med_mesh->getNumberOfNodes();
// 
//     if (FDEBUG) { fDebugPos << "n = " << n << std::endl; }
//     array->alloc(n, nComp);
// 
//     const std::vector<int> & id = mesh->getNodeID();
// 
//     for (int iComp = 0; iComp < nComp; iComp++) {
//       int iVar = vars_id[iComp];
//       unsigned int k = 0;
//       for (k=0; (k<n) && (k < id.size()); ++k) {
//         if (FDEBUG) { fDebugPos << "k " << k << " " << std::endl; }
//         int iN = id[k];
//         if (FDEBUG) { fDebugPos << "iN " << iN << " " << std::endl; }
//         const Node & N = mesh->node(iN);
//         const unsigned int dof_nr = N.dof_number(0,iVar,0);
//         if (FDEBUG) {
//           fDebugPos << "current " << k << " "
//                     << system->current_solution(dof_nr) << "\n";
//         }
//         array->setIJ(k,iComp,system->current_solution(dof_nr));
//       }
//       for (; k<n; ++k) array->setIJ(k,iComp,0.0);
//     }
//   } else { // typeF == ON_CELLS
// 
//     long n = med_mesh->getNumberOfCells();
//     array->alloc(n, nComp);
//     const DofMap& dof_map = system->get_dof_map();
//     std::vector<dof_id_type> dof_indices_var;
//     const std::vector<int> & id = mesh->getElemID();
// 
//     for (int iComp = 0; iComp < nComp; iComp++) {
//       int iVar = vars_id[iComp];
// 
//       unsigned int k = 0;
//       for (k=0; (k<n) && (k < id.size()); ++k) {
//         int iE = id[k];
//         const Elem * E = mesh->elem(iE);
//         dof_map.dof_indices(E, dof_indices_var, iVar);
//         dof_id_type dof_index = dof_indices_var[0];
//         if ((system->solution->first_local_index() <= dof_index) &&
//             (dof_index < system->solution->last_local_index()))
//           array->setIJ(k,iComp, (*system->solution)(dof_index));
//       }
// 
//       for (; k<n; ++k) array->setIJ(k,iComp,0.0);
//     }
//   }
// 
// 
//   if (FDEBUG) fDebug.flush();
//   f->setArray(array);
//   array->decrRef();
//   f->checkCoherency();
//   if (FDEBUG) { fDebugPos << "ok  ici" << std::endl; }
// 
//   return f;
// }
// 
// // ========================================
// void _LibMeshProblem::terminate(
// ) { // ======================================
//   delete equation_systems;
//   delete mesh;
// }
// 
// // =================================================
// void _LibMeshProblem::setBoundaryConditionType(
//   int id,
//   const char * typeBC
// ) { // =============================================
// 
//   std::string s(typeBC);
//   if (_BCTypes.find(s) == _BCTypes.end())
//     return;
// 
//   equation_systems->setBCType(id, _BCTypes[s]);
// }
// 
// 
// 
// void _LibMeshProblem::setAnalyticBoundaryValues(int bc_id, const char *s) {
//   equation_systems->setBC(bc_id, s);
// }
// 
// void _LibMeshProblem::setFieldBoundaryValues
// (int bc_id, const MEDCouplingFieldDouble *f) {
//   if (f->getTypeOfField() == ON_CELLS)
//     equation_systems->setBC(bc_id, f);
//   else
//     equation_systems->setNodeBC(bc_id, f);
// }
// 
// // *******************************************************************
// //  **********************  sources **********************************
// 
// // ===================================================================
// // This function set the anaytic expression for source
// void _LibMeshProblem::setAnalyticSource(
//   const char *code
// ) {// ================================================================
//   equation_systems->setSource(med_mesh, code);
// }
// 
// // ===================================================================
// // This function set the  MEDCouplingFieldDouble *f into the source
// void _LibMeshProblem::setSource(
//   const MEDCouplingFieldDouble *f         // numerical source
// ) { // ===============================================================
//   f;
// }
// 
// // ===================================================================
// // This function set the average source vector
// void _LibMeshProblem::setAverageSources(
//   int k,
//   double val[]
// ) { // ===============================================================
//   equation_systems->setAverageSources(k,val);
//   return;
// }
// 
// 
// 
// /*
// void _LibMeshProblem::setNAverageSources(int n_csources) {
//   equation_systems->setNAverageSources(n_csources);
//   return;
// }
// */
// 
// boundary_id_type _LibMeshProblem::defineBoundary(
//   const MEDCouplingUMesh * b
// ) {
//   boundary_id_type bc_id = equation_systems->addBC(b);
//   mesh->boundary_info->build_node_list_from_side_list();
//   return bc_id;
// }
// 
// void printMesh(const MEDCouplingMesh *m) {
//   if (FDEBUG) {
//     fDebugPos << "mesh name " << m->getName() << std::endl;
//     fDebugPos << " nodes : " << m->getNumberOfNodes() << std::endl;
//     fDebugPos << " cells : " << m->getNumberOfCells() << std::endl;
//   }
// }
// 
// void printField(const MEDCouplingFieldDouble *f) {
//   if (FDEBUG) {
//     fDebugPos << "field name " << f->getName() << std::endl;
//     printMesh(f->getMesh());
//     const DataArrayDouble * arrayf = f->getArray();
//     fDebug.array(arrayf, "values");
//   }
// }
// 
// 
// // ==================================================================
// MEDCouplingFieldDouble * _LibMeshProblem::getValuesOnBoundary(
//   int id,
//   const char *vName
// ) const {
//   LibMeshFunction * fct = equation_systems->getBC(id);
//   if (fct == NULL)  return NULL;
// 
//   System* system;
//   std::vector<int> vars_id;
// 
//   if (!(system = _getVarIds(vName,vars_id)))    return NULL;
// 
//   int nComp = vars_id.size();
//   std::cerr << "nComp = " << nComp << std::endl;
// 
//   if (FDEBUG) fDebugPos << " _LibMeshProblem::getValuesOnBoundary "
//                           << std::endl;
// 
//   MEDCouplingFieldDouble * f = MEDCouplingFieldDouble::New(ON_CELLS);
//   f->setMesh(fct->getSupport());
//   f->setName(vName);
//   //   f->setNature(ParaMEDMEM::ConservativeVolumic);
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
//       const Elem * E = mesh->elem(elemID);  // libmesh element
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
//   f->setName(vName);
//   array->decrRef();
//   f->checkCoherency();
//   printField(f);
//   std::ofstream fout("map.dat");
//   fct->printOn(fout);
//   return f;
// }
// 
// MEDCouplingFieldDouble * _LibMeshProblem::getValuesOnBoundary
// (int id,  std::vector<char *> vName) const {
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
//   MEDCouplingFieldDouble * f
//     = MEDCouplingFieldDouble::New(ON_CELLS);
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
// }
// 
// MEDCouplingFieldDouble * _LibMeshProblem::getValuesOnBoundary_nodes
// (int id,  std::vector<char *> vName) const {
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
