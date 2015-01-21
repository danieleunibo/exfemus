#include "LibMeshUtil.h"
#include "LibMeshFunction.h"
// #include "libmesh/libmesh.h"
// #include "libmesh/point.h"
// #include "libmesh/elem.h"
#define FDEBUG 0
#include "MGFE_conf.h"
#include "LibMeshUtil.h"
#include "LibMeshFunction.h"
// #include "libmesh/libmesh.h"
// #include "libmesh/point.h"
// #include "libmesh/elem.h"
#define FDEBUG 0
// #include "Debug.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingUMesh.hxx"
#include <iomanip>
#include <limits>
#include "ParallelMeshExtended.h"





// using namespace ParaMEDMEM;
// using namespace libMesh;
// ======================================================
LibMeshFunction::~LibMeshFunction() {
  if (_field) _field->decrRef();
}

// =======================================================
// This function computes the mesh interface
// to set in the interface-function. The field in the
// interface-function is set by set_mesh_femus_interface
void LibMeshFunction::set_mesh_interface(
  const bool node_only,   //  
  const MGMesh * mesh,               // package-mesh
  const ParaMEDMEM::MEDCouplingUMesh * support   // med-mesh
) { // ==============================================================
  // Femus-mesh
  _mesh_femus = mesh;                                 // mesh
  int nlevels=_mesh_femus->_NoLevels-1;                  // top level
  int n_nodes=_mesh_femus->_NoNodes[nlevels];    // top level  n nodes
//   int n_elements_top=_mesh_femus->_NoElements[0][nlevels]; // top level  n elements
  int dim= _mesh_femus->_dim;                          // dimension//
   int n_elements=0;
   for(int ilev=0;ilev<nlevels+1;ilev++) n_elements +=_mesh_femus->_NoElements[0][ilev]; // top level  n elements
// int n_nod_el=(dim==DIMENSION)?NDOF_P:NDOF_PB;
  // Med-mesh
  _support_med = support;                               // mesh
  int dim_med=_support_med->getMeshDimension();         // dimension
  int n_elements_med =_support_med->getNumberOfCells();// n elements
  int n_nodes_med= _support_med->getNumberOfNodes();   //  n nodes
//   int n_nod_el_med=(dim_med==DIMENSION)?NDOF_FEM:NDOF_FEMB;
  
  std::set<const Node *, nodeCompare> s1, s2;
  std::vector<NodeExt> n1,n2;

  // Corresponding nodes -----------------------------------------------------
  _nodeID.clear();  // clea map node: FEMus-mesh -> MED-mesh
  CleanSetOfNodes(s2);
  // med-mesh   d=coords s1=map n=Node_vec
  const ParaMEDMEM::DataArrayDouble * d =
    _support_med->getCoords();// Med-mesh coordinates
  DataArrayToSetOfNodes(d, s1, n1);               // Med-mesh s1=map n1=Node_vec
  xyzToSetOfNodes(dim,n_nodes,mesh->_xyz, s2,n2); //femus-mesh s2=map n2=Node_vec
  commonPoints(_nodeID, s2, s1);                  // Common points 
  if(_nodeID.size()!=s1.size()) std::cout<< " not all nods were found ! \n";
  s1.clear();  n1.clear();// clean
  s2.clear();  n2.clear();// clean
  CleanSetOfNodes(s2);
 
//
  if (!node_only) {
    //  Corresponding element map ----------------------------------------------------------
    _elemID.clear();  // clea map node: FEMus-mesh -> MED-mesh
    std::set<const Node *, nodeCompare>  s1m2; CleanSetOfNodes(s1m2);
    std::vector<NodeExt> n1m2;
    double *xyz_m; xyz_m=new double [n_elements_med*dim_med];// coord center vect
    getSupportBarycenterCoords(dim_med,n_elements_med,d,xyz_m);// Med-mesh baricenter coords
    centerToSetOfNodes(dim_med,n_elements_med,xyz_m, s1m2,n1m2);   // Med-mesh s1=map n1=Node_
 n1m2.clear();
//  d->decrRef();

    std::set<const Node *, nodeCompare>  s2m2; /* CleanSetOfNodes(s2m2);*/
    std::vector<NodeExt> n2m2;
//      const int init_lev_el=_mesh_femus->_off_el[0][_mesh_femus->_NoLevels*0+nlevels];
    // Corresponding elements (through baricenter
    double *xyz1_m; xyz1_m=new double [n_elements*dim];       // MGMesh coord centers
    getMGMeshBarycenterCoords(dim,nlevels,n_nodes,n_elements,xyz1_m);               // MGMesh baricenter coords
    centerToSetOfNodes2(*_mesh_femus,dim,n_elements,nlevels,xyz1_m, s2m2,n2m2);// MGMesh baricenter nods
    n2m2.clear();
    commonPoints(_elemID, s2m2,s1m2);                  // common element ->_elemID
     if(_elemID.size()!=s1m2.size()) std::cout<< " not all elements were found ! \n";
    // clean
    delete []xyz1_m; delete []xyz_m;
//     s1.clear();  n1.clear();
    s2m2.clear(); /* n2m2.clear();*/// clean
    s1m2.clear(); /* n1m2.clear();*/// clean
     CleanSetOfNodes(s2);
  } else {




    // Corresponding faces
//     std::set<const NodeExt *, nodeCompare> s2;

//     const MeshBase::const_element_iterator end_el
//       = _mesh_femus->active_local_elements_end();
//     for (MeshBase::const_element_iterator el
//          = _mesh_femus->active_local_elements_begin() ;
//          el != end_el; ++el) {
//       const Elem* elem = *el;
//       for (unsigned int side =0 ; side < elem->n_sides(); side++)
//         if (elem->neighbor(side) == NULL) {
//           NodeExt * q = new NodeExt(elem->build_side(side)->centroid(),
//                                     elem->id(), side);
//           s2.insert(q);
//         };
  }
//
//     commonPoints(_faceID, s1, s2);
//     CleanSetOfNodes(s2);
//    }
//

//   s2m2.clear();  n2m2.clear();// clean
  return;
}

// ============================================================================
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
void LibMeshFunction::set_mesh_interface_nodeID(
  const ParallelMeshExtended * mesh,                      ///< Femus-mesh  (in)
  const ParaMEDMEM::MEDCouplingUMesh * support,   ///< med-mesh             (in)
  const int order_cmp                ///< only linear pt (1)        (in)
) { // ========================================================================
   // Femus-mesh
  _mesh_femus = mesh;                                 // mesh
  int nlevels=_mesh_femus->_NoLevels-1;                  // top level
  int n_nodes=_mesh_femus->_NoNodes[nlevels];    // top level  n nodes
//   int n_elements_top=_mesh_femus->_NoElements[0][nlevels]; // top level  n elements
  int dim= _mesh_femus->_dim;                          // dimension//
   int n_elements=0;
   for(int ilev=0;ilev<nlevels+1;ilev++) n_elements +=_mesh_femus->_NoElements[0][ilev]; // top level  n elements
// int n_nod_el=(dim==DIMENSION)?NDOF_P:NDOF_PB;
  // Med-mesh
  _support_med = support;                               // mesh
//   int dim_med=_support_med->getMeshDimension();         // dimension
//   int n_elements_med =_support_med->getNumberOfCells();// n elements
//   int n_nodes_med= _support_med->getNumberOfNodes();   //  n nodes
//   int n_nod_el_med=(dim_med==DIMENSION)?NDOF_FEM:NDOF_FEMB;
  
  std::set<const Node *, nodeCompare> s1, s2;
  std::vector<NodeExt> n1,n2;

  // Corresponding nodes -----------------------------------------------------
  _nodeID.clear();  // clea map node: FEMus-mesh -> MED-mesh
  CleanSetOfNodes(s2);
  // med-mesh   d=coords s1=map n=Node_vec
  const ParaMEDMEM::DataArrayDouble * d =
    _support_med->getCoords();// Med-mesh coordinates
  DataArrayToSetOfNodes(d, s1, n1);               // Med-mesh s1=map n1=Node_vec
//   xyzToSetOfNodes(dim,n_nodes,mesh->_xyz, s2,n2); //femus-mesh s2=map n2=Node_vec
  
  switch(order_cmp){
    case 2: 
    xyzToSetOfNodes_onBd(dim,n_nodes,mesh->_xyz,mesh->_bc_id,s2,n2);
    break;
    case 1: 
    xyzToSetOfNodes_onBd_lin(dim,n_nodes,mesh->_xyz,mesh->_bc_id,s2,n2);
    break;
    default:
      std::cout<< " \n LibMeshFunction::set_mesh_interface_nodeID:"<<
      " order not implemented \n";
      abort();
  }
//   xyzToSetOfNodes_onbd(dim,n_nodes,mesh->_xyz,mesh->_bc_id, s2,n2); //femus-mesh s2=map n2=Node_vec
  commonPoints(_nodeID, s2, s1);                  // Common points 
  if(_nodeID.size()!=s1.size()) std::cout<< " not all nods were found ! \n";
  s1.clear();  n1.clear();// clean
  s2.clear();  n2.clear();// clean
  CleanSetOfNodes(s2);
  
  
  

  return;
}












// ============================================================================
// This function computes the analytical expression over the mesh
// to set in the interface-function. The mesh1, mesh2 in the interface
// function are set by set_mesh_femus_interface
void LibMeshFunction::set_analytic_field_interface(
  const char *symbolic_eq,     // symbolic function
  int nComp             // number of componenents
) { // ========================================================================
  if (_field) _field->decrRef();
  //ON CELLS does not work (due to baryc in quad 9 or hex 27)
  ParaMEDMEM::TypeOfField type = ParaMEDMEM::ON_NODES;
  int dim=_support_med->getSpaceDimension();
  std::vector<std::string> vars(dim);
  if (dim > 0)    vars[0] = "x";
  if (dim > 1)    vars[1] = "y";
  if (dim > 2)    vars[2] = "z";

  _field = _support_med->fillFromAnalytic3(type, nComp, vars,symbolic_eq);
  _field->setName(symbolic_eq);
  _field->checkCoherency();
//   std::ostream out("prova.field");
  std::cout << "LibMeshFunction::set_analytic_field_interface \n";
 printOn(std::cout);
//   if (FDEBUG) {
//     int i,j;
//     for (i=0; i<_field->getNumberOfTuples(); i++) {
//       fDebugPos << " " << i << ": ";
//       for (j=0; j<nComp; j++)   fDebug << " " << _field->getIJ(i,j);
//       fDebug << std::endl;
//     }
//   }
  return;
}
// ==========================================================================
void LibMeshFunction::setField(
  const ParaMEDMEM::MEDCouplingFieldDouble *f
) {// =======================================================================
   if (_field) _field->decrRef();
//
  ParaMEDMEM::TypeOfField type = f->getTypeOfField();
  _field = ParaMEDMEM::MEDCouplingFieldDouble::New(type);
  *_field=*f;

  _field->setName(f->getName());
  _field->checkCoherency();
   printOn(std::cout);
  
 /* 
  
  if (_field) _field->decrRef();
   _field=static_cast<const ParaMEDMEM::MEDCouplingFieldDouble *> (f);
   _field->checkCoherency();
//   std::ostream out("prova.field");
  std::cout << "LibMeshFunction::setField \n";
 printOn(std::cout);*/
  
// //
//   ParaMEDMEM::TypeOfField type = f->getTypeOfField();
//   _field = ParaMEDMEM::MEDCouplingFieldDouble::New(type);
//   _field->setMesh(_support_med);
//   ParaMEDMEM::DataArrayDouble *array = ParaMEDMEM::DataArrayDouble::New();
//   int nTuples;
//   int nComp;
// 
//   const ParaMEDMEM::DataArrayDouble * arrayf = f->getArray(), *coord1, *coord2;
// //   if (FDEBUG) fDebug.array(arrayf, "arrayf");
// 
//   nTuples = _support_med->getNumberOfCells();
//   nComp = f->getNumberOfComponents();
//   array->alloc(nTuples, nComp);
// 
//   coord1 = _support_med->getBarycenterAndOwner();
//   coord2 = f->getMesh()->getBarycenterAndOwner();
// 
//   array->fillWithValue(std::numeric_limits<double>::max());
//   std::map<int, int> corresp;
//   commonPoints(corresp, coord1, coord2);
// 
//   for (std::map<int, int>::iterator it = corresp.begin();
//        it != corresp.end(); it++) {
//     for (int j=0; j<nComp; j++) {
// //       if (FDEBUG) fDebugPos << " " << it->first << " -> " << it->second
// // 		<< " val = " << arrayf->getIJ(it->first, j) << std::endl;
//       array->setIJ(it->second, j, arrayf->getIJ(it->first, j));
//     }
//   }
//   _field->setArray(array);
//   _field->setName(f->getName());
//   array->decrRef();
//   _field->checkCoherency();
 return;
}

// ==========================================================================
void LibMeshFunction::set_field_source(
  const ParaMEDMEM::MEDCouplingFieldDouble *f
) {// =======================================================================
  if (_field) _field->decrRef();
//
  ParaMEDMEM::TypeOfField type = f->getTypeOfField();
  _field = ParaMEDMEM::MEDCouplingFieldDouble::New(type);
  *_field=*f;
//   _field->setMesh(_support_med);
//   ParaMEDMEM::DataArrayDouble *array = ParaMEDMEM::DataArrayDouble::New();
//   int nTuples;
//   int nComp;
//
//   const ParaMEDMEM::DataArrayDouble * arrayf = f->getArray(), *coord1, *coord2;
// //   if (FDEBUG) fDebug.array(arrayf, "arrayf");
//
//   nTuples = _support_med->getNumberOfCells();
//   nComp = f->getNumberOfComponents();
//   array->alloc(nTuples, nComp);
//
//   coord1 = _support_med->getBarycenterAndOwner();
//   coord2 = f->getMesh()->getBarycenterAndOwner();
//
//   array->fillWithValue(std::numeric_limits<double>::max());
//   std::map<int, int> corresp;
//   commonPoints(corresp, coord1, coord2);
//
//   for (std::map<int, int>::iterator it = corresp.begin();
//        it != corresp.end(); it++) {
//     for (int j=0; j<nComp; j++) {
// //       if (FDEBUG) fDebugPos << " " << it->first << " -> " << it->second
// //              << " val = " << arrayf->getIJ(it->first, j) << std::endl;
//       array->setIJ(it->second, j, arrayf->getIJ(it->first, j));
//     }
//   }
//   _field->setArray(array);
  _field->setName(f->getName());
//   array->decrRef();
  _field->checkCoherency();
//   printOn(std::cout);
}

// ==========================================================================
ParaMEDMEM::MEDCouplingFieldDouble * LibMeshFunction::getField(
  char const* /*name */) {// =====================================================
//   if (FDEBUG) fDebugPos << "LibMeshFunction::getField _field = " << _field << std::endl;
//   if (!_field) return NULL;
//   if (FDEBUG) fDebugPos << "LibMeshFunction::getField _field = " << _field << std::endl;
//   _field->incrRef();
//   if (FDEBUG) fDebugPos << "LibMeshFunction::getField " << name << std::endl;
//   _field->setName(name);
  return _field;
}
// ==========================================================================
void LibMeshFunction::eval(
  int id,
  std::vector<double> & val
) {// =======================================================================
  if (_field == NULL) {for (int i=0; i<(int)val.size(); i++) val[i]=.0;}
  else {
    for (int i=0; i<(int)val.size(); i++) val[i]=_field->getIJ(_elemID[id],i);
  }
}

// ========================================================================
void LibMeshFunction::eval_elem(
  int iel,
  std::vector<double> & val
) { // ====================================================================
  if (_field == NULL) {for (int i=0; i<(int)val.size(); i++) val[i] =0.0;}
  else {
//     int mid = _faceID[std::pair<int, int>(id, side)];
    for (int i = 0; i<(int)val.size(); i++) val[i] = _field->getIJ(iel, i);
  }
}

// ========================================================================
void LibMeshFunction::eval_node(
  int node,
  std::vector<double> & val
) { // ====================================================================
  if (_field == NULL) {for (int i=0; i<(int)val.size(); i++) val[i] =0.0;}
  else {
//     int mid = _faceID[std::pair<int, int>(id, side)];
    for (int i = 0; i<(int)val.size(); i++) val[i] = _field->getIJ(node, i);
  }
}

// ========================================================================
void LibMeshFunction::eval(
  int id,
  int side,
  std::vector<double> & val
) { // ====================================================================
  if (_field == NULL) {for (int i=0; i<(int)val.size(); i++) val[i] =0.0;}
  else {
    int mid = _faceID[std::pair<int, int>(id, side)];
    for (int i = 0; i<(int)val.size(); i++) val[i] = _field->getIJ(mid, i);
  }
}
// =============================================================
// This function prints the interface function
void LibMeshFunction::printOn(
  std::ostream & out
) const { // ===================================================

  // set up title --------------------------------------------------------
  out << std::endl
//       << "=================================================="
      << std::endl
      << "LibMeshFunction printOn" << std::endl << std::endl;
  // ---------------------------------------------------------------------
  // ---------------- Mapping ---------------------------------------------
  //  node map MGMesh <-> MEDCoupling ------------------------------------
  out << "node correspondance (MGMesh <-> MEDCoupling)" << std::endl;
  std::map<int, int>::const_iterator itN = _nodeID.begin();
  for (; itN != _nodeID.end(); itN++)
    out << "  " << std::setw(4) << itN->first
        << " <-> " << std::setw(4) << itN->second << std::endl;
  out << std::endl;
  // element map MGMesh <-> MEDCoupling ----------------------------------
  out << "element correspondance (MGMesh <-> MEDCoupling)" << std::endl;
  std::map<int, int>::const_iterator itE = _elemID.begin();
  for (; itE != _elemID.end(); itE++)
    out << "  " << std::setw(4) << itE->first
        << " <-> " << std::setw(4) << itE->second << std::endl;
  out << std::endl;
  // face map    MGMesh <-> MEDCoupling ----------------------------------
  out << "face correspondance (MGMesh <-> MEDCoupling)" << std::endl;
  std::map<std::pair<int,int>, int>::const_iterator itF = _faceID.begin();
  for (; itF != _faceID.end(); itF++)
    out << " (" << std::setw(4) << itF->first.first
        << ", " << std::setw(4) << itF->first.second << ")"
        << " <-> " << std::setw(4) << itF->second << std::endl;
  out << std::endl;

  // ---------------------------------------------------------------------
  // ----------------  field ---------------------------------------------
  if (!_field) {
    out << "field= NULL" << std::endl;
    out << "==================================================" << std::endl;
    return;
  }
  ParaMEDMEM::TypeOfField type = _field->getTypeOfField();
  const ParaMEDMEM::DataArrayDouble * d;
  int n;
  const ParaMEDMEM::MEDCouplingUMesh * m
    = dynamic_cast<const ParaMEDMEM::MEDCouplingUMesh *>(_field->getMesh());

  out << "name : " << _field->getName() << std::endl;

  switch (type) {
  case ParaMEDMEM::ON_NODES:
    out << "type : ON_NODES" << std::endl;
    d = m->getCoords();
    d->incrRef();
    n = d->getNumberOfTuples();
    out << " " <<n << " nodes";
    break;
  case ParaMEDMEM::ON_CELLS:
    out << "type : ON_CELLS" << std::endl;
//     d = m->getBarycenterAndOwner();
    n = d->getNumberOfTuples();
    out << " " << n << " cells";
    break;
  default :
    std::cout<< "Error type"; abort();
  }

  ParaMEDMEM::DataArrayDouble * v = _field->getArray();
  out << " (" << v->getNumberOfTuples() << ")" << std::endl;
  int ncf=_field->getNumberOfComponents();
  int nc = d->getNumberOfComponents();
  for (int i=0; i<n; i++) {
    for (int j=0; j<nc; j++) out << (j == 0 ? "f( " : ", ")
                                   << std::setw(10) << std::fixed << d->getIJ(i,j);
    out << ") = ";
    for (int j=0; j<ncf; j++)  out << (j == 0 ? "( " : ", ")
                                     << std::setw(10) << std::fixed << v->getIJ(i,j);
    out << ")  ";
    out << std::endl;
  }
  out << std::endl << std::endl;
  out << "==================================================" << std::endl;

  //d->decrRef();
}
// ======================================================================
// This function get the barycenter for Med-element
void  LibMeshFunction::getSupportBarycenterCoords(
  int dim_med,
  int n_elements_med,
  const ParaMEDMEM::DataArrayDouble * d ,  // Med-mesh coords      <-
  double *xyz_m                            // Med-Mesh baricenters ->
) { // ==================================================================
  int n_nod_el_med=(dim_med==DIMENSION)?NDOF_FEM:NDOF_FEMB;
  std::vector<int> nodes1; //node of elems
  double x_m[DIMENSION];//tmp buffer
  for (int ielem=0; ielem<  n_elements_med; ielem++) {
    for (int idim=0; idim<DIMENSION; idim++) x_m[idim]=0.;
    _support_med->getNodeIdsOfCell(ielem, nodes1);
    for (int inode=0; inode< n_nod_el_med; inode++) {
      for (int idim=0; idim<dim_med; idim++) {
        double a=d->getIJ(nodes1[inode], idim);
        x_m[idim] += a;
      }  // idim
    } // end inode
    nodes1.clear();
    for (int idim=0; idim<dim_med; idim++){
      xyz_m[ielem+idim*n_elements_med]=x_m[idim]/n_nod_el_med;
//       std::cout << ielem+idim*n_elements_med<< " " << xyz_m[ielem+idim*n_elements_med]  <<"\n";
    }
  }
  return;
}

// ====================================================
void LibMeshFunction::getMGMeshBarycenterCoords(
  int dim,
  int nlevels,
  int ,
  int n_elements,
  double *xyz_m
) { // =================================================================
  int n_nod_el=(dim==DIMENSION)?NDOF_FEM:NDOF_FEMB;
  double x_m[DIMENSION];//tmp buffer
 int Level=nlevels;
  const int  offset   =_mesh_femus->_NoNodes[_mesh_femus->_NoLevels-1];
  //  2bB element interpolation over the fine mesh -----------------------

//   int ielemn=0;
  for (int iproc=0; iproc<_mesh_femus->_n_subdom; iproc++) {
    int nel_b=_mesh_femus->_off_el[0][_mesh_femus->_NoLevels*iproc+Level];
    int nel_e=_mesh_femus->_off_el[0][Level+_mesh_femus->_NoLevels*iproc+1];
    for (int iel=nel_b; iel <nel_e; iel++) {
      
      for (int idim=0; idim<dim; idim++) x_m[idim]=0.;
      for (int  inode=0; inode<n_nod_el; inode++)    {
        //get the global node number
        int el_conn = _mesh_femus->_el_map[0][iel*n_nod_el+inode];
        // get the element coordinades
        for (int  idim=0; idim<dim; idim++) {
          x_m[idim] += _mesh_femus->_xyz[el_conn+idim*offset];
        }
      }
    for (int idim=0; idim<dim; idim++){
      xyz_m[iel+idim*n_elements]=x_m[idim]/n_nod_el;
//       std::cout <<iel+idim*n_elements<< " " << xyz_m[iel+idim*n_elements]  <<"\n";
    }
//     ielemn++;
  } // ---- end iel -------------------------
} // 2bB end interpolation over the fine mesh -
return;
}
