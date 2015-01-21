#include "LibMeshUtil.h"
#include "MEDCouplingMemArray.hxx"
#include "MGMesh.h"
// #include "Debug.hxx"

using namespace ParaMEDMEM;

const double nodeCompare::eps = 1e-16;
// ================================================================================
// This function computes the common point map:  extNode(2mesh) -> Node(1mesh)
void commonPoints(
  std::map<int, int> & corresp,                   // correspondence map  ->
  const std::set<const Node *, nodeCompare> & p,  // Nodes from 1 mesh   <-
  const std::set<const Node *, nodeCompare> & q   // Nodes form 2 mesh   <-
){// ==============================================================================
  // map 1mesh -> 2 mesh
  std::set<const Node *>::const_iterator pIt, qIt;
  // Find the nodes
  int icount=0; int jcount=0;
  for (qIt = q.begin(); qIt != q.end(); ++qIt) {
//      int t1=(*qIt)->id();
//      std::cout << icount <<" " << t1 << "\n ";
//      icount++;
//     const Node * nodeq=qIt.first;
    if ((pIt = p.find(*qIt)) != p.end()) {
//       const Node * nodep=pIt;
     
//       int t2=(*pIt)->id();
//       int r=0;
//       corresp[(*qIt)->_id] = (*pIt)->_id;
       corresp[ (*pIt)->_id] =(*qIt)->_id;
    }
   else{
      int a=(*qIt)->id();
//       std::cout << icount <<" " << t1 << "\n ";
//        std::cout << icount <<" " << t1 << "\n ";
      std::cout<< "commonPoints : point not found "<<a; abort();
     
     
  }
  } 
  return;
}
// =================================================================================
// This function computes the  common point map: extNode(2mesh) -> Node(1mesh)
void commonPoints(
  std::map<std::pair<int, int>, int> & corresp,     // correspondence map  ->
  const std::set<const Node *, nodeCompare> & p,    // Nodes from 1 mesh   <-
  const std::set<const NodeExt *, nodeCompare> & q  // Nodes form 2 mesh   <-
){ // ==============================================================================
  std::set<const Node *>::const_iterator pIt;
  std::set<const NodeExt *>::const_iterator qIt;

  for (qIt = q.begin(); qIt != q.end(); ++qIt) {
    if ((pIt = p.find(*qIt)) != p.end()) {
      const NodeExt * q = *qIt;
      corresp[std::pair<int, int>(q->id(), q->id2())] = (*pIt)->id();
    }
  }
  return;
} 

// ===========================================================================
// This function computes the  common point map: Node(2mesh) -> Node(1mesh)
void commonPoints(
  std::map<int, int> & corresp,
  const ParaMEDMEM::DataArrayDouble * d1, 
  const ParaMEDMEM::DataArrayDouble * d2)
{ // ==========================================================================
  std::vector<NodeExt> n1, n2;
  std::set<const Node *, nodeCompare> s1, s2;
  
  std::cout << "d1 = " << d1 << " d2 = " << d2 << std::endl;
  DataArrayToSetOfNodes(d1, s1, n1);
  DataArrayToSetOfNodes(d2, s2, n2);
  commonPoints(corresp, s1, s2);
  return;
}
// =======================
// This function find the node id in the map
const Node * findNode(
  const std::set<const Node *, nodeCompare> & p, 
  int id
){
  for (std::set<const Node *>::const_iterator pIt = p.begin();   pIt != p.end(); ++pIt)
    if ((*pIt)->id() == id) {
      return *pIt;
    }
  return NULL;
}
// =======================================================
// This function set the med-mesh point into ExtNode victor n
void DataArrayToSetOfNodes(
  const DataArrayDouble * d,                 // coords   <-
  std::set<const Node *, nodeCompare> & s,   // (Node map, compare) ->
  std::vector<NodeExt> & extNode_vect               // extNode vector   -> 
){  // ==================================================
  
  // From med-mesh coordinates
  int ndim = d->getNumberOfTuples();
  int dim = d->getNumberOfComponents();
  // setting node ext
  extNode_vect.clear(); if (ndim == 0)  return;  extNode_vect.resize(ndim);
  
//   unsigned int i, idim;
  for (int i=0; i<ndim; i++) {
    // set  coordinates to Node ext
    for (int idim=0; idim<dim; idim++)         extNode_vect[i](idim) = d->getIJ(i, idim);
    for (int idim=dim; idim<DIMENSION; idim++) extNode_vect[i](idim) = 0.0;
    extNode_vect[i].set_id() = i;
    //     n[i].processor_id(345);
    // insert ext Node to s structure
    s.insert(&(extNode_vect[i]));
  }
  return;
}


// ============================================================================
/// This function set the boundary xyz coord vector into ExtNode vector n
void xyzToSetOfNodes_onBd(
  const int dim,                           ///< dimension                 (in) 
  const int n_nodes,                       ///< number of nodes           (in)
  const double  xyz[],                     ///< coord vector              (in)
  const std::vector<int> bc_flag,          ///< boundary flag >9          (in)
//   int n_bd_nodes,                    ///< number of  bd nodes found (out)
  std::set<const Node *, nodeCompare> & s, ///< (Node map, compare)       (out)
  std::vector<NodeExt> & extNode_vect      ///< extNode vector            (out) 
){  // ==================================================
   
  // From med-mesh coordinates
//   int ndim = d->getNumberOfTuples();
//   int dim = d->getNumberOfComponents();
  // setting node ext
  extNode_vect.clear(); if (n_nodes == 0)  return;  extNode_vect.resize(n_nodes);
//     n_bd_nodes=0;
//   unsigned int i, idim;
  for (int i=0; i<n_nodes; i++) {
    // set  coordinates to Node ext 
    if(fabs(bc_flag[i]) >9){
//       n_bd_nodes++;
    for (int idim=0; idim<dim; idim++)         extNode_vect[i](idim) = xyz[i+idim*n_nodes];// d->getIJ(i, idim);
    for (int idim=dim; idim<DIMENSION; idim++) extNode_vect[i](idim) = 0.0;
    extNode_vect[i].set_id() = i;
    //     n[i].processor_id(345);
    // insert ext Node to s structure
    s.insert(&(extNode_vect[i]));
    }
  }
  return;
}
 
// ============================================================================
/// This function set the boundary xyz coord vector into ExtNode vector n
void xyzToSetOfNodes_onBd_lin(
  const int dim,                           ///< dimension                 (in) 
  const int n_nodes,                       ///< number of nodes           (in)
  const double  xyz[],                     ///< coord vector              (in)
  const std::vector<int> bc_flag,          ///< boundary flag >9          (in)
//   int n_bd_nodes,                    ///< number of  bd nodes found (out)
  std::set<const Node *, nodeCompare> & s, ///< (Node map, compare)       (out)
  std::vector<NodeExt> & extNode_vect      ///< extNode vector            (out) 
){  // ==================================================
   
  // From med-mesh coordinates
//   int ndim = d->getNumberOfTuples();
//   int dim = d->getNumberOfComponents();
  // setting node ext
  extNode_vect.clear(); if (n_nodes == 0)  return;  extNode_vect.resize(n_nodes);
//     n_bd_nodes=0;
//   unsigned int i, idim;
  for (int i=0; i<n_nodes; i++) {
    // set  coordinates to Node ext 
    if(bc_flag[i] >9){
//       n_bd_nodes++;
    for (int idim=0; idim<dim; idim++)         extNode_vect[i](idim) = xyz[i+idim*n_nodes];// d->getIJ(i, idim);
    for (int idim=dim; idim<DIMENSION; idim++) extNode_vect[i](idim) = 0.0;
    extNode_vect[i].set_id() = i;
    //     n[i].processor_id(345);
    // insert ext Node to s structure
    s.insert(&(extNode_vect[i]));
    }
  }
  return;
}



// =======================================================
// This function set the xyz vector into ExtNode victor n
void xyzToSetOfNodes(
  const int dim,
  const int n_nodes,
  const double  xyz[],                 // coords   <-
  std::set<const Node *, nodeCompare> & s,   // (Node map, compare) ->
  std::vector<NodeExt> & extNode_vect               // extNode vector   -> 
){  // ==================================================
  
  // From med-mesh coordinates
//   int ndim = d->getNumberOfTuples();
//   int dim = d->getNumberOfComponents();
  // setting node ext
  extNode_vect.clear(); if (n_nodes == 0)  return;  extNode_vect.resize(n_nodes);
  
//   unsigned int i, idim;
  for (int i=0; i<n_nodes; i++) {
    // set  coordinates to Node ext
    for (int idim=0; idim<dim; idim++)         extNode_vect[i](idim) = xyz[i+idim*n_nodes];// d->getIJ(i, idim);
    for (int idim=dim; idim<DIMENSION; idim++) extNode_vect[i](idim) = 0.0;
    extNode_vect[i].set_id() = i;
    //     n[i].processor_id(345);
    // insert ext Node to s structure
    s.insert(&(extNode_vect[i]));
  }
  return;
}


  
// =======================================================
// This function set the xyz vector into ExtNode victor n
void centerToSetOfNodes(
  const int dim_med,
  const int n_element_med,
  const double  xyz_m[],                 // coords   <-
  std::set<const Node *, nodeCompare> & s,   // (Node map, compare) ->
  std::vector<NodeExt> & extNode_vect               // extNode vector   -> 
){  // ==================================================
  
  // From med-mesh coordinates
//   int ndim = d->getNumberOfTuples();
//   int dim = d->getNumberOfComponents();
  // setting.clear() node ext
  extNode_vect.clear(); if (n_element_med == 0)  return; 
  extNode_vect.resize(n_element_med);
//   s.clear();
//   unsigned int i, idim;
  for (int i=0; i<n_element_med; i++) {
    // set  coordinates to Node ext
    for (int idim=0; idim<dim_med; idim++)         extNode_vect[i](idim) = xyz_m[i+idim*n_element_med];// d->getIJ(i, idim);
    for (int idim=dim_med; idim<DIMENSION; idim++) extNode_vect[i](idim) = 0.0;
    extNode_vect[i]._id = i;
    extNode_vect[i]._id2 = i;
//     extNode_vect[i].processor_id(345);
    // insert ext Node to s structure
    s.insert(&(extNode_vect[i]));
    
  }
   return;
}

// =======================================================
// This function set the xyz vector into ExtNode victor n
void centerToSetOfNodes2(
  const MGMesh &mesh_femus,
  const int dim,
  const int n_element,
  const int Level,
  const double  xyz_m[],                 // coords   <-
  std::set<const Node *, nodeCompare> & s,   // (Node map, compare) ->
  std::vector<NodeExt> & extNode_vect               // extNode vector   -> 
){  // ==================================================
  
  // From med-mesh coordinates
//   int ndim = d->getNumberOfTuples();
//   int dim = d->getNumberOfComponents();
  // setting.clear() node ext
  extNode_vect.clear(); if (n_element == 0)  return; 
  extNode_vect.resize(n_element);
//   s.clear();
//   unsigned int i, idim;
//   int Level=nlevels;
 
  //  2bB element interpolation over the fine mesh -----------------------

  int ielemn=0;
  for (int iproc=0; iproc<mesh_femus._n_subdom; iproc++) {
    int nel_b =mesh_femus._off_el[0][mesh_femus._NoLevels*iproc+Level];
    int nel_e= mesh_femus._off_el[0][Level+mesh_femus._NoLevels*iproc+1];
    for (int iel=nel_b; iel <nel_e; iel++) {
     
     // set  coordinates to +init_lev_elNode ext
    for (int idim=0; idim<dim; idim++)         extNode_vect[iel](idim) = xyz_m[iel+idim*n_element];// d->getIJ(i, idim);
    for (int idim=dim; idim<DIMENSION; idim++) extNode_vect[iel](idim) = 0.0;
    extNode_vect[iel]._id = iel;
    extNode_vect[iel]._id2 = iel;
//     extNode_vect[i].processor_id(345);
    // insert ext Node to s structure
    s.insert(&(extNode_vect[iel]));
    
  } // ---- end iel -------------------------
  } // end iproc
   return;
}




// // ===================================================================
void CleanSetOfNodes(std::set<const Node *, nodeCompare> & s){
  typename std::set<const Node *, nodeCompare>::iterator itS;
  for (itS = s.begin(); itS != s.end(); itS++) {
    delete *itS;
  }
   s.clear();
}