#include "LibMeshUtil.h"
#include "MEDCouplingMemArray.hxx"
#include "Debug.h"

using namespace ParaMEDMEM;

const double nodeCompare::eps = 1e-12;

void commonPoints(std::map<int, int> & corresp,
		  const std::set<const libMesh::Node *, nodeCompare> & p, 
		  const std::set<const libMesh::Node *, nodeCompare> & q)
{
  std::set<const libMesh::Node *>::const_iterator pIt, qIt;

  for (qIt = q.begin(); qIt != q.end(); ++qIt) {
    if ((pIt = p.find(*qIt)) != p.end()) {
      corresp[(*qIt)->id()] = (*pIt)->id();
    }
  } 
}

void commonPoints(std::map<std::pair<int, int>, int> & corresp,
		  const std::set<const libMesh::Node *, nodeCompare> & p, 
		  const std::set<const NodeExt *, nodeCompare> & q)
{
  std::set<const libMesh::Node *>::const_iterator pIt;
  std::set<const NodeExt *>::const_iterator qIt;

  for (qIt = q.begin(); qIt != q.end(); ++qIt) {
    if ((pIt = p.find(*qIt)) != p.end()) {
      const NodeExt * q = *qIt;
      corresp[std::pair<int, int>(q->id(), q->id2())] = (*pIt)->id();
    }
  }
} 


void commonPoints(std::map<int, int> & corresp,
		  const DataArrayDouble * d1, 
		  const DataArrayDouble * d2)
{
  std::vector<NodeExt> n1, n2;
  std::set<const libMesh::Node *, nodeCompare> s1, s2;
  
  fDebugPos << "d1 = " << d1 << " d2 = " << d2 << std::endl;
  DataArrayToSetOfNodes(d1, s1, n1);
  DataArrayToSetOfNodes(d2, s2, n2);
  commonPoints(corresp, s1, s2);
}

const libMesh::Node * findNode(const std::set<const libMesh::Node *, nodeCompare> & p, int id)
{
  for (std::set<const libMesh::Node *>::const_iterator pIt = p.begin(); 
       pIt != p.end(); ++pIt)
    if ((*pIt)->id() == id) return *pIt;
  return NULL;
}

void DataArrayToSetOfNodes(const DataArrayDouble * d, 
 			   std::set<const libMesh::Node *, nodeCompare> & s,
 			   std::vector<NodeExt> & n)
{
  int ndim = d->getNumberOfTuples();
  int dim = d->getNumberOfComponents();

  n.clear();
  if (ndim == 0)
    return;
  
  n.resize(ndim);  
  unsigned int i, idim;
  for (i=0; i<ndim; i++) {
    
    for (idim=0; idim<dim; idim++)
      n[i](idim) = d->getIJ(i, idim);
    for (; idim<LIBMESH_DIM; idim++)
      n[i](idim) = 0.0;
    n[i].set_id() = i;
    n[i].processor_id(345);
    
    s.insert(&(n[i]));
  }
  
}

