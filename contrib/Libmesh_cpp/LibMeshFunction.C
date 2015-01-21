#include "LibMeshUtil.h"
#include "LibMeshFunction.h"
#include "libmesh/libmesh.h"
#include "libmesh/point.h"
#include "libmesh/elem.h"
#define FDEBUG 0
#include "Debug.h"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingUMesh.hxx"
#include <iomanip>
#include <limits>

using namespace ParaMEDMEM;
using namespace libMesh;

LibMeshFunction::~LibMeshFunction()
{
  if (_f) _f->decrRef();
}

void LibMeshFunction::setSupport(const MeshBase * mesh, 
				 const MEDCouplingUMesh * support ) 
{
  _mesh = mesh;
  _support = support;

  std::set<const Node *, nodeCompare> s1, s2;
  std::vector<NodeExt> n;

  // Corresponding nodes
  
  const DataArrayDouble * d = _support->getCoords();
  DataArrayToSetOfNodes(d, s1, n);
  const MeshBase::const_node_iterator end_nl 
    = _mesh->active_nodes_end();
  for (MeshBase::const_node_iterator nl = _mesh->active_nodes_begin(); 
       nl != end_nl; nl++) s2.insert(*nl);
  _nodeID.clear();
  commonPoints(_nodeID, s1, s2);
  s1.clear();

  d = _support->getBarycenterAndOwner();
  DataArrayToSetOfNodes(d, s1, n);

  if (_mesh->mesh_dimension() == _support->getMeshDimension()) {
    // Corresponding elements
    std::set<const Node *, nodeCompare> s2;
    const MeshBase::const_element_iterator end_el 
      = _mesh->active_local_elements_end();
    for (MeshBase::const_element_iterator el 
	   = _mesh->active_local_elements_begin() ; 
	 el != end_el; ++el)
      {
	const Elem* elem = *el;
	s2.insert(new NodeExt(elem->centroid(), elem->id(), 0));
      }
    commonPoints(_elemID, s1, s2);
    CleanSetOfNodes(s2);
  }
  else {
    // Corresponding faces
    std::set<const NodeExt *, nodeCompare> s2;
    const MeshBase::const_element_iterator end_el 
      = _mesh->active_local_elements_end();
    for (MeshBase::const_element_iterator el 
	   = _mesh->active_local_elements_begin() ; 
	 el != end_el; ++el)
      {
	const Elem* elem = *el;
	for (unsigned int side =0 ; side < elem->n_sides(); side++)
	  if (elem->neighbor(side) == NULL) {
	    NodeExt * q = new NodeExt(elem->build_side(side)->centroid(), 
				      elem->id(), side);
	    s2.insert(q);
	  };
      }

    commonPoints(_faceID, s1, s2);
    CleanSetOfNodes(s2);  
  }
  //d->decrRef();
}

void LibMeshFunction::setCode(const char *code, int nComp)
{
  if (_f) _f->decrRef();

  TypeOfField type = ON_CELLS;
  std::vector<std::string> vars(_support->getSpaceDimension());
  if (_support->getSpaceDimension() > 0)
    vars[0] = "x";
  if (_support->getSpaceDimension() > 1)
    vars[1] = "y";
  if (_support->getSpaceDimension() > 2)
    vars[2] = "z";

  _f = _support->fillFromAnalytic3(type, nComp, vars, code);
  _f->setName(code);

  if (FDEBUG) {
    int i,j;
    for (i=0; i<_f->getNumberOfTuples(); i++) {
      fDebugPos << " " << i << ": ";
      for (j=0; j<nComp; j++)
	     fDebug << " " << _f->getIJ(i,j);
      fDebug << std::endl;
    }
  }
}

void LibMeshFunction::setField(const MEDCouplingFieldDouble *f)
{
  if (_f) _f->decrRef();

  TypeOfField type = f->getTypeOfField();
  _f = MEDCouplingFieldDouble::New(type);
  _f->setMesh(_support);
  DataArrayDouble *array = DataArrayDouble::New();
  int nTuples;
  int nComp;

  const DataArrayDouble * arrayf = f->getArray(), *coord1, *coord2;
  if (FDEBUG) fDebug.array(arrayf, "arrayf");

  nTuples = _support->getNumberOfCells();
  nComp = f->getNumberOfComponents();
  array->alloc(nTuples, nComp);

  coord1 = _support->getBarycenterAndOwner();
  coord2 = f->getMesh()->getBarycenterAndOwner();

  array->fillWithValue(std::numeric_limits<double>::max());
  std::map<int, int> corresp;
  commonPoints(corresp, coord1, coord2);

  for (std::map<int, int>::iterator it = corresp.begin();
       it != corresp.end(); it++) {
    for (int j=0; j<nComp; j++) {
      if (FDEBUG) fDebugPos << " " << it->first << " -> " << it->second 
		<< " val = " << arrayf->getIJ(it->first, j) << std::endl;
      array->setIJ(it->second, j, arrayf->getIJ(it->first, j));
    }
  }
  _f->setArray(array);
  _f->setName(f->getName());
  array->decrRef();
  _f->checkCoherency();
}

MEDCouplingFieldDouble * LibMeshFunction::getField(char const* name)
{
  if (FDEBUG) fDebugPos << "LibMeshFunction::getField _f = " << _f << std::endl;
  if (!_f) return NULL;
  if (FDEBUG) fDebugPos << "LibMeshFunction::getField _f = " << _f << std::endl;
  _f->incrRef();
  if (FDEBUG) fDebugPos << "LibMeshFunction::getField " << name << std::endl;
  _f->setName(name);
  return _f;
}

void LibMeshFunction::eval (int id, std::vector<double> & val)
{
  if (_f == NULL) {
    for (int i = 0; i<val.size(); i++)
      val[i] = 0.0;
  }
  else {
    for (int i = 0; i<val.size(); i++)
      val[i] = _f->getIJ(_elemID[id], i);
  }
}

void LibMeshFunction::eval(int id, int side, std::vector<double> & val)
{
  if (_f == NULL) {
    for (int i = 0; i<val.size(); i++)
      val[i] = 0.0;
  }
  else {
    int mid = _faceID[std::pair<int, int>(id, side)];
    for (int i = 0; i<val.size(); i++)
      val[i] = _f->getIJ(mid, i);
  }
}

void LibMeshFunction::printOn(std::ostream & out) const
{
  out << std::endl 
      << "==================================================" << std::endl
      << "LibMeshFunction printOn" << std::endl << std::endl;

  out << "node correspondance (libMesh <-> MEDCoupling)" << std::endl;
  std::map<int, int>::const_iterator itN = _nodeID.begin();
  for (; itN != _nodeID.end(); itN++)
    out << "  " << std::setw(4) << itN->first 
	<< " <-> " << std::setw(4) << itN->second << std::endl;
   out << std::endl;

  out << "element correspondance (libMesh <-> MEDCoupling)" << std::endl;
  std::map<int, int>::const_iterator itE = _elemID.begin();
  for (; itE != _elemID.end(); itE++)
    out << "  " << std::setw(4) << itE->first 
	<< " <-> " << std::setw(4) << itE->second << std::endl;
 
  out << std::endl;

  out << "face correspondance (libMesh <-> MEDCoupling)" << std::endl;
  std::map<std::pair<int,int>, int>::const_iterator itF = _faceID.begin();
  for (; itF != _faceID.end(); itF++)
    out << " (" << std::setw(4) << itF->first.first 
	<< ", " << std::setw(4) << itF->first.second << ")"
	<< " <-> " << std::setw(4) << itF->second << std::endl;
 
  out << std::endl;

  if (!_f) {
    out << "NULL" << std::endl;
    out << "==================================================" << std::endl;
    return;
  }
  TypeOfField type = _f->getTypeOfField();
  const DataArrayDouble * d;
  int n;
  const MEDCouplingUMesh * m 
    = dynamic_cast<const MEDCouplingUMesh *>( _f->getMesh());

  out << "name : " << _f->getName() << std::endl;

  switch (type) {
  case ON_NODES:
    out << "type : ON_NODES" << std::endl;
    d = m->getCoords();
    d->incrRef();
    n = d->getNumberOfTuples();
    out << " " <<n << " nodes";
    break;
  case ON_CELLS:
    out << "type : ON_CELLS" << std::endl;
    d = m->getBarycenterAndOwner();
    n = d->getNumberOfTuples();
    out << " " << n << " cells";
    break;
  }

  DataArrayDouble * v = _f->getArray();
  out << " (" << v->getNumberOfTuples() << ")" << std::endl;

  int nc = d->getNumberOfComponents();
  for (int i=0; i<n; i++) {
    for (int j=0; j<nc; j++)
      out << (j == 0 ? "f( " : ", ") 
	  << std::setw(10) << std::fixed << d->getIJ(i,j);
    out << ") = ";
    for (int j=0; j<nc; j++)
      out << (j == 0 ? "( " : ", ") 
	  << std::setw(10) << std::fixed << v->getIJ(i,j);
    out << std::endl;
  }
  out << std::endl << std::endl;
  out << "==================================================" << std::endl;

  //d->decrRef();
}
