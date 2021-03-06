#define FDEBUG 0
#include <set>
#include "libmesh/equation_systems.h"
#include "libmesh/node.h"
#include "libmesh/boundary_info.h"
#include "EquationSystemsExtended.h"
#include "Debug.h"
#include "LibMeshUtil.h"
#include "LibMeshFunction.h"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "BCTypes.h"
#include "MEDCouplingRemapper.hxx"

using namespace libMesh;
using namespace ParaMEDMEM;

// ========================================================================
/// Constructor 
EquationSystemsExtended::EquationSystemsExtended(MeshBase & mesh, int nComp) : 
  libMesh::EquationSystems(mesh), _nComp(nComp), _source(new LibMeshFunction) 
{
  int _nAverageSources=1;
 _AverageSources=new double[_nAverageSources];
 for(int i=0;i<_nAverageSources;i++) _AverageSources[i]=0.;
}

// ========================================================================
EquationSystemsExtended::~EquationSystemsExtended() {
  if (_source) delete _source;

  std::map<boundary_id_type, LibMeshFunction *>::iterator it = _bc.begin();
  for (; it != _bc.end(); it++) {
    if (it->second) delete it->second;
    it->second = NULL;
  }
}
// ========================================================================
LibMeshFunction * EquationSystemsExtended::getSource(){
  return _source;
}
// ========================================================================
void EquationSystemsExtended::setSource
(const MEDCouplingUMesh * medmesh, const char *s) {
  MeshBase & mesh = get_mesh();

  _source->setSupport(&mesh, medmesh);
  _source->setCode(s, _nComp);
  if (FDEBUG) fDebugPos << "setSource (string " << s << ")" << std::endl;
  if (FDEBUG) _source->printOn(fDebug);
}
// ========================================================================
void EquationSystemsExtended::eraseSource() {
  if (_source) {
    delete _source;
    _source = NULL;
  }
}
// ========================================================================
void EquationSystemsExtended::setBCType(boundary_id_type id, BCType type)
{
  if (FDEBUG) fDebugPos << "id " << id << " BCType : " << type << std::endl;
  _bc_type[id] = type;
}
// =======================================================================
BCType EquationSystemsExtended::getBCType(boundary_id_type id){
  return (_bc_type.find(id) == _bc_type.end()) ? Undefined : _bc_type[id];
}
// =======================================================================
void EquationSystemsExtended::setBC(
  boundary_id_type id, 
  const MEDCouplingFieldDouble *f
) {// =======================================================================
  
  LibMeshFunction * fct = getBC(id);
  if (fct == NULL) return;

  MeshBase & mesh = get_mesh();
  std::map<std::pair<int, int>, int>::iterator itF = fct->FaceID().begin();
  for (; itF != fct->FaceID().end(); itF++) {
    int elem_id = itF->first.first;
    int side = itF->first.second;
    mesh.boundary_info->add_side(elem_id, side, id);
  }
 
  fct->setField(f);
  if (FDEBUG) {
     fDebugPos << "setBC (type " << _bc_type[id] 
	           << " field " << f->getName() << ")" << std::endl;
     fct->printOn(fDebug);
  }
}

// ===========================================================
void EquationSystemsExtended::setNodeBC(
  boundary_id_type id, 
  const MEDCouplingFieldDouble *f
){//==========================================================
  
  LibMeshFunction * fct = getBC(id);
  if (fct == NULL) return;
   
  if (FDEBUG) fDebug.array(f->getArray(), "ici f 3");
  
  MEDCouplingFieldDouble * ff=MEDCouplingFieldDouble::New(ON_CELLS);
  const MEDCouplingUMesh * mesh = fct->getSupport();
  ff->setMesh(mesh);
  int nCells = mesh->getNumberOfCells();
  int nComp = f->getNumberOfComponents();
  DataArrayDouble * d = DataArrayDouble::New();
  const DataArrayDouble * df = f->getArray();
  d->alloc(nCells, nComp);

  const DataArrayInt * c = mesh->getNodalConnectivity();
  const DataArrayInt * cI = mesh->getNodalConnectivityIndex();

  for (int iCell=0; iCell < nCells; iCell++) {

    INTERP_KERNEL::NormalizedCellType t
      = mesh->getTypeOfCell (iCell);
    int n_nodes_elem;

    switch (t) {
    case INTERP_KERNEL::NORM_SEG2:
      n_nodes_elem = 2;
      break;
          case INTERP_KERNEL::NORM_SEG3:
      n_nodes_elem = 3;
      break;
    case INTERP_KERNEL::NORM_TRI3:
      n_nodes_elem = 3;
      break;
    case INTERP_KERNEL::NORM_TETRA4:
      n_nodes_elem = 4;
      break;
    case INTERP_KERNEL::NORM_QUAD4:
      n_nodes_elem = 4;
      break;
        case INTERP_KERNEL::NORM_QUAD8:
      n_nodes_elem = 8;
      std::cout << "\n I am quad8  \n";
      break;  
      case INTERP_KERNEL::NORM_QUAD9:
      n_nodes_elem = 9;
      std::cout << "\n I am quad9  \n";
      break;  
    case INTERP_KERNEL::NORM_TETRA10:
      n_nodes_elem = 10;
      std::cout << "\n I am tetra10  \n";
      break;  
      case INTERP_KERNEL::NORM_HEXA8:
      n_nodes_elem = 8;
      std::cout << "\n I am hexa8  \n";
      break;
     case INTERP_KERNEL::NORM_HEXA20:
      n_nodes_elem = 20;
      std::cout << "\n I am hexa20  \n";
      break;
      case INTERP_KERNEL::NORM_HEXA27:
      n_nodes_elem = 27;
      std::cout << "\n I am hexa27  \n";
      break;
    default:
      n_nodes_elem = 0;
    }

    int cstart = cI->getIJ(iCell, 0);
    //   int cend = cI->getIJ(iCell+1, 0);

    double v = 0.0;
    for (int i=0; i<n_nodes_elem; i++) {
      double x = df->getIJ(c->getIJ(cstart+1+i, 0), 0);
      v += x;
    }
    v /= n_nodes_elem;
    d->setIJ(iCell, 0, v);
  }

  ff->setArray(d);
  ff->checkCoherency();

  if (FDEBUG) fDebugPos << "setNodeBC (type " << _bc_type[id] 
	    << " field " << ff->getName() << ")" << std::endl;
  setBC(id, ff); 
}

// =========================================================================
void EquationSystemsExtended::setBC(
  boundary_id_type id, 
  const char *s
) { // ======================================================================
  
  LibMeshFunction * fct = getBC(id);
  if (fct == NULL) return;
  fct->setCode(s, _nComp);

  MeshBase & mesh = get_mesh();
  std::map<std::pair<int, int>, int>::iterator itF = fct->FaceID().begin();
  for (; itF != fct->FaceID().end(); itF++) {
    int elem_id = itF->first.first;
    int side = itF->first.second;
    mesh.boundary_info->add_side(elem_id, side, id);
  }
  if (FDEBUG) fDebugPos << "setBC (type " << _bc_type[id] 
	    << " string " << s <<")" << std::endl;
  if (FDEBUG) fct->printOn(fDebug);
}


// =======================================================================
LibMeshFunction * EquationSystemsExtended::getBC(
  boundary_id_type id
) { // ===================================================================
  return _bc.find(id) == _bc.end() ? NULL : _bc[id];
}

// =======================================================================
void EquationSystemsExtended::eraseBC(
  boundary_id_type id
) {// =======================================================================
  _bc.erase(id);
}
// =======================================================================
boundary_id_type EquationSystemsExtended::addBC(
  const MEDCouplingUMesh * b
) {// =======================================================================
 
  unsigned int i = 0;
  for (i=0; i<_bc.size() + 1; i++) if(_bc.find(i) == _bc.end()) break;
  
  MeshBase & mesh = get_mesh();
  LibMeshFunction *f = new LibMeshFunction;
  f->setSupport(&mesh, b);
  _bc[i] = f;

  return boundary_id_type(i);
}
