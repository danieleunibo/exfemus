#ifndef __LIBMESH_ANALYTICFUNCTION__
#define __LIBMESH_ANALYTICFUNCTION__

#include <iostream>
#include <map>
#include <vector>

// #include "libmesh/mesh_base.h"
// #include "libmesh/mesh.h"

#include "InterpKernelExprParser.hxx"
#include "MEDCouplingRefCountObject.hxx"


namespace ParaMEDMEM {
  class MEDCouplingFieldDouble;
  class MEDCouplingUMesh;
}


class LibMeshFunction {

 public:
  LibMeshFunction() : 
    _support(NULL),  _f(NULL), _mesh(NULL)
  {
  }
  ~LibMeshFunction();
  void eval(int, std::vector<double> & v);
  void eval(int, int, std::vector<double> & v);

  void setSupport(const libMesh::MeshBase * mesh,
		  const ParaMEDMEM::MEDCouplingUMesh * support);
  void setCode(const char *code, int nComp);
  void setField(const ParaMEDMEM::MEDCouplingFieldDouble *f);
  std::map<std::pair<int, int>, int> & FaceID() { return _faceID; }
  std::map<int, int> & NodeID() { return _nodeID; }
  std::map<int, int> & ElemID() { return _elemID; }

  ParaMEDMEM::MEDCouplingFieldDouble * getField(const char * vName);
  const ParaMEDMEM::MEDCouplingUMesh * getSupport() 
  {
    return _support;
  }

  void printOn(std::ostream & out) const;

  
 protected:

  const libMesh::MeshBase * _mesh;
  const ParaMEDMEM::MEDCouplingUMesh * _support;
  ParaMEDMEM::MEDCouplingFieldDouble *_f;

  std::map<int, int> _nodeID;
  std::map<std::pair<int, int>, int> _faceID;
  std::map<int, int> _elemID;
};

#endif

