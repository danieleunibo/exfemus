#ifndef _EQNEXTENDED__
#define _EQNEXTENDED__

#include <map>
#include "libmesh/id_types.h"
#include "libmesh/equation_systems.h"
#include "Debug.h"

namespace ParaMEDMEM {
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
}

class LibMeshFunction;

#include "BCTypes.h"

class EquationSystemsExtended : public libMesh::EquationSystems
{
  protected:
  // system 
  int _nComp;
  
  // bc to use in assembling
public:
  std::map<libMesh::boundary_id_type, LibMeshFunction *> _bc;
  std::map<libMesh::boundary_id_type, BCType> _bc_type;
protected:
  // sources to use in assembling
  LibMeshFunction *_source;  // field source
  int _nAverageSources;     // constant sources
  double *_AverageSources;  // constant source vector
  
public:
  // Constructor-Desctructor
  EquationSystemsExtended(libMesh::MeshBase & mesh, int nComp);
  virtual ~EquationSystemsExtended();

  int getNumberOfComponents() const { return _nComp; }
  
  // source  to use in assembly --------------------------------------------
  void setSource(const ParaMEDMEM::MEDCouplingUMesh * mesh,const char *s);
  void eraseSource();
  LibMeshFunction * getSource();
  
  // double _AverageSources[_nAverageSources]
  // setAverageSources ->   _AverageSources   <- getAverageSources
  void setAverageSources(int n_AverageSources,double Average[]){
    _nAverageSources=n_AverageSources;
    delete []_AverageSources; _AverageSources=new double[_nAverageSources];
    for(int i=0;i<_nAverageSources;i++) _AverageSources[i]=Average[i];
    return;}
  double getAverageSources(int i){return _AverageSources[i];}
   
  // boundary conditions to use in assembly -------------------------------
  libMesh::boundary_id_type addBC(const ParaMEDMEM::MEDCouplingUMesh * b);
  void setBCType (libMesh::boundary_id_type id, BCType type);
  void setBC     (libMesh::boundary_id_type id, 
		  const char *s);
  void setBC     (libMesh::boundary_id_type id, 
		  const ParaMEDMEM::MEDCouplingFieldDouble *f);
  void setNodeBC(libMesh::boundary_id_type id, 
		 const ParaMEDMEM::MEDCouplingFieldDouble *f);
  LibMeshFunction * getBC(libMesh::boundary_id_type id);
  BCType getBCType(libMesh::boundary_id_type id);
  void eraseBC(libMesh::boundary_id_type id);

};

#endif
