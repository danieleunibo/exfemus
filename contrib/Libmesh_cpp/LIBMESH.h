#ifndef __LIBMESH__
#define __LIBMESH__

#include <mpi.h>
#include <vector>
#include "MEDCouplingRefCountObject.hxx"
#include "BCTypes.h"

namespace ParaMEDMEM {
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
}

namespace libMesh {
  class EquationSystems;
  class LinearImplicitSystem;
}

class ParallelMeshExtended;
class _LibMeshProblem;
class LibMeshFunction;

class LIBMESH {

public:
  LIBMESH();
  LIBMESH(MPI_Comm comm);
  ~LIBMESH();
  void terminate();
  void setType(const std::string & pbName);
  void setMesh(const std::string & meshFileName);
  void solve();
  ParaMEDMEM::MEDCouplingFieldDouble * getOutputField(const std::string & vName) const;
  void saveOutputField(const std::string & vName, const std::string & prefix) const;
  void setSource(const ParaMEDMEM::MEDCouplingFieldDouble * f);
  void setAnalyticSource(const std::string & f);
  void setAverageSources(int n,double val[]);
  std::vector<std::string> getBoundaryNames();
  std::string getBoundaryName(int i);
  void setAnalyticBoundaryValues(const std::string & name,const std::string & typeBC,const std::string & f);
  void setFieldBoundaryValues(const std::string & name,const std::string & bcType,const ParaMEDMEM::MEDCouplingFieldDouble * bcField);
  ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary(const std::string & name, const std::string & vName) const;
  ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary(const std::string & name,  std::vector<char *> vName) const;
  ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary_nodes(const std::string & name,  std::vector< char *>) const;
  ParaMEDMEM::MEDCouplingFieldDouble * getInputFieldTemplate(const std::string &name);
  void setDebugFile(const std::string & s);
  void endDebug();
  
protected:
  // Problem  -------------------------------------------------
  //data
  MPI_Comm _comm;
  _LibMeshProblem * _problem;                 
  ParaMEDMEM::MEDCouplingUMesh * _mesh;
  // functions ----------------------------------------------
  void setBCType(const std::string & bcName, const std::string & typeBC);
  void setMesh(const ParaMEDMEM::MEDCouplingUMesh * m);
  
  
  // Boundary -----------------------------------------------------
  // data 
   struct sBC{
    int id; 
    std::string name;
    ParaMEDMEM::MEDCouplingUMesh * support;
  };
  std::vector<sBC> _bc;
  bool _local_MPI_Init;
  // functions
  int defineBoundary(const ParaMEDMEM::MEDCouplingUMesh * support);
  int searchBoundaryCondition(const std::string &) const;
  void setAnalyticBCValues(const std::string & bcName, const std::string & f);
  void setFieldBCValues(const std::string & bcName,const ParaMEDMEM::MEDCouplingFieldDouble * f);
  
 }; 


#endif
