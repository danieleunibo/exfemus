#ifndef ___LIBMESHPROBLEM__
#define ___LIBMESHPROBLEM__

#include <mpi.h>
#include <map>
#include "BCTypes.hxx"

namespace ParaMEDMEM {
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
}
// class LibMeshFunction;

// namespace libMesh {
//   class System;
// }

// class EquationSystemsExtended;
// class ParallelMeshExtended;


class _LibMeshProblem {

public:
  // Constructor-destructor
  _LibMeshProblem(const char *pbName, int nComp, MPI_Comm comm);
  virtual ~_LibMeshProblem();
//   virtual void terminate();
//  
//   virtual void setMesh(const ParaMEDMEM::MEDCouplingUMesh * mesh);
//   const ParaMEDMEM::MEDCouplingUMesh * getMesh() const { return med_mesh; }
//   virtual libMesh::boundary_id_type  defineBoundary(const ParaMEDMEM::MEDCouplingUMesh * bc);
//   
//   virtual void solve ();
//  
//   virtual void setAnalyticSource(const char * f);
//   virtual void setSource(const ParaMEDMEM::MEDCouplingFieldDouble * f);
//   virtual void setAverageSources(int n, double val[]);
//    
//   virtual void setBoundaryConditionType(int, const char * typeBC);
//   virtual void setAnalyticBoundaryValues(int, const char * f);
//   virtual void setFieldBoundaryValues(int, const ParaMEDMEM::MEDCouplingFieldDouble * f);
//   
//   virtual ParaMEDMEM::MEDCouplingFieldDouble * getOutputField(const char * vName) const = 0;
//   ParaMEDMEM::MEDCouplingFieldDouble * getBoundaryValues(const ParaMEDMEM::MEDCouplingUMesh * bc) const;
//    ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary(int, const char *) const;
//   ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary(int,  std::vector< char *>) const;
//   ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary_nodes(int,  std::vector< char *>) const;
//  
 
protected:

  MPI_Comm _comm;                              // communicator
//   libMesh::LibMeshInit *_initLibMesh;
  std::string problemName;
   
  const ParaMEDMEM::MEDCouplingUMesh * med_mesh;
  ParallelMeshExtended * mesh;
 
  EquationSystemsExtended * equation_systems;
  virtual void _initSystem() = 0;
 
  int _nComp;
//   libMesh::System * _getVarIds( std::vector<char *> name, std::vector<int> &) const;
//    libMesh::System * _getVarIds( const char * name, std::vector<int> &) const;
//   std::map<std::string, BCType> _BCTypes;
//   
//   std::vector <void (*)(libMesh::EquationSystems&)> _post_processing_functions;
//   ParaMEDMEM::MEDCouplingFieldDouble *  _getOutputField(const char * vName,ParaMEDMEM::TypeOfField) const;
};


#endif
