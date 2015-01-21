#ifndef __ParallelMeshExtended__
#define __ParallelMeshExtended__

#include "MEDCouplingUMesh.hxx"

#include <mpi.h>
#include "libmesh/mesh.h"
#include "libmesh/parallel_mesh.h"
#include <vector>

class ParallelMeshExtended : public libMesh::ParallelMesh
{
public:
  ParallelMeshExtended(const ParaMEDMEM::MEDCouplingUMesh & mesh, 
	               MPI_Comm comm = MPI_COMM_NULL);
  virtual ~ParallelMeshExtended();

  const std::vector<int> & getNodeID() const { return _node_id; }
  const std::vector<int> & getProcID() const { return _proc_id; }
  const std::vector<int> & getElemID() const { return _elem_id; }

  int meshDimension() { return _meshDim; }

protected:
  const ParaMEDMEM::MEDCouplingUMesh & _mesh;
  MPI_Comm _comm;
  std::vector<int> _node_id, _proc_id, _elem_id;
  int _meshDim;
};

#endif
