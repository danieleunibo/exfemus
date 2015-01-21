#ifndef __ParallelMeshExtended__
#define __ParallelMeshExtended__

#include <mpi.h>
#include "MGMesh.h"
#include "MGFE_conf.h"        // FEM conf
// #include "libmesh/parallel_mesh.h"
#include <vector>
namespace ParaMEDMEM {
  class MEDCouplingUMesh;
  class MEDCouplingFieldDouble;
  class MEDLoader;
}


class ParallelMeshExtended : public MGMesh {
  protected:
  int _num_meshes;
  std::map <int, int> _bd_group;  // group number <-> number of pts
//   std::vector< const ParaMEDMEM::MEDCouplingUMesh * >  _meshes;
  const ParaMEDMEM::MEDCouplingUMesh & _mesh;
  MPI_Comm                             _comm;
 
  std::vector<int>                     _node_id;
  std::vector<int>                     _proc_id;
  std::vector<int>                     _elem_id;
//   std::vector<int>                     _bd_elem_id;
  int                                  _meshDim;

  
public:
  
   std::vector<int>                     _bc_id;
  std::vector<int>                     _mat_id;
  
  // Constructor-Destructor
  ParallelMeshExtended(const ParaMEDMEM::MEDCouplingUMesh & mesh,
                       const ParallelObjectM &comm,
                       MGUtils& mgutils,
                       MGGeomEl& geomel,
                       MPI_Comm comm1 = MPI_COMM_NULL
                      );
//   ParallelMeshExtended(const ParaMEDMEM::MEDCouplingUMesh & mesh,
// 	               MPI_Comm comm = MPI_COMM_NULL);
  virtual ~ParallelMeshExtended();

  const std::vector<int> & getNodeID() const {return _node_id;}
  const std::vector<int> & getProcID() const {return _proc_id;}
  const std::vector<int> & getElemID() const {return _elem_id;}
  void read_bc_id(int Level);
  void read_mat_id(int Level);
  void print_med(int Level,std::string filename);
  int meshDimension() {return _meshDim; }
  

};

#endif
