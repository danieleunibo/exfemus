#ifndef __MeshExtendedM__
#define __MeshExtendedM__

#include "MGMesh.h"
#include <vector>
#include "Solverlib_conf.h" 

class MeshExtended : public MGMesh {
public:
   std::vector<int>                      _bc_id;
   std::vector<int>                     _mat_id;
  
  // Constructor-Destructor
  MeshExtended(        const ParallelObjectM &comm,
                       MGUtils& mgutils,
                       MGGeomEl& geomel
                      );
  
  virtual ~MeshExtended();

  void read_bc_id(int Level);
  void read_mat_id(int Level);
  
#ifdef HAVE_MED
  void print_med(int Level,std::string filename);
#endif
  
};

#endif
