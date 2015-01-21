#ifndef __FEMUS__
#define __FEMUS__

#include <mpi.h>
#include <vector>
#include "MEDCouplingRefCountObject.hxx"
#include "Femus_BCTypes.h"

namespace ParaMEDMEM {
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
}

// namespace libMesh {
//   class EquationSystems;
//   class LinearImplicitSystem;
// }

class ParallelMeshExtended;
// class _LibMeshProblem;
class LibMeshFunction;

class MGFemusInit;
class EquationSystemsExtended;
class MGSystem;
class MGUtils;
class MGMesh;
class MGGeomEl;
class MGFEMap;
class MGTimeLoop;



class FEMUS {
  // ==========================================================================
  //=========   data  =========================================================
protected:

  //data communication (defined in Constructor) -------------------------------
  MGFemusInit *    _start;          // start function
  MPI_Comm         _comm;           // communicator
  bool             _local_MPI_Init; // initial mpi flag

  // param and file data (defined in init(.,.)) ------------------------------------------
  MGUtils *     _mg_utils;     // param and files

  // fem ------- (defined in init(.,.))  --------------------------------------------------
  MGGeomEl *                 _mg_geomel;    // set element type
  MGFEMap  *                 _mg_femap;     // fem map

  // data meshes --------------------------------------------------------------------------
  int                                    _num_mgmesh;   // num of MGmeshes
  std::vector<ParallelMeshExtended*>     _mg_mesh;      // FEMus-mesh
  ParaMEDMEM::MEDCouplingUMesh *         _med_mesh;     // Med-mesh
  // system data ---------------------------------------------------------------------------
  //   _LibMeshProblem * _problem;
  std::vector<EquationSystemsExtended*>  _mg_equations_map; // system
  MGTimeLoop *                           _mg_time_loop;     // transient

  // Interface mesh  ----------------------------------------------------------
  struct interface_mesh {                        ///< interface mesh structure
    int id;                                       ///< name identity (int name)
    std::string name;                             ///< name
    ParaMEDMEM::MEDCouplingUMesh * support;       ///< Med-submesh
  };
  std::vector<interface_mesh> _interface_mesh_vect;///< interface-mesh vector

 // User interface bc --------------------------------------------------------- 
  struct UsersBC {                      ///< interface user bc structure
    std::string name;                            ///< name
    std::string type;                            ///< type: Marked (UnMarked)
    int from_cmp;                                ///< UsersBC initial cmp
    int n_cmp;                                   ///< UsersBC n of cmp
    int order_cmp;                               ///< UsersBC order cmp
    bool isAnalytic;                             ///< true (false)
    std::string equation;                        ///< equation string
    ParaMEDMEM::MEDCouplingFieldDouble * field;  ///< numerical field
    bool on_nodes;                             ///< node field (true)
  };
  std::vector<UsersBC>  _boundary;      ///< interface user bc vector
 
  // ==========================================================================
  //=========  pubblic functions  ============================================
public:

  FEMUS();
  FEMUS(MPI_Comm comm);
  void init_param(MGUtils &mgutils );
  void init_fem(MGGeomEl & mggeomel,MGFEMap & mgfemap);
  void clean();
  ~FEMUS();
  void terminate();

  void solve();

//   void setType(const std::string & pbName);
//   void setMesh(const std::string & meshFileName);
  void setMesh();
  void setMedMesh(const std::string & dataFile);
  void setSystem(const std::vector<NS_FIELDS> & pbName);
  
  void init_interfaces(
    const std::string & medfilename,
    const int index_mgmesh=0,
    const int index_medmesh=0);

  void setSource(const ParaMEDMEM::MEDCouplingFieldDouble * f);
  void setAnalyticSource(const std::string & f);
  void setAverageSources(int n,double val[]);

  void setAnalyticBoundaryValues(
    const std::string & name,
      int n_cmp,
    const std::string & 
    typeBC,const std::string & f);
  
  void setFieldBoundaryValues(
    const std::string & name,
    int n_cmp,
    const std::string & bcType,
    const ParaMEDMEM::MEDCouplingFieldDouble * bcField);
  

  ParaMEDMEM::MEDCouplingFieldDouble * getOutputField(const std::string & vName) const;
  void saveOutputField(const std::string & vName, const std::string & prefix) const;
  ParaMEDMEM::MEDCouplingFieldDouble * getInputFieldTemplate(const std::string &name);

  ParaMEDMEM::MEDCouplingUMesh * get_mesh_med() { return _med_mesh;}

  // interface vector function (_interface_mesh_vect)
  ParaMEDMEM::MEDCouplingUMesh * get_interface_mesh(int i) { 
    return _interface_mesh_vect[i].support;
  }
  int get_interface_byname(const std::string &name) { return  search_idxBC_byid(search_idBC_byname(name));}
  int get_interface_byname(int name) { return  search_idxBC_byid(name);}

  std::vector<std::string> get_interface_Names();

  // ============================================================================
/// This function gets the values of the variable with names in
/// std::vector<char *> vName on boundary with identity id
ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary(// field (out) 
  const std::string & nameBC,                    ///< boundary name (char*) (in)
  const std::string & systemName                 ///< system name      (in)
);
 // ============================================================================
/// This function gets all the values on boundary with identity id
ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary(  // field (out) 
  const std::string & nameBC,                     ///< boundary name (char*) (in)
 const std::string & systemName,                  ///< system name           (in)
  int iCmp                                        ///< component             (in)
);
    
//   ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary(const std::string & name,  std::vector<char *> vName) ;
//   ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary_nodes(const std::string & name,  std::vector< char *>) const;

//   void setDebugFile(const std::string & s);
//   void endDebug();
  
  // UsersBC structure functions
  void set_param_BCzone(); 
  void  set_param_BCzone(FEMUS & P1);
  void init_BCzone(std::vector<std::string> &bcNames);
  void set_interface_userbc();

protected:
//      void setMesh(const ParaMEDMEM::MEDCouplingUMesh * m);
  // name -> id (location in the interface vector) (_interface_mesh_vect)
  int search_idxBC_byid( // vector index (return)
    const int  inameBC          //  id boundary name (int)   <-
  )const;
  int search_idBC_byname(        // id name (return)
    const std::string & nameBC          //  boundary name  (in)
  ) const;
//    int searchBoundaryCondition(const std::string &) const;
  std::string get_interface_Name(int i);

  // function from FEMUS class to EquationSystemsExtended
  void setBCType(
    const std::string & bcName, 
    const std::string & typeBC);
  
  void setAnalyticBCValues(
    const std::string & bcName,
      int n_cmp, const std::string & f);
  
  void setFieldBCValues(
    const std::string & bcName,
    int n_cmp,  const ParaMEDMEM::MEDCouplingFieldDouble * f);
};

#endif
