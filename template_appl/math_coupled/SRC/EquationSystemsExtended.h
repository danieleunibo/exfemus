#ifndef _EQNEXTENDED__
#define _EQNEXTENDED__

#include <map>
// #include "libmesh/id_types.h"
// #include "libmesh/equation_systems.h"
#include "Femus_BCTypes.h"
#include "MGEquationsSystem.h"
// #include "Debug.hxx"
#include "ParallelMeshExtended.h"

typedef int boundary_id_type;

namespace ParaMEDMEM {
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
}
class LibMeshFunction;


// ==========================================================
// ====================                  ====================

class EquationSystemsExtended : public MGEquationsSystem {

  // data ------------------------------------------------------
public:
  // function (med/femus mesh,field) in the boundary
  std::map<boundary_id_type, LibMeshFunction *> _interfaceFunMap; // function boundary map
//   std::map<boundary_id_type, FemusBCType> _bc_type;       // type boundary map
protected:

  // system
  int _nComp; //< number of components
  // sources to use in assembling
  LibMeshFunction *_source;  // function source (med/femus mesh,field)
  int _nAverageSources;     // constant sources
  double *_AverageSources;  // constant source vector

  MPI_Comm _comm;                              // communicator
//   libMesh::LibMeshInit *_initLibMesh;
  std::string _problemName;

  const ParaMEDMEM::MEDCouplingUMesh * _med_mesh;
  ParallelMeshExtended * _mesh;
//   std::map<std::string, FemusBCType> _BCTypes;


  // functions --------------------------------------


//   EquationSystemsExtended * equation_systems;
  void _initSystem();

//   int _nComp;
  void _getVarIds(std::vector<char *> name, std::vector<int> &) const;
  void _getVarIds(const char * name, std::vector<int> &) const;
//   MGEquationsSystem * _getVarIds( std::vector<char *> name, std::vector<int> &) const;
//   MGEquationsSystem * _getVarIds( const char * name, std::vector<int> &) const;


//   std::vector <void (*)(MGEquationsSystem&)> _post_processing_functions{
//   return *NULL;}
  ParaMEDMEM::MEDCouplingFieldDouble *  _getOutputField(const char */* vName*/,ParaMEDMEM::TypeOfField) const {
    return NULL;
  }


  // =============================================================


public:

  // Constructor-Desctructor (similar ot MGSystem)
  EquationSystemsExtended(
    MGUtils& mgutils_in,
//                                 MGSystem& mgphys_in,
    ParallelMeshExtended& mgmesh_in,
    MGFEMap& mgfemap_in,
    int nComp,
    int npoint_data,
    int ncell_data
  );
//   // Constructor-Desctructor
//   EquationSystemsExtended(libMesh::MeshBase & mesh, int nComp);
  virtual ~EquationSystemsExtended();

  int getNumberOfComponents() const { return _nComp; }

  // source  to use in assembly --------------------------------------------
  void setSource(const ParaMEDMEM::MEDCouplingUMesh * mesh,const char *s);
  void eraseSource();
  LibMeshFunction * getSource();

  // double _AverageSources[_nAverageSources]
  // setAverageSources ->   _AverageSources   <- getAverageSources
  void setAverageSources(int n_AverageSources,double Average[]) {
    _nAverageSources=n_AverageSources;
    delete []_AverageSources; _AverageSources=new double[_nAverageSources];
    for(int i=0; i<_nAverageSources; i++) _AverageSources[i]=Average[i];
    return;
  }
  double getAverageSources(int i) {return _AverageSources[i];}
  void set_meshext(ParallelMeshExtended & mgmesh_ext) {_mesh= &mgmesh_ext;   }

//=============================================================================
  ///  This functon add an interface-function (med-mesh,femus-mesh, field)
  ///  in the interface-function boundary map (_interfaceFunMap):
  ///  the field in the interface-function is not assigned here  
  int  addfunction_BC(                      // map position            (return)
    const boundary_id_type  bd_name_int,    ///< boundary id (int)     (in)
    const ParaMEDMEM::MEDCouplingUMesh * b, ///< med-submesh           (in)
    const int from_cmp,                     ///< initial id component  (in)
    const int n_cmp,                        ///< n components          (in)
    const int order_cmp,                    ///< order component (2=quad;lin=1)
    const bool on_nodes                     ///< values on nodes (true)
  ); // =======================================================================
//   void setBCType (
//     boundary_id_type id,
//     FemusBCType type);
// ============================================================================  
  void setBC(
    boundary_id_type id,
    int n_cmp,
    const char *s);
  void setBC(
    boundary_id_type id,
    int n_cmp,
    const ParaMEDMEM::MEDCouplingFieldDouble *f);

  void setNodeBC(
    boundary_id_type id,
    int n_cmp,
    const ParaMEDMEM::MEDCouplingFieldDouble *f);

  void set_mesh_femus(ParallelMeshExtended & mg_mesh_femus_in) {
    _mesh=&mg_mesh_femus_in;
    return;
  }
  void set_mesh_med(ParaMEDMEM::MEDCouplingUMesh  & mg_mesh_med_in) {
    _med_mesh=&mg_mesh_med_in;
    return;
  }

  LibMeshFunction * get_interface_function_BC(boundary_id_type id);
//   FemusBCType getBCType(boundary_id_type id);
  void erase_interface_function_BC(boundary_id_type id);









  // ==============  _Libmesh ===============================================


  void terminate() {}

  void setMesh(const ParaMEDMEM::MEDCouplingUMesh * /*mesh*/) {
    return;
  }
  // ==========================================================================
  const ParaMEDMEM::MEDCouplingUMesh * getMesh() const {
    return _med_mesh;
  }
// ============================================================================
/// This function gives the identity to interfaces-meshes
/// by adding a new interface-function
//   int  defineIdInterface(                      ///< local id (return)
//     const boundary_id_type  bd_name_int,       ///< global name id (int)
//     const ParaMEDMEM::MEDCouplingUMesh * bc    ///< bc mesh
//   );


  void solve() {
    return;
  }

  void setAnalyticSource(const char * /*f*/) {
    return;
  };

  void setSource(
    const ParaMEDMEM::MEDCouplingUMesh * mesh,
    const ParaMEDMEM::MEDCouplingFieldDouble * f);
//    virtual void setSource(const ParaMEDMEM::MEDCouplingUMesh * mesh,const char *s);

//   virtual void setAverageSources(int n, double val[]);

  void setBoundaryConditionType(int, const char * typeBC);

  void setAnalyticBoundaryValues(int, const char * /*f*/) {
    return;
  };
  void setFieldBoundaryValues(
    int,
    int n_cmp,
    const ParaMEDMEM::MEDCouplingFieldDouble * f);

  ParaMEDMEM::MEDCouplingFieldDouble * getOutputField(const char */* vName*/) const {
    return NULL;
  }
  ParaMEDMEM::MEDCouplingFieldDouble * getBoundaryValues(const ParaMEDMEM::MEDCouplingUMesh * /*bc*/) const {
    return NULL;
  }
  
// ============================================================================
/// This function gets the  the value of the variable with id number
///  "variable_id" on nodes on the boundary with identity "id" in the
///  system "system_name"
  ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary_nodes(
    int id,                            // int boundary identity   (in)
    const char *system_name,           // system name             (in)
    int         variable_id            // variable identity       (in)
  );// ========================================================================
// ============================================================================
/// This function gets the all the values on nodes on the boundary
/// with identity "id" in the system "system_name"
  ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary_nodes(
    int id,                           // int boundary identity   (in)
    const char *system_name           // system name             (in)
  ) ;// ========================================================================
  
  
ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary(
  int id,
  std::vector<char *> vName
) ;
//   ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary_nodes(int, const char *);
//   ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary_nodes(int);
//   ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary(int,  const char *) const;


  // =============================================================






};

#endif
