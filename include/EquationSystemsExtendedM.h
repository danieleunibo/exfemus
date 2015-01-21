#ifndef _EQNEXTENDED__
#define _EQNEXTENDED__

#include <map>

#include "MGEquationsSystem.h"
#include "Solverlib_conf.h"

#ifdef HAVE_MED
namespace ParaMEDMEM {
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
}

class InterfaceFunctionM;
#endif

class MeshExtended;


// ==========================================================
// ==========================================================

class EquationSystemsExtendedM : public MGEquationsSystem {

// data ------------------------------------------------------
protected:
  // mesh data
  const MeshExtended                 * _mg_mesh;  ///< MG-mesh
  
#ifdef HAVE_MED
  const ParaMEDMEM::MEDCouplingUMesh * _med_mesh; ///< MED-mesh
  // system
  std::map<int,InterfaceFunctionM *> _interfaceFunMap;///< MG-Med interface map
#endif

public:
// ============================================================================
// Constructor-Desctructor
// ============================================================================
  EquationSystemsExtendedM(   ///< Constructor function
    MGUtils& mgutils_in,             ///<  MGUtils class       (in)
    MeshExtended& mgmesh_in,         ///<  MeshExtended class  (in)
    MGFEMap& mgfemap_in,             ///<  MGFEMap class       (in)
    int npoint_data,                 ///<  point data          (in)
    int ncell_data                   ///<  data cell           (in)
  );
// ============================================================================
  virtual ~EquationSystemsExtendedM();     ///< Desctructor function

// ============================================================================
  void set_mesh_mg(
    MeshExtended & mg_mesh_femus_in
  ) { _mg_mesh=&mg_mesh_femus_in; return;}
  const MeshExtended * getMeshMG() const {return _mg_mesh;}
// ==================  Mesh ===================================================
// ============================================================================

#ifdef HAVE_MED
  const ParaMEDMEM::MEDCouplingUMesh * getMeshMED() const {return _med_mesh;}
    void set_mesh_med(
    ParaMEDMEM::MEDCouplingUMesh  & mg_mesh_med_in
  ) {_med_mesh=&mg_mesh_med_in; return;}





// ============================================================================
// ==================  interfaceFunMap ========================================
// ============================================================================
/// This function returns the pointer of the  interface with identity id

  InterfaceFunctionM * get_interface_fun(
    int id                                   ///< interface identity
  ) {
    return _interfaceFunMap.find(id) ==
           _interfaceFunMap.end() ? NULL :_interfaceFunMap[id];
  }

// ============================================================================
  /// This function erases the interface with identity id
  void erase_interface_fun(
    int  id                                  ///< interface identity
  )  { _interfaceFunMap.erase(id); return;}
//=============================================================================
  ///  This functon add an interface-function (med-mesh,femus-mesh, field)
  ///  in the interface-function boundary map (_interfaceFunMap):
  ///  the field in the interface-function is not assigned here
  void  add_interface_fun(                      // map position            (return)
    const int  interface_name,
    const int  bd_name_int,    ///< boundary id (int)     (in)
    const ParaMEDMEM::MEDCouplingUMesh * b, ///< med-submesh           (in)
    const bool on_nodes,                     ///< values on nodes (true)
    const int order_cmp=2                    ///< order component (2=quad;lin=1)
  );


// ============================================================================
// ==================  functions===============================================
// ============================================================================
// ============================================================================
  void setBC(
    int name,
    int n_cmp,
    const char *s);
  void setBC(
    int name,
    int n_cmp,
    const ParaMEDMEM::MEDCouplingFieldDouble *f);

// =========================================================================
/// This function sets the value from first_cmp to end_cmp of the field
/// on the old solution x_old of the interface id
  // =========================================================================
  void write_Boundary_value(
    int id_name,                ///< identity interface name (in)
    std::string mgsystem_name, ///< system name          (in)
    int n_cmp,                 ///< from variable system (in)
    int first_cmp=0            ///< to variable system   (in)
  );

// ============================================================================
// This function gets the  the value of the variable with id number
//  "variable_id" on nodes on the boundary with identity "id" from
//  variable_start to n_variable in the
//  system "system_name"
  ParaMEDMEM::MEDCouplingFieldDouble *getValuesOnBoundary_nodes(
    int id,                            // int boundary identity   (in)
    const char *system_name,           // system name             (in)
    int         n_cmp,            //  first variable       (in)
    int         first_cmp=0            // n variables       (in)
  );
  
// ============================================================================
// This function gets the  the value of the variable with id number
//  "variable_id" on elements on the boundary with identity "id" from
//  variable_start to n_variable in the
//  system "system_name"
  ParaMEDMEM::MEDCouplingFieldDouble *getValuesOnBoundary_elem(
    int id,                            // int boundary identity   (in)
    const char *system_name,           // system name             (in)
    int         n_cmp,            //  first variable       (in)
    int         first_cmp=0            // n variables       (in)
  );
#endif

};

#endif
