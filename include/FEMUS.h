#ifndef __FEMUS__
#define __FEMUS__

#include <vector>
#include "Solverlib_conf.h"

#ifdef HAVE_MED
namespace ParaMEDMEM {
class MEDCouplingUMesh;
class MEDCouplingFieldDouble;
}
#endif

class MeshExtended;
class MGFemusInit;
class EquationSystemsExtendedM;
class MGSystem;
class MGUtils;
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
//   int                                    /*_num_mgmesh*/;   // num of MGmeshes
  MeshExtended*                          _mg_mesh;      // FEMus-mesh
#ifdef HAVE_MED
  ParaMEDMEM::MEDCouplingUMesh *         _med_mesh;     // Med-mesh
#endif
  // system data ---------------------------------------------------------------------------
  //   _LibMeshProblem * _problem;
  EquationSystemsExtendedM*           _mg_equations_map; // system
  MGTimeLoop *                            _mg_time_loop;     // transient


  // ==========================================================================
  //=========  public functions  ============================================
public:

  // Constructor-Destructor
  FEMUS();
  FEMUS(MPI_Comm comm);
  void init_param(MGUtils &mgutils);
  void init_fem(MGGeomEl & mggeomel,MGFEMap & mgfemap);
  ~FEMUS();
//   void terminate();

  void setMesh();
  void setSystem(const std::vector<FIELDS> & pbName,
      int n_data_points=0,
  int n_data_cell=0
  );

 //====================================================================
/// This function sets up the intial set
  void solve_setup(
    int        & t_in,                 ///< initial time iteration
    double     &  time                 ///< actual time
  );
  //=============================================================================
// This function solves the problem
//   void solve();

//=============================================================================
// This function solves one step  for transient problems
  void solve_onestep(
    const int  & t_in,                 ///< initial time iteration
    const int  & t_step,               ///< actual time iteration
    const int  & print_step,            ///< print every
    double     &  time,                ///< actual time
    double     &  dt                   ///< step time
  ) ;
  
  
   const  MeshExtended & get_MGMesh(){return *_mg_mesh;};
   const EquationSystemsExtendedM& get_MGExtSystem(){return *_mg_equations_map;};
#ifdef HAVE_MED
  void setMedMesh(const std::string & dataFile);
// // ========================================================================
// This function gets the Group mesh from the med file (medfile_name)
// and sets the id and the mesh to the interfaces function
// through the index of the interface-functions
// (from EquationSystemsExtendedM)
  void init_interface(
    const int interface_name,
    const int interface_id,
    int order_cmp,
    const std::string & medfile_name, // medfile name    (in)
    bool on_nodes=true,
    const int index_medmesh=0                      // med-mesh index  (in)

  ) ;

  void init_interface(
    const int interface_name,
    int interface_id,
    int order_cmp,
    const std::string & medfile_name, // medfile name    (in)
    const std::string & medfile_name_old, // medfile name    (in)
    const FEMUS & P_old,
    const int index_mgmesh,           // mgmesh index    (in)
    const int index_medmesh=0                      // med-mesh index  (in)
  ) ;


  // =========================================================================
  /// This function sets the value from first_cmp to end_cmp of the field
  /// on the old solution x_old of the interface id
  void write_Boundary_value(
    int id_boundary_name   ,    ///< identity interface name (in)
    std::string mgsystem_name, ///< system name          (in)
    int n_cmp,             ///< from variable system (in)
    int first_cmp =0             ///< to variable system   (in)

  );


// ============================================================================
  /// This function sets the field of interface function (interface_id)
  /// with an analitic expression
  void setAnalyticSource(
    int interface_name,
    int n_cmp,
    const std::string & bcExpression      // boundary symbolic expr
  );
  /// This function sets the field of interface function (interface_id)
  void setFieldSource(
    int interface_name,
    int n_cmp,
    const ParaMEDMEM::MEDCouplingFieldDouble * srcField);
  /// This function sets the field of interface function (interface_id)
  /// with an analitic expression
  void setAnalyticBoundaryValues(
    int interface_name,
    int n_cmp,
    const std::string & f);
  /// This function sets the field of interface function (interface_id)
  void setFieldBoundaryValues(
    int interface_name,
    int n_cmp,
    const ParaMEDMEM::MEDCouplingFieldDouble * bcField);

 
// // ============================================================================
/// This function gets all the values on boundary with identity id
  ParaMEDMEM::MEDCouplingFieldDouble * getValuesOnBoundary(  // field (out)
    int interface_name,                  ///< boundary name (char*) (in)
    const std::string & systemName,              ///< system name           (in)
    int n_cmp,                                    ///< component             (in)
    int first_cmp=0                                     ///< component             (in)
  );
  
#endif



};

#endif