#ifndef __LIBMESH_ANALYTICFUNCTION__
#define __LIBMESH_ANALYTICFUNCTION__

#include <iostream>
#include <map>
#include <vector>

// #include "libmesh/mesh_base.h"
// #include "libmesh/mesh.h"

#include "InterpKernelExprParser.hxx"
#include "MEDCouplingRefCountObject.hxx"




namespace ParaMEDMEM {
class MEDCouplingFieldDouble;
class MEDCouplingUMesh;
}
class ParallelMeshExtended;

class LibMeshFunction {

protected:         // data interface function
  // meshes
  const MGMesh                       * _mesh_femus;   // MGMesh
  const ParaMEDMEM::MEDCouplingUMesh * _support_med;  // MEDCoupling mesh
  // field
  ParaMEDMEM::MEDCouplingFieldDouble *_field;         // field
  // mapping
  std::map<int, int>                 _nodeID;// node map MGMesh-MEDCoupling
  std::map<std::pair<int, int>, int> _faceID;// face map MGMesh-MEDCoupling
  std::map<int, int>                 _elemID;// elem map MGMesh-MEDCoupling
  
public:
  // Constructors- Destructors --------------------------------------
  LibMeshFunction() :   ///< Constructor ============================
     _mesh_femus(NULL),       // Femus-mesh
    _support_med(NULL),     // med-sub/mesh 
    _field(NULL)           // function
 
    {         }
  ~LibMeshFunction();       ///< Destructor =========================
  
  // evaluation functions -------------------------------------------
  void eval(            ///< Evaluation function ====================
    int id,                  // index     
    std::vector<double> & v);// field vector    
  void eval(             ///< Evaluation function ===================
    int id,                  // index <-
    int side,                // side  <-
    std::vector<double> & v);// field vector ->
  void eval_elem(
  int iel,
  std::vector<double> & val);
  
  void eval_node(
  int node,
  std::vector<double> & val
) ;
  
  // set functions ------------------------------------------------------------
// ============================================================================
  void set_mesh_interface( ///< Setting the interface-meshes 
    const bool nodes_only,                        ///< nodes only flag 
    const MGMesh * mesh,                          ///< MGMesh     
    const ParaMEDMEM::MEDCouplingUMesh * support  ///< MedCouplng <-
  );
  
// ============================================================================ 
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
  void set_mesh_interface_nodeID(
   const int search_only_bc,                      ///< only bd needed (1)  (in) 
  const ParallelMeshExtended * mesh,             ///< Femus-mesh           (in)
  const ParaMEDMEM::MEDCouplingUMesh * support   ///< med-mesh             (in)
);
  
  
// ============================================================================ 
/// This function computes the mesh interface
/// to set in the interface-function. The field in the
/// interface-function is set by set_mesh_femus_interface
  void set_mesh_interface_nodeID(
  const ParallelMeshExtended * mesh,             ///< Femus-mesh           (in)
  const ParaMEDMEM::MEDCouplingUMesh * support,  ///< med-mesh             (in)
  const int search_only_linear                   ///< linear only flag     (in)
); 
 
  
  
  
  void set_analytic_field_interface(///< Setting the field (analytic)=
    const char *code,          // symbolic expression  <-
    int nComp);                // n of components      <-
  void setField(    ///< Setting the field (MEDfield)=================
    const ParaMEDMEM::MEDCouplingFieldDouble *f); // field  <-
   void set_field_source(    ///< Setting the field (MEDfield)=================
    const ParaMEDMEM::MEDCouplingFieldDouble *f); // field  <-
  
  // get functions ---------------------------------------------------
  std::map<std::pair<int, int>, int> & FaceID() { return _faceID; }//=
  std::map<int, int> & NodeID() { return _nodeID; } // ===============
  std::map<int, int> & ElemID() {return _elemID;}   //================
  void getSupportBarycenterCoords( // Get barycenter MEDmesh =========
    int dim_med,                           // dimension         <-
    int n_elements_med,                    // n of elements     <-
    const ParaMEDMEM::DataArrayDouble * d, // point coords      <-
    double *xyz_m);                        // barycenter vector ->
  void getMGMeshBarycenterCoords(// Get barycenter MGMesh    =========
    int dim,                 // dimension         <-
    int nlevels,             // mesh level          <-
    int n_nodes,             // n of elements     <-
    int n_elements,          // n of nodes     <-
    double *xyz_m);          // barycenter vector ->
  ParaMEDMEM::MEDCouplingFieldDouble * getField(
    const char * vName);
  const ParaMEDMEM::MEDCouplingUMesh * getSupport(
  ) {return _support_med;}

  void printOn(std::ostream & out) const;


};

#endif

