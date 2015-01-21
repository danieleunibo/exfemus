#ifndef __FEMUSee_UTIL__
#define __FEMUSee_UTIL__

// #include "libmesh/node.h"
#include <set>
#include <vector>
#include <memory>
#include "Domain_conf.h"
#include <math.h>
#include <map>
#include "type_vector.h"
namespace ParaMEDMEM {
  class DataArrayDouble;
}
class MGMesh;



// class Point {
//  public:
// 
//   /// Constructor.  By default sets all entries to 0.  Gives the point 0 in
//   Point(const double x=0., const double y=0., const double z=0){
//    _coords[0]=x;
//    _coords[1]=y;
//     _coords[3]=z;
//   }
//   /// Copy-constructor.
// //   Point (const Point& p);
// //   /// Copy-constructor.
// //   Point (const TypeVector<double>& p);
//   /// Empty.
//   ~Point() {}
// 
// //   /**
// //    * @returns a key associated with this point.  Useful for sorting.
// //    */
// //   int key() const;
// 
// 
// const double  & operator () (const unsigned int i) const{
//   return _coords[i];
// }
// 
//  double _coords[DIMENSION];
//  protected:
// 
//   /**
//    * Make the derived class a friend
//    */
//   friend class Node;
// };



 // ======================================================
class Node : public Point{

public:
  int  _id;
 /// ==============================
  //  Constructor.
  explicit  Node  (const double x=0,
	 const double y=0,
	 const double z=0,
	 const int id = -1);
  //  Copy-constructor.
  Node (const Node& n);
  /// Copy-constructor from a \p Point.  Optionally assigned the \p id.
  explicit Node (const Point& p,
	         const int id = -1);
  /// Destructor.
  ~Node ();
  /// \returns the \p id for this \p DofObject
  int id () const{ return _id;}
  /// \returns the \p id for this \p DofObject as a writeable reference.
  int & set_id (){  return _id;}  
  /// Sets the \p id for this \p DofObject
  void set_id (const int dofid) { this->set_id() = dofid; }
  
  /// Assign to a node from a point
  Node& operator= (const Point& p);

  /// Builds a \p Node and returns an \p std::auto_ptr<Node> 
  static std::auto_ptr<Node> build (const Node& n);

  /// Builds a \p Node from \p Point p and returns an \p std::auto_ptr<Node>
  static std::auto_ptr<Node> build (const Point& p,const int id);

  /// Builds a \p Node from specified points and returns an \p std::auto_ptr<Node>
  static std::auto_ptr<Node> build (const double x,
			      const double y,
			      const double z,
			      const int id);

  /// @returns \p true if the node is active.  An active node is
  // defined as one for which \p id() is not \p Node::invalid_id.
  bool active () const;
  /// @returns \p true if this node equals rhs, false otherwise.
  bool operator ==(const Node& rhs) const;
  /// Prints relevant information about the node
  void print_info (std::ostream& os=std::cout) const;
//   void set_id(int id_in){id=id_in;}
  /// Prints relevant information about the node to a string.
  std::string get_info () const;
  
};

// ------------------------------------------------------------
// global Node functions
inline std::ostream& operator << (std::ostream& os, const Node& n){
  n.print_info(os);  return os;
}
//------------------------------------------------------
// ========================================================
// Inline functions
inline Node::Node (const double x,
	    const double y,
	    const double z,
	    const int dofid) :
  Point(x,y,z){
  this->set_id() = dofid;
}
// =========================================
inline Node::Node (const Node& n) :
  Point(n)
//   ,
//   DofObject(n),
//   ReferenceCountedObject<Node>()
  {
}
// ========================================
inline Node::Node (
  const Point& p,
  const int dofid) :
  Point(p)
  {
  // optionally assign the id.  We have
  // to do it like this otherwise
  // Node n = Point p would erase
  // the id!
  if (dofid != -1)  this->set_id() = dofid;
}
//==============================================
inline Node::~Node (){}
// ==============================================
inline Node & Node::operator= (const Point& p){
  (*this)(0) = p(0);
#if DIMENSION > 1
  (*this)(1) = p(1);
#endif
#if DIMENSION > 2
  (*this)(2) = p(2);
#endif

  return *this;
}
// ==============================================
inline std::auto_ptr<Node> Node::build(const Node& n){
  std::auto_ptr<Node> ap(new Node(n));
  return ap;
}
inline std::auto_ptr<Node> Node::build
      (const Point& p, const int id){
  std::auto_ptr<Node> ap(new Node(p,id));
  return ap;
}
// ==============================================
inline std::auto_ptr<Node> Node::build(const double x,
			  const double y,
			  const double z,
			  const int id){
  std::auto_ptr<Node> ap(new Node(x,y,z,id));
  return ap;
}
// ====================================
inline bool Node::active () const{
  return (this->id() != -1);
}


 
// class Node:public Point{
//   public:
//   int _id;
//   
//   Node (const double x,
//              const double y,
//             const double z,
//             const int dofid) :
//   Point(x,y,z){ set_id(dofid);}
//  
//  void set_id(int id_in){_id=id_in;}
//  int id()const{return _id;} 
// };
//   
  
  
// }
// 
// class 
// Node::Node (const double x,
//             const double y,
//             const double z,
//             const int dofid) :
//   Point(x,y,z)
// {
//   this->set_id() = dofid;
// }


// Extended nodes ==============================================
class NodeExt : public Node {
// 
 public:
   NodeExt(): Node(0.0,0.0,0.0,0), _id2(0) { }
   NodeExt(const Point &n, int id, int id2): 
     Node(n._coords[0], n._coords[1], n._coords[2], id), _id2(id2) {
   }
   int & set_id2() { return _id2; }
   int id2() const { return _id2; }
// protected:
   int _id2;
//  private:
  NodeExt(const NodeExt &q) : Node(q), _id2(q._id2) {}
 };
// ---------------------------------------------------
// Structure node compare
 struct nodeCompare {
  //  data   
  static const double eps;
  // function
  bool operator() (const Node * p, const Node * q) const { 
     unsigned i;  double d;
     // distance(p,q) in C_infinity 
     for (i=0; i<DIMENSION; i++) {
       d = (*p)(i) - (*q)(i);
       if (fabs(d) > eps) break;
     }
     
     return (i<DIMENSION && d < - eps);
  }
 };
// ----------------- end structure nodeCompare -------------------- 
// ===================================================================
void commonPoints(std::map<std::pair<int, int>, int > & corresp,
		  const std::set<const Node *, nodeCompare> & p, 
		  const std::set<const NodeExt *, nodeCompare> & q);
// ===================================================================
void commonPoints(std::map<int, int> & corresp,
		  const std::set<const Node *, nodeCompare> & p, 
		  const std::set<const Node *, nodeCompare> & q);
// ===================================================================
void commonPoints(std::map<int, int> & corresp,
		  const ParaMEDMEM::DataArrayDouble * d1, 
		  const ParaMEDMEM::DataArrayDouble * d2);
// ===================================================================
const Node * findNode
  (const std::set<const Node *, nodeCompare> & p, int id);

// ===================================================================
// This function sets   
void DataArrayToSetOfNodes(
  const ParaMEDMEM::DataArrayDouble * d,    // coords   <-
  std::set<const Node *, nodeCompare> & s,  // (Node map, compare) ->
  std::vector<NodeExt> & n                  // ext nod             -> 
);
// ==================================================================
void xyzToSetOfNodes(
  const int dim,
  const int ndim,
  const double  xyz[],                       // coords   <-
  std::set<const Node *, nodeCompare> & s, //,   // (Node map, compare) ->
  std::vector<NodeExt> & extNode_vect        // extNode vector   -> 
);

// ===========================================================================
/// This function set the boundary xyz coord vector into ExtNode vector n
  void xyzToSetOfNodes_onBd(
  const int dim,                           ///< dimension                 (in) 
  const int n_nodes,                       ///< number of nodes           (in)
  const double  xyz[],                     ///< coord vector              (in)
  const std::vector<int> bc_flag,          ///< boundary flag >9          (in)
//   int n_bd_nodes,                    ///< number of  bd nodes found (out)
  std::set<const Node *, nodeCompare> & s, ///< (Node map, compare)       (out)
  std::vector<NodeExt> & extNode_vect      ///< extNode vector            (out) 
);
  
// ============================================================================
/// This function set the boundary xyz coord vector into ExtNode vector n
void xyzToSetOfNodes_onBd_lin(
  const int dim,                           ///< dimension                 (in) 
  const int n_nodes,                       ///< number of nodes           (in)
  const double  xyz[],                     ///< coord vector              (in)
  const std::vector<int> bc_flag,          ///< boundary flag >9          (in)
//   int n_bd_nodes,                    ///< number of  bd nodes found (out)
  std::set<const Node *, nodeCompare> & s, ///< (Node map, compare)       (out)
  std::vector<NodeExt> & extNode_vect      ///< extNode vector            (out) 
);
  
// xyzToSetOfNodes_onbd(dim,n_nodes,mesh->_xyz,mesh->_bc_id, s2,n2);

// ===================================================================
// This function set the xyz vector into ExtNode victor n
void centerToSetOfNodes(
  const int dim_med,
  const int n_element_med,
  const double  xyz_m[],                 // coords   <-
  std::set<const Node *, nodeCompare> & s,   // (Node map, compare) ->{
  std::vector<NodeExt> & extNode_vect               // extNode vector   ->  
);
// ===================================================================
// This function set the xyz vector into ExtNode victor n
void centerToSetOfNodes2(
  const MGMesh &mesh_femus,
  const int dim_med,
  const int n_element_med,
   const int init_lev_el,
  const double  xyz_m[],                 // coords   <-
  std::set<const Node *, nodeCompare> & s,   // (Node map, compare) ->{
   std::vector<NodeExt> & extNode_vect               // extNode vector   ->  
);
void CleanSetOfNodes(std::set<const Node *, nodeCompare> & s);

#endif
