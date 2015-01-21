#ifndef __LIBMESH_UTIL__
#define __LIBMESH_UTIL__

// #include "libmesh/node.h"
#include <set>
#include <vector>

namespace ParaMEDMEM {
  class DataArrayDouble;
}

// class NodeExt : public libMesh::Node {
// 
// public:
//   NodeExt(): Node(0.0,0.0,0.0,0), _id2(0) {};
//   NodeExt(const libMesh::Point &n, int id, int id2): 
//     Node(n(0), n(1), n(2), id), _id2(id2) {
//   }
//   int & set_id2() { return _id2; }
//   int id2() const { return _id2; }
// protected:
//   int _id2;
// private:
//   //  NodeExt(const NodeExt &q) : Node(q), _id2(q._id2) {}
// };

// struct nodeCompare {
//   
//   static const double eps;
//   bool operator() (const libMesh::Node * p, const libMesh::Node * q) const {
//     
//     unsigned i;
//     double d;
// 
//     for (i=0; i<LIBMESH_DIM; i++) {
//       d = (*p)(i) - (*q)(i);
//       if (fabs(d) > eps) break;
//     }
//   
//     return i<LIBMESH_DIM && d < - eps;
//   }
// };
// 
// void commonPoints(std::map<std::pair<int, int>, int > & corresp,
// 		  const std::set<const libMesh::Node *, nodeCompare> & p, 
// 		  const std::set<const NodeExt *, nodeCompare> & q);
// 
// void commonPoints(std::map<int, int> & corresp,
// 		  const std::set<const libMesh::Node *, nodeCompare> & p, 
// 		  const std::set<const libMesh::Node *, nodeCompare> & q);
// 
// void commonPoints(std::map<int, int> & corresp,
// 		  const ParaMEDMEM::DataArrayDouble * d1, 
// 		  const ParaMEDMEM::DataArrayDouble * d2);
// 
// const libMesh::Node * findNode
//   (const std::set<const libMesh::Node *, nodeCompare> & p, int id);
// 
// 
// void DataArrayToSetOfNodes(const ParaMEDMEM::DataArrayDouble * d, 
//  			   std::set<const libMesh::Node *, nodeCompare> & s,
//  			   std::vector<NodeExt> & n);
// 
// template <typename T>
// void CleanSetOfNodes(std::set<const T *, nodeCompare> & s)
// {
//   typename std::set<const T *, nodeCompare>::iterator itS;
//   for (itS = s.begin(); itS != s.end(); itS++)
//     delete *itS;
//   s.clear();
// }


#endif
