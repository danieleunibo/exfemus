#ifndef __FEMBCTYPES_HXX__
#define __FEMBCTYPES_HXX__

#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingFieldDouble.hxx"

// enum FemusBCType {
//   Undefined = -1,
//   UnMark = 0,
//   Mark = 1
// //   Neumann,
// //   Displacement,
// //   Force
// };

struct BCdata {
  std::string name;
  std::string type;
  bool isAnalytic;
  std::string equation;
  ParaMEDMEM::MEDCouplingFieldDouble * field;
};

#endif
