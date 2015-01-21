#ifndef __mgsosadj0_h__
#define __mgsosadj0_h__

#include "Equations_conf.h"

#ifdef NSAD_EQUATIONS

// Local Includes
#include "MGSolverDA.h"
// Forwarded classes
class MGSystem;

class MGSolNSAD : public MGSolDA {

  public:

  MGSolNSAD(const MGSystem & mgsys,
	  const uint NoLevels_in,
          std::string eqname_in="NSAD",
          std::string varname_in="lambda");

  ~MGSolNSAD(){};

  void ic_read(double xp[],double u_value[]);

  void bc_read(double xp[],double normal[],int u_value[]);

  void GenMatRhs(const double /*time*/,const uint Level,const int mode);

  void GenRhsB(const double time, const uint Level);


};


#endif

#endif
