#include <stdexcept>
#include "LIBMESH.hxx"
#include <iostream>
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <set>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <mpi.h>
#include <libmesh/libmesh.h>
#include <cstdlib>
#include "Debug.hxx"
#include "Timer.hxx"

#define FDEBUG 0

using namespace ParaMEDMEM;

int main(int argc, char **argv)
{

  Timer T;
  bool result;

  fDebug.open("test3", MPI_COMM_NULL);

  MEDCouplingUMesh * mesh; 

  std::string dataFile = "MESH/Sphere.med";
  double good = 0.509397;

  if (argc > 1)
    dataFile = argv[1];
  if (argc > 2)
    good = std::strtod(argv[2], NULL);
  
  std::string prefix;
  size_t p = dataFile.find_last_of('/');
  size_t q = dataFile.find_last_of('.');
  p = (p == std::string::npos) ? 0 : p+1;
  prefix = dataFile.substr(p, q-p);

  MEDCouplingFieldDouble * sol = NULL, * bdy = NULL;
    
  T.start();
  LIBMESH P;
  fDebug << "Elapsed (initialize)         : " << T.elapsed() << std::endl;
    
  T.start();
  P.setType("Poisson");
  fDebug << "Elapsed (problem type)       : " << T.elapsed() << std::endl;
    
  T.start();
  P.setMesh(dataFile);
  fDebug << "Elapsed (set mesh)           : " << T.elapsed() << std::endl;
  
  std::vector<std::string> bcNames = P.getBoundaryNames();
  int nBC = bcNames.size();

  T.start();
  P.setAnalyticSource("0");
  fDebug << "Elapsed (set source)         : " << T.elapsed() << std::endl;

  T.start();
  std::string vBC;
  for (int i=0; i<nBC; i++) {
    if (i==0)
      vBC = "if(abs(z)<0.5,1.0-4*z*z,0.0)";
    else
      vBC = "0.0";
    P.setAnalyticBoundaryValues(bcNames[i], "Dirichlet", vBC);
  }
  fDebug << "Elapsed (set boundary cond.) : " << T.elapsed() << std::endl;

    
  T.start();
  P.solve();
  fDebug << "Elapsed (solve)              : " << T.elapsed() << std::endl;
    
  T.start();
  sol = P.getOutputField("u");
  fDebug << "Elapsed (get output)         : " << T.elapsed() << std::endl;
  
  T.start();
  std::string s = "RESU/test3_sol_";
  s += prefix;
  s += ".med";
  MEDLoader::WriteField(s.c_str(), sol, true);

  for (int i=0; i<nBC; i++) {
    bdy = P.getValuesOnBoundary(bcNames[i], "u");
    std::ostringstream s;
    s << "RESU/test3_bdy_" << prefix << "_" << i << ".med";
    MEDLoader::WriteField(s.str().c_str(), bdy, true);
  }

  P.terminate();
  fDebug << "Elapsed (terminate)          : " << T.elapsed() << std::endl;
  
  if (sol) {
    unsigned int n = 2*sol->getNumberOfTuples()/3-1;
    double v = sol->getIJ(n,0);
    result = fabs(v - good) < 1e-5;
 
    std::ostringstream s;
    s << "TestLaplace3: result at node " << n << " ) = "
      << sol->getIJ(n,0) << " expected " << good;
    std::cerr << s.str() << std::endl;

    sol->decrRef();
  }
  return result ? 0 : 1;
}
