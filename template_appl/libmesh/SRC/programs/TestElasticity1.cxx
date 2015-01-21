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
#include "Debug.hxx"
#include <mpi.h>
#include <libmesh/libmesh.h>
#include <cstdlib>

using namespace ParaMEDMEM;

int main(int argc, char **argv)
{
 
  bool result;

  fDebug.open("test1E", MPI_COMM_NULL);
  MEDCouplingUMesh * mesh, *mBC; 

  std::string dataFile = "MESH/Sphere.med";
  double good = 0.124878;

  if (argc > 1)
    dataFile = argv[1];
  if (argc > 2)
    good = std::strtod(argv[2], NULL);
  
  std::string prefix;
  size_t p = dataFile.find_last_of('/');
  size_t q = dataFile.find_last_of('.');
  p = (p == std::string::npos) ? 0 : p+1;
  prefix = dataFile.substr(p, q-p);
  
  MEDCouplingFieldDouble * sol = NULL;
  
  LIBMESH P;
  P.setType("Elasticity");
  P.setMesh(dataFile);

  std::vector<std::string> bcNames = P.getBoundaryNames();

  {
    int i;
    for (i=0; i<bcNames.size(); i++) {
      if (bcNames[i] == "In")
	P.setAnalyticBoundaryValues(bcNames[i], "Force", "5000*z*KVec");
    }
  }
  
  P.solve();
  sol = P.getOutputField("u");
  
  std::string s = "RESU/test1E_sol_";
  s += prefix;
  s += ".med";
  MEDLoader::WriteField(s.c_str(), sol, true);

  P.terminate();

  if (sol) {
    // unsigned int n = sol->getNumberOfTuples()/2-1;
    // result = fabs(sol->getIJ(n,0) - good) < 1e-5;
    // std::ostringstream s;
    // s << "TestLaplace1: result at node " << n << " = ";
    // s << sol->getIJ(n,0) << " expected " << good;
    // s << std::endl;
    // std::cerr << s.str() << std::endl;

    std::vector<const MEDCouplingFieldDouble *> vsol(1);
    vsol[0] = sol;
    MEDCouplingFieldDouble::WriteVTK("RESU/test1E_sol.vtk", vsol);

    sol->decrRef();
  }
  result = 1;
  return result ? 0 : 1;
}
