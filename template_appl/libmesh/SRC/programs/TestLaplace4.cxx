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

struct sBC {
  std::string name;
  std::string type;
  bool isAnalytic;
  std::string equation;
  MEDCouplingFieldDouble * field;
};

void getInput(LIBMESH &, std::vector<sBC> & boundary);

int main(int argc, char **argv)
{

  Timer T;
  bool result;

  fDebug.open("test4", MPI_COMM_NULL);

  std::string dataFile = "MESH/square2.med";
  double good = -0.0524589;

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
  
  T.start();
  P.setAnalyticSource("0");
  fDebug << "Elapsed (set source)         : " << T.elapsed() << std::endl;

    
  T.start();
  std::vector<sBC> boundary;
  getInput(P, boundary);
  for (int i=0; i<boundary.size(); i++) {
    sBC & b = boundary[i];
      
    if (b.isAnalytic)
      P.setAnalyticBoundaryValues(b.name, b.type, b.equation);
    else 
      P.setFieldBoundaryValues(b.name, b.type, b.field);
  }

  fDebug << "Elapsed (set boundary cond.) : " << T.elapsed() << std::endl;

  T.start();
  P.solve();
  fDebug << "Elapsed (solve)              : " << T.elapsed() << std::endl;
    
  T.start();
  sol = P.getOutputField("u");
  fDebug << "Elapsed (get output)         : " << T.elapsed() << std::endl;
  
  T.start();
  std::string s = "RESU/test4_sol_";
  s += prefix;
  s += ".med";
  MEDLoader::WriteField(s.c_str(), sol, true);

  for (int i=0; i<boundary.size(); i++) {
    bdy = P.getValuesOnBoundary(boundary[i].name, "u");
    std::ostringstream s;
    s << "RESU/test4_bdy_" << prefix << "_" << i << ".med";
    MEDLoader::WriteField(s.str().c_str(), bdy, true);
  }

  P.terminate();
  fDebug << "Elapsed (terminate)          : " << T.elapsed() << std::endl;
  
  if (sol) {
    unsigned int n = sol->getNumberOfTuples()/2-1;
    result = fabs(sol->getIJ(n,0) - good) < 1e-5;
      std::ostringstream s;
      s << "TestLaplace4: result at node " << n 
	<< " = " << sol->getIJ(n,0) << " expected " << good;
      s << std::endl;
      std::cerr << s.str() << std::endl;

    sol->decrRef();
  }

  return result ? 0 : 1;
}

void getInput(LIBMESH &P, std::vector<sBC> & boundary)
{

  std::string inFieldFile = "MESH/in_square2_2.med";

  std::vector<std::string> InFieldNames 
    = MEDLoader::GetAllFieldNames(inFieldFile.c_str());
  std::vector<std::string> InMeshNames 
    = MEDLoader::GetMeshNames(inFieldFile.c_str());

  MEDCouplingFieldDouble *infield = MEDLoader::ReadField
    (ON_NODES, inFieldFile.c_str(), InMeshNames[0].c_str(),
     0, InFieldNames[0].c_str(), -1, -1);

  if (infield->getNature() == NoNature)
    infield->setNature(ConservativeVolumic);

  std::vector<std::string> bcNames = P.getBoundaryNames();
  int nBC = bcNames.size();
  boundary.resize(nBC);
  
  for (int i=0; i<nBC; i++) {
    boundary[i].name = bcNames[i];
  }
  
  boundary[0].type = "Dirichlet";
  boundary[0].isAnalytic = true;
  boundary[0].equation = "0.25-y*y";
  boundary[0].field = NULL;

  boundary[1].type = "Neumann";
  boundary[1].isAnalytic = true;
 
  boundary[2].type = "Dirichlet";
  boundary[2].isAnalytic = false;
  boundary[2].field = infield;

  boundary[3].type = "Dirichlet";
  boundary[3].isAnalytic = true;
  boundary[3].equation = "0.0";
  boundary[3].field = NULL;
 }
