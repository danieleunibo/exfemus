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
#include <vector>
#include <mpi.h>
#include <libmesh/libmesh.h>
#include <cstdlib>
#include "Debug.hxx"
#include "Timer.hxx"

#define FDEBUG 0

using namespace ParaMEDMEM;

void fileNameSplit(int rank,
		   std::string & fileName,
		   std::string & baseName,
		   std::string & dirName,
		   std::string & suffix);


int main(int argc, char **argv)
{
 
  bool result;
  MPI_Init (&argc, &argv);

  Timer T;

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  fDebug.open("testp1_", MPI_COMM_WORLD);

  MEDCouplingUMesh * mesh; 

  std::string dataFile = "Sphere";
  std::vector<double> good(size);
  good[0] = 0.0186867;
  good[1] = 0.0600085;
  good[2] = 0.474626;
  good[3] = 0.343559;
  for (int i=4; i<size; i++) good[i] = -1;

  if (argc > 1)
    dataFile = argv[1];
  if (argc > 5) {
    good[rank] = std::strtod(argv[rank+2], NULL);
  }
  std::string dirName, baseName, suffix;
  fileNameSplit(rank, dataFile, baseName, dirName, suffix);
  
  MEDCouplingFieldDouble * sol = NULL;


  {
    T.start();
    LIBMESH P(MPI_COMM_WORLD);
    fDebug << "Elapsed (initialize)         : " << T.elapsed() << std::endl;
     
    T.start();
    P.setType("Poisson");
    fDebug << "Elapsed (problem type)       : " << T.elapsed() << std::endl;

    T.start();
    try {
      P.setMesh(dataFile);
    }
    catch (...) {
      MPI_Finalize();
      return 1;
    }
    fDebug << "Elapsed (set mesh)           : " << T.elapsed() << std::endl;

    T.start();
    P.setAnalyticSource("0");
    fDebug << "Elapsed (set source)         : " << T.elapsed() << std::endl;
    
    T.start();
    std::vector<std::string> bcNames = P.getBoundaryNames();
    int nBC = bcNames.size();
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
    
    std::ostringstream s;
    s << "testp1_sol_" << baseName << "_" << rank+1 << suffix;
    MEDLoader::WriteField(s.str().c_str(), sol, true);
    fDebug << "Elapsed (write solution)     : " << T.elapsed() << std::endl;

    T.start();
    P.terminate();
    fDebug << "Elapsed (terminate)          : " << T.elapsed() << std::endl;
        
    if (sol) {
      unsigned int n = 2*sol->getNumberOfTuples()/3-1;
      result = good[rank] >= 0 
	? fabs(sol->getIJ(n,0) - good[rank]) < 1e-5 
	: true;
      std::ostringstream s;
      s << "TestPLaplace1: result at node " << n << " ) = " << sol->getIJ(n,0);
      if (good[rank] < 0) 
	s << " unknown";
      else
	s << " expected " << good[rank];
      s << std::endl;
      std::cerr << s.str() << std::endl;
      
      sol->decrRef();
    }
  }   
  MPI_Finalize();
  
  std::cerr << "end" << std::endl;
  return result ? 0 : 1;
} 


void fileNameSplit(int rank, 
		   std::string & fileName,
		   std::string & baseName,
		   std::string & dirName,
		   std::string & suffix)
{
  size_t p  = fileName.find_last_of('/');
  p = (p == std::string::npos) ? 0 : p+1;
 
  size_t q  = fileName.find_last_of('.');

  if (q == std::string::npos)
    baseName = fileName.substr(p);
  else {
    suffix = fileName.substr(q, std::string::npos);
    baseName = fileName.substr(p, q-p);
  }

  dirName = fileName.substr(0,p);

/*  
  std::ostringstream s;
  s << dirName << baseName << "_" << rank+1 << suffix;
  fileName = s.str();
*/
}
