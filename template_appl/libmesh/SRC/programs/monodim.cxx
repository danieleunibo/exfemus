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

// structure BCondition
// boundary can be analitic or field
struct sBC {
  std::string name;  // name
  std::string type;  // type
  bool isAnalytic;                // analitic
  std::string equation;           // equation
  MEDCouplingFieldDouble * field; //field
};

void getInput(LIBMESH &, std::vector<sBC> & boundary);

int main(int argc, char **argv){

  Timer T;
  bool result;

  fDebug.open("test2", MPI_COMM_NULL);

  std::string dataFile = "MESH/line.med";
//   double good = -0.0530494;

  if (argc > 1)
    dataFile = argv[1];
//   if (argc > 2)
//     good = std::strtod(argv[2], NULL);
  
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
  P.setType("monodim");
  fDebug << "Elapsed (problem type)       : " << T.elapsed() << std::endl;
    
  T.start();
  P.setMesh(dataFile);

  fDebug << "Elapsed (set mesh)           : " << T.elapsed() << std::endl;
  
  T.start();
  P.setAnalyticSource("10.");
  fDebug << "Elapsed (set source)         : " << T.elapsed() << std::endl;

    
  T.start();
  // get boundary condition <- getInput 
  std::vector<sBC> boundary;  getInput(P, boundary);
  
  for (int i=0; i<boundary.size(); i++) {
    sBC & b = boundary[i];
    // bc are 1) analitic 2) numeric ( MEDCouplingFieldDouble::field)
    if (b.isAnalytic) P.setAnalyticBoundaryValues(b.name, b.type, b.equation);
    else                 P.setFieldBoundaryValues(b.name, b.type, b.field);
  }

  fDebug << "Elapsed (set boundary cond.) : " << T.elapsed() << std::endl;

  
  
  
   // time loop
  for (int itime=0; itime<20; itime++) {

    // Solving --------------------------------------------------
    //   6. solving
    T.start();  
    P.solve();   
    fDebug << "\n solve at "<<  itime <<" : t= " << T.elapsed() << std::endl;
    // Postprocessing --------------------------------------------------
    // 7. solution storage -----------------------------------------
    T.start();
    sol = P.getOutputField("u");
//     fDebug << "\n 7. get output : t= " << T.elapsed() << std::endl;
    // 8. solution write
    T.start();
    std::string s = "RESU/monodim_sol_";  s += prefix;  s += ".med";
    std::cout << " File "<< s <<"\n";
    std::ostringstream  s2;
    s2<<"RESU/monodim_sol_" << prefix<<"."<<itime<<".med";
//         s2<<"RESU/monodim_sol_" << prefix<<".med";
    MEDLoader::WriteField(s2.str().c_str(), sol, false);

  } // end time
  
  P.terminate();
  fDebug << "Elapsed (terminate)          : " << T.elapsed() << std::endl;
  
  if (sol) {
//     unsigned int n = sol->getNumberOfTuples()/2-1;
//     result = fabs(sol->getIJ(n,0) - good) < 1e-5;
//       std::ostringstream s;
//       s << "Monodim2: result at node " << n 
// 	<< " = " << sol->getIJ(n,0) << " expected " << good;
//       s << std::endl;
//       std::cerr << s.str() << std::endl;

    sol->decrRef();
  }

  return result ? 0 : 1;
}

// =========================================================
// get input boundary condition
void getInput(LIBMESH &P, std::vector<sBC> & boundary)
{

  // Field
  std::string inFieldFile = "MESH/in_square2_2.med";

  std::vector<std::string> InFieldNames = MEDLoader::GetAllFieldNames(inFieldFile.c_str());
  std::vector<std::string> InMeshNames  = MEDLoader::GetMeshNames(inFieldFile.c_str());

  MEDCouplingFieldDouble *infield = 
     MEDLoader::ReadField(ON_NODES, inFieldFile.c_str(), InMeshNames[0].c_str(),
                          0, InFieldNames[0].c_str(), -1, -1);

  // boundary condition construction   
  std::vector<std::string> bcNames = P.getBoundaryNames();
  int nBC = bcNames.size();
  boundary.resize(nBC);
  
  for (int i=0; i<nBC; i++) {
    boundary[i].name = bcNames[i];
    boundary[i].type = "Neumann";
    boundary[i].isAnalytic = true;
  }


  return;
 }
 
