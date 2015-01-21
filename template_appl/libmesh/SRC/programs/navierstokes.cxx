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
#include "_LibMeshProblem.hxx"
using namespace ParaMEDMEM;

struct sBC {
  std::string name;
  std::string type;
  bool isAnalytic;
  std::string equation;
  MEDCouplingFieldDouble * field;
//   MEDCouplingFieldDouble * field1;
};

// user boundary condition function
void getInput(LIBMESH &, std::vector<sBC> & boundary);


// ==================================================
int main(
  int argc,
  char **argv
) { // ==============================================

  Timer T;         // time counter
  bool result;     // logic variable
  fDebug.open("navierstokes", MPI_COMM_NULL); // file log -> test2.log

  // Reading data file ----------------------------------------------
  std::string dataFile = "MESH/squarequad.med";// Initial dataFile
//   double good = -0.0530494;                    // Initial result check data
//   if(argc > 1)  dataFile = argv[1];            // dataFile from argv[1]
//   if(argc > 2)  good = std::strtod(argv[2], NULL); // Check data  from argv[2]

  // only filename ("square2") ------------------------------
  std::string prefix;
  size_t p = dataFile.find_last_of('/');
  size_t q = dataFile.find_last_of('.');
  p = (p == std::string::npos) ? 0 : p+1;
  prefix = dataFile.substr(p, q-p);
  // --------------------------------------------------------
  // coupling fields sol bdy  -------------------------------
  MEDCouplingFieldDouble *sol = NULL, *solv = NULL ,*solp = NULL , *bdy = NULL;//,*bdyv = NULL;//
  // **********************************************************
  // Preprocessing
  // 1- initialize -> LIBMESH P
  T.start();
  LIBMESH P;
  fDebug << "\n 1- initialize: t= " << T.elapsed() << std::endl;
  // **********************************************************
  // 2. set problem  -> P.setType
  T.start();
  P.setType("navierstokes");
  fDebug << "\n 2- set problem: t= " << T.elapsed() << std::endl;
  // **********************************************************
  // 3. set mesh  -> P.setType
  T.start();
  P.setMesh(dataFile);
  fDebug << "\n 3- set mesh: t= " << T.elapsed() << std::endl;
  // **********************************************************
  // 4. set source  ->  P.setAnalyticSource("0")
  T.start();
  P.setAnalyticSource("IVec*0.0 + JVec*0.0");
  fDebug << "\n 4- set source: t= " << T.elapsed() << std::endl;

  // **********************************************************
  // 5. boundary cds.  ->  P.setAnalyticBoundaryValues or P.setFieldBoundaryValues
  T.start();
  std::vector<sBC> boundary;   //  boundary condition vect
  getInput(P, boundary);       //  get boundary conditions
  for(int i=0; i<boundary.size(); i++) {
    sBC & b = boundary[i];
    if(b.isAnalytic)
      P.setAnalyticBoundaryValues(b.name, b.type, b.equation);
    else
      P.setFieldBoundaryValues(b.name, b.type, b.field);
  }
  fDebug << "\n 5. set boundary cond. : t= " << T.elapsed() << std::endl;
  
  
  // **********************************************************
  for(int itime=0; itime<20; itime++) { // time loop

    //   6. solving
    T.start();
    P.solve();
    fDebug << "\n 6. solve : t= " << T.elapsed() << std::endl;
    // **********************************************************
    // 7. solution storage
    T.start();
    sol = P.getOutputField("IVec u +JVec v +KVec w ");
//     solv = P.getOutputField("v");
    // solp = P.getOutputField("p");
    fDebug << "\n 7. get output : t= " << T.elapsed() << std::endl;
    // **********************************************************
    // 8. solution write
    T.start();  
   std::string s = "RESU/ns_";  s += prefix;  s += ".med";
   std::cout << " File "<< s <<"\n";
   std::ostringstream  s2; 
   s2<<"RESU/ns_" << prefix<<"."<<itime<<".med";
   MEDLoader::WriteField(s2.str().c_str(), sol, false);
//    MEDLoader::WriteField(s2.str().c_str(), solv, false);
   // MEDLoader::WriteField(s.c_str(), solp, true);
   
  } // end time
//   // boundary write
//   for(int i=0; i<boundary.size(); i++) {
//     bdy = P.getValuesOnBoundary(boundary[i].name, "u");
//     //bdyv = P.getValuesOnBoundary(boundary[i].name, "v");
//     std::ostringstream s;  s << "RESU/testNavier_" << prefix << "_" << i << ".med";
//     MEDLoader::WriteField(s.str().c_str(), bdy, false);
//     //MEDLoader::WriteField(s.str().c_str(), bdyv, false);
//   }

  P.terminate();
  fDebug << "Elapsed (terminate)          : " << T.elapsed() << std::endl;

  if(sol) {
//     unsigned int n = sol->getNumberOfTuples()/2-1;
//     result = fabs(sol->getIJ(n,0) - good) < 1e-5;
//     std::ostringstream s;
//     s << "convdiff: result at node " << n
//       << " = " << sol->getIJ(n,0) << " expected " << good;
//     s << std::endl;
//     std::cerr << s.str() << std::endl;
    std::vector<const MEDCouplingFieldDouble *> vsol(2);
    vsol[0] = sol;
    vsol[1] = solv;
    // vsol[2] = solv;
    MEDCouplingFieldDouble::WriteVTK("RESU/testNavier_prova.vtk", vsol);
    sol->decrRef();
  }


  return 0;
}

// get user boundary conditions
void getInput(LIBMESH &P, std::vector<sBC> & boundary) {

  // boundary
  std::vector<std::string> bcNames = P.getBoundaryNames();
  int nBC = bcNames.size();
  boundary.resize(nBC);

  for(int i=0; i<nBC; i++) {
    boundary[i].name = bcNames[i];
    boundary[i].type = "Neumann";
    boundary[i].isAnalytic = true;

    if(bcNames[i]=="Group_1") {
      boundary[i].type = "Dirichlet";
      boundary[i].isAnalytic = true;
      boundary[i].field = NULL;
      boundary[i].equation = "IVec*0.0 + JVec*0.0";
    }
    if(bcNames[i]=="Group_2") {
      boundary[i].type = "Dirichlet";
      boundary[i].isAnalytic = true;
      boundary[i].field = NULL;
      boundary[i].equation = "IVec*0.0 + JVec*10.0";
    }
    if(bcNames[i]=="Group_3") {
      boundary[i].type = "Dirichlet";
      boundary[i].isAnalytic = true;
      boundary[i].field = NULL;
      boundary[i].equation = "IVec*0.0 + JVec*0.0";
    }
  }

  return;
}

