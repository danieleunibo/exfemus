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
  MEDCouplingFieldDouble * field1;
};

void getInput(LIBMESH &, std::vector<sBC> & boundary);



// ==================================================
int main(
  int argc,
  char **argv
) { // ==============================================

  Timer T;         // time counter
  bool result;     // logic variable
  fDebug.open("navierstokes3D", MPI_COMM_NULL); // file log -> test2.log

  // Reading data file ----------------------------------------------
  std::string dataFile = "MESH/cube_hex20.med";// Initial dataFile
  if(argc > 1)  dataFile = argv[1];            // dataFile from argv[1]

  // only filename ("square2") ------------------------------
  std::string prefix;
  size_t p = dataFile.find_last_of('/');
  size_t q = dataFile.find_last_of('.');
  p = (p == std::string::npos) ? 0 : p+1;
  prefix = dataFile.substr(p, q-p);
  // --------------------------------------------------------
  // coupling field sol  -------------------------------
   std::vector<const MEDCouplingFieldDouble *> vsol(3);

   // **********************************************************
  // Preprocessing
  // 1- initialize -> LIBMESH P
  T.start();
  LIBMESH P;
  fDebug << "\n 1- initialize: t= " << T.elapsed() << std::endl;
  // **********************************************************
  // 2. set problem  -> P.setType
  T.start();
  P.setType("navierstokes3D");
  fDebug << "\n 2- set problem: t= " << T.elapsed() << std::endl;
  // **********************************************************
  // 3. set mesh  -> P.setType
  T.start();
  P.setMesh(dataFile);
  fDebug << "\n 3- set mesh: t= " << T.elapsed() << std::endl;
  // **********************************************************
  // 4. set source  ->  P.setAnalyticSource("0")
  T.start();
//   P.setAnalyticSource("IVec*0.0 + JVec*0.0 +KVec *0.0");
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
    vsol[0] = P.getOutputField("u");
    vsol[1] = P.getOutputField("v");
    vsol[2] = P.getOutputField("w");
    // solp = P.getOutputField("p");
    fDebug << "\n 7. get output : t= " << T.elapsed() << std::endl;
    
    // **********************************************************
    // 8. solution write
    T.start();  
   std::string s = "RESU/ns_";  s += prefix;  s += ".med"; std::cout << " File "<< s <<"\n";
   std::ostringstream  s2; 
   s2<<"RESU/ns_" << prefix<<"."<<itime;
   MEDCouplingFieldDouble::WriteVTK((s2.str()+".vtk").c_str(), vsol);
  } // end time
//   // boundary write
//   for(int i=0; i<boundary.size(); i++) {
//     bdy = P.getValuesOnBoundary(boundary[i].name, "u");
//     //bdyv = P.getValuesOnBoundary(boundary[i].name, "v");
//     std::ostringstream s;  s << "RESU/testNavier_" << prefix << "_" << i << ".med";
//     MEDLoader::WriteField(s.str().c_str(), bdy, false);
//     //MEDLoader::WriteField(s.str().c_str(), bdyv, false);
//   }
  // Final ending -----------------------------------------------------
  // End MPI and problem solution 
  P.terminate();
  fDebug << "Elapsed (terminate)          : " << T.elapsed() << std::endl;
   // Print final solution in MED format  -----------------------------------
  if(vsol[0]) {
//     MEDCouplingFieldDouble::WriteVTK("RESU/testNavier_prova.vtk", vsol);
    MEDLoader::WriteField("RESU/navierstokes.med", vsol[0], true);
    MEDLoader::WriteField("RESU/navierstokes.med", vsol[1], false);
    MEDLoader::WriteField("RESU/navierstokes.med", vsol[2], false);
   }
 T.end();  
  vsol.clear();
  return 0;
}



//========================================================================
// This function gets the user boundary conditions
void getInput(
  LIBMESH &P,                      // Problem pointer     <-
  std::vector<sBC> & boundary      // boundary conds      ->
) {

  // boundary
  std::vector<std::string> bcNames = P.getBoundaryNames();
  int nBC = bcNames.size();
  boundary.resize(nBC);

  for(int i=0; i<nBC; i++) {
    boundary[i].name = bcNames[i];
    boundary[i].type = "Neumann";
    boundary[i].isAnalytic = true;

    if(bcNames[i]=="lateral") {
      boundary[i].type = "Dirichlet";
      boundary[i].isAnalytic = true;
      boundary[i].field = NULL;
      boundary[i].equation = "IVec*0.0 + JVec*0.0+ KVec*0.0";
    }
    if(bcNames[i]=="bottom") {
      boundary[i].type = "Dirichlet";
      boundary[i].isAnalytic = true;
      boundary[i].field = NULL;
      boundary[i].equation = "IVec*0.0 + JVec*0.0+ KVec*1.0";
    }
//     if(bcNames[i]=="Group_3") {
//       boundary[i].type = "Dirichlet";
//       boundary[i].isAnalytic = true;
//       boundary[i].field = NULL;
//       boundary[i].equation = "IVec*0.0 + JVec*0.0+ KVec*0.0";
//     }
  }

  return;
}


