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
#include "BCTypes.hxx"
#define FDEBUG 0

using namespace ParaMEDMEM;

void getInput(LIBMESH &, std::vector<BCdata> & boundary);

int main(int argc, char **argv) {

  Timer T;         // time counter
  bool result;     // logic variable
  fDebug.open("convdiff", MPI_COMM_NULL); // file log -> test2.log

  // Reading data file ----------------------------------------------
  std::string dataFile =  "MESH/cubeup.med"; // Initial dataFile
  if (argc > 1)  dataFile = argv[1];           // dataFile from argv[1]

  // only filename ("square2") ------------------------------
  std::string prefix;
  size_t p = dataFile.find_last_of('/');
  size_t q = dataFile.find_last_of('.');
  p = (p == std::string::npos) ? 0 : p+1;
  prefix = dataFile.substr(p, q-p);
  // --------------------------------------------------------
  // coupling fields sol bdy  ------------------------------
  MEDCouplingFieldDouble *sol = NULL, *bdy = NULL;

  // Preprocessing  ------------------------------
  // 1- initialize -> LIBMESH P
  T.start(); 
  LIBMESH P;  
  fDebug << "\n 1- initialize: t= " << T.elapsed() << std::endl;
  // 2. set problem  -> P.setType
  T.start();
  P.setType("convdiff");
  fDebug << "\n 2- set problem: t= " << T.elapsed() << std::endl;
  // 3. set mesh  -> P.setType
  T.start();  
  P.setMesh(dataFile); 
  fDebug << "\n 3- set mesh: t= " << T.elapsed() << std::endl;
  // 4. set source  ->  P.setAnalyticSource("0")
  T.start();  
  P.setAnalyticSource("1.");  
  fDebug << "\n 4- set source: t= " << T.elapsed() << std::endl;
  // 5. boundary cds.  ->  P.setAnalyticBoundaryValues or P.setFieldBoundaryValues
  T.start();  
  std::vector<BCdata> boundary;   getInput(P, boundary);
  for (int i=0; i<boundary.size(); i++) {
    BCdata & b = boundary[i];
    if (b.isAnalytic)
      P.setAnalyticBoundaryValues(b.name, b.type, b.equation);
    else
      P.setFieldBoundaryValues(b.name, b.type, b.field);
  }
  fDebug << "\n 5. set boundary cond. : t= " << T.elapsed() << std::endl;

  // time loop
  for (int itime=0; itime<10; itime++) {

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
    std::string s = "RESU/convdiff_sol";  s += prefix;  s += ".med";
    std::cout << " File "<< s <<"\n";
    std::ostringstream  s2;
    s2<<"RESU/convdiff_sol_" << prefix<<"."<<itime;
    MEDLoader::WriteField((s2.str()+".med").c_str(), sol, true);

    if (sol) {
      std::vector<const MEDCouplingFieldDouble *> vsol(1);
      vsol[0] = sol;
      MEDCouplingFieldDouble::WriteVTK((s2.str()+".vtk").c_str(), vsol);
      sol->decrRef();
  }

  } // end time
  // boundary write
//    for(int i=0; i<boundary.size(); i++) {
//      bdy = P.getValuesOnBoundary(boundary[i].name, "u");
//      std::ostringstream s;  s << "RESU/testconvdiff_bdy_bob_" << prefix << "_" << i << ".med";
//      MEDLoader::WriteField(s.str().c_str(), bdy, true);
//    }

  P.terminate();
  fDebug << "Elapsed (terminate)          : " << T.elapsed() << std::endl;


  return 0;
}


void getInput(LIBMESH &P, std::vector<BCdata> & boundary) {

//   // Field
// //   std::string inFieldFile = "MESH/in_square2_2.med";
// //    std::string inFieldFile = "RESU/test2_bdy_square2_2_pass.med";
//
// //   std::vector<std::string> InFieldNames = MEDLoader::GetAllFieldNames(inFieldFile.c_str());
// //   std::vector<std::string> InMeshNames  = MEDLoader::GetMeshNames(inFieldFile.c_str());
//
// //   MEDCouplingFieldDouble *infield =
// //      MEDLoader::ReadField(ON_NODES, inFieldFile.c_str(), InMeshNames[0].c_str(),
// //                           0, InFieldNames[0].c_str(), -1, -1);
//
//   // boundary condition construction
//   std::vector<std::string> bcNames = P.getBoundaryNames();
//   int nBC = bcNames.size();
//   boundary.resize(nBC);
//
//   for (int i=0; i<nBC; i++) {
//     boundary[i].name = bcNames[i];
//     boundary[i].type = "Neumann";
//     boundary[i].isAnalytic = true;
//
//
//   if( bcNames[i]=="inlet"){
// //   boundary[i].type = "Dirichlet";
// //   boundary[i].isAnalytic = true;
// //   boundary[i].field = NULL;
// //   boundary[i].equation = "1";
//   }
//
//   if( bcNames[i]=="lateral"){
//   boundary[i].type = "Dirichlet";
//   boundary[i].isAnalytic = true;
// //   boundary[i].field = infield ;
// //   boundary[i].field = NULL;
//   boundary[i].equation = "0";
//  }
//   if( bcNames[i]=="share"){
//  //   boundary[i].type = "Dirichlet";
// //   boundary[i].isAnalytic = true;
// //   boundary[i].equation = "2.0";
// //   boundary[i].field = NULL;
//     }
//   }
//
//

  /// 3D  ////////////

  // Field
//   std::string inFieldFile = "MESH/in_square2_2.med";
//    std::string inFieldFile = "RESU/test2_bdy_square2_2_pass.med";

//   std::vector<std::string> InFieldNames = MEDLoader::GetAllFieldNames(inFieldFile.c_str());
//   std::vector<std::string> InMeshNames  = MEDLoader::GetMeshNames(inFieldFile.c_str());

//   MEDCouplingFieldDouble *infield =
//      MEDLoader::ReadField(ON_NODES, inFieldFile.c_str(), InMeshNames[0].c_str(),
//                           0, InFieldNames[0].c_str(), -1, -1);

  // boundary condition construction
  std::vector<std::string> bcNames = P.getBoundaryNames();
  int nBC = bcNames.size();
  boundary.resize(nBC);

  for (int i=0; i<nBC; i++) {
    boundary[i].name = bcNames[i];
    boundary[i].type = "Neumann";
    boundary[i].isAnalytic = true;
    boundary[i].equation = "0.";

    if (bcNames[i]=="top") {
      boundary[i].type = "Dirichlet";
      boundary[i].isAnalytic = true;
      boundary[i].field = NULL;
      boundary[i].equation = "2.";
    }
//
   if( bcNames[i]=="lat1" || bcNames[i]=="lat2"){
   boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
//   boundary[i].field = infield ;
//   boundary[i].field = NULL;
  boundary[i].equation = "0.";
  }
//   if( bcNames[i]=="share"){
//    boundary[i].type = "Dirichlet";
//   boundary[i].isAnalytic = true;
//   boundary[i].equation = "2.0";
//   boundary[i].field = NULL;
//     }
  }


  return;
}

