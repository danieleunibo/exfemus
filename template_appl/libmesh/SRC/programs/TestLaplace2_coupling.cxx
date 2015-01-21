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
void getInput2(LIBMESH &,LIBMESH &, std::vector<sBC> & boundary);

int main(int argc, char **argv){

  Timer T;
  bool result;

  fDebug.open("test2", MPI_COMM_NULL);

//   std::string dataFile = "MESH/square2.med";
//    std::string dataFile = "MESH/rhs.med";
   std::string dataFile = "MESH/cubeup.med";
   
   
  double good = -0.0530494;

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
  P.setAnalyticSource("40*x*x*x*x");
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
    
  T.start();
  P.solve();
  fDebug << "Elapsed (solve)              : " << T.elapsed() << std::endl;
    
  T.start();
  sol = P.getOutputField("u");
  fDebug << "Elapsed (get output)         : " << T.elapsed() << std::endl;
  
  T.start();
  std::string s = "RESU/test2_sol_";
  s += prefix;
  s += ".med";
//   MEDLoader::WriteField(s.c_str(), sol, true);

  for (int i=0; i<boundary.size(); i++) {
    bdy = P.getValuesOnBoundary(boundary[i].name, "u");
    std::ostringstream s;
    s << "RESU/test2_bdy_" << prefix << "_" << i << ".med";
//     MEDLoader::WriteField(s.str().c_str(), bdy, true);
  }


  
//   //////////////////////////////////////////////////
// // // // //   MESH 2
  ////////////////////////////////////////////
//   ////////////////////////////////////////////
  
    MEDCouplingFieldDouble * sol1 = NULL, * bdy1 = NULL;
//  dataFile = "MESH/lhs.med";
 dataFile = "MESH/cubedown.med";
 
  T.start();
  LIBMESH P1;
  fDebug << "Elapsed (initialize)         : " << T.elapsed() << std::endl;
    
  T.start();
  P1.setType("Poisson");
  fDebug << "Elapsed (problem type)       : " << T.elapsed() << std::endl;
    
  T.start();
  P1.setMesh(dataFile);

  fDebug << "Elapsed (set mesh)           : " << T.elapsed() << std::endl;
  
  T.start();
  P1.setAnalyticSource("40*x*x*x*x");
  fDebug << "Elapsed (set source)         : " << T.elapsed() << std::endl;

    
  T.start();
  // get boundary condition <- getInput 
  std::vector<sBC> boundary1;  getInput2(P1, P,boundary1);
  
  for (int i=0; i<boundary1.size(); i++) {
    sBC & b = boundary1[i];
    // bc are 1) analitic 2) numeric ( MEDCouplingFieldDouble::field)
    if (b.isAnalytic) P1.setAnalyticBoundaryValues(b.name, b.type, b.equation);
    else                 P1.setFieldBoundaryValues(b.name, b.type, b.field);
  }

  fDebug << "Elapsed (set boundary cond.) : " << T.elapsed() << std::endl;
    
  T.start();
  P1.solve();
  fDebug << "Elapsed (solve)              : " << T.elapsed() << std::endl;
    
  T.start();
  sol1 = P1.getOutputField("u");
  fDebug << "Elapsed (get output)         : " << T.elapsed() << std::endl;
  
  T.start();
  std::string s1 = "RESU/test2_sol1_";
  s1 += prefix;
  s1 += ".med";
//   MEDLoader::WriteField(s1.c_str(), sol1, true);

  for (int i=0; i<boundary1.size(); i++) {
    bdy1 = P1.getValuesOnBoundary(boundary1[i].name, "u");
    std::ostringstream s1;
    s1 << "RESU/test2_bdy1_" << prefix << "_" << i << ".med";
//     MEDLoader::WriteField(s1.str().c_str(), bdy, true);
  }

   P.terminate();
   P1.terminate(); 
  
  
  //   //////////////////////////////////////////////////
// // // // //  END  MESH 2
  ////////////////////////////////////////////
//   ////////////////////////////////////////////
  
  
  
  fDebug << "Elapsed (terminate)          : " << T.elapsed() << std::endl;
  
  if (sol) {
    unsigned int n = sol->getNumberOfTuples()/2-1;
    result = fabs(sol->getIJ(n,0) - good) < 1e-5;
      std::ostringstream s;
      s << "TestLaplace2: result at node " << n 
	<< " = " << sol->getIJ(n,0) << " expected " << good;
      s << std::endl;
      std::cerr << s.str() << std::endl;

    sol->decrRef();
  }

  return result ? 0 : 1;
}

// =========================================================
// get input boundary condition
void getInput(LIBMESH &P, std::vector<sBC> & boundary)
{
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
    
    
  if( bcNames[i]=="top"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
  boundary[i].field = NULL;
  boundary[i].equation = "1.";
  }
    
  if( bcNames[i]=="lat1"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
//   boundary[i].field = infield ;
//   boundary[i].field = NULL;
  boundary[i].equation = "0.";
 }
  if( bcNames[i]=="share"){
 //   boundary[i].type = "Dirichlet";
//   boundary[i].isAnalytic = true;
//   boundary[i].equation = "2.0";
//   boundary[i].field = NULL;
    } 
  }
  

  return;
 }
 // =========================================================
// get input boundary condition
void getInput2(LIBMESH &P,LIBMESH &P_OLD, std::vector<sBC> & boundary)
{
  
//   MEDCouplingFieldDouble *bdy=NULL;
//   // Field
// //   std::string inFieldFile = "MESH/in_square2_2.med";
// //    std::string inFieldFile = "RESU/test2_bdy_square2_2_pass.med";
// bdy = P_OLD.getValuesOnBoundary("share", "u");
// // const char func[]="u*20+2";
// // bdy->applyFunc(func);
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
//   if( bcNames[i]=="outlet"){
//   boundary[i].type = "Dirichlet";
//   boundary[i].isAnalytic = false;
// //   boundary[i].field = NULL;
// //   boundary[i].equation = "2";
//   boundary[i].field = bdy ;
//   }
    
//   if( bcNames[i]=="lateral"){
//   boundary[i].type = "Dirichlet";
//   boundary[i].isAnalytic = true;
//   //   boundary[i].field = bdy ;
//   boundary[i].field = NULL;
//   boundary[i].equation = "0.0";
//  }
//   if( bcNames[i]=="share"){
// //   boundary[i].type = "Dirichlet";
// //   boundary[i].isAnalytic = true;
// //   boundary[i].equation = "2.0";
// //   boundary[i].field = NULL;
//     } 
//   }
  
  
  ////// 3D ////////////////
  
   
  MEDCouplingFieldDouble *bdy=NULL;
  bdy = P_OLD.getValuesOnBoundary("share", "u");
// const char func[]="u*20+2";
// bdy->applyFunc(func);
  
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
    
    
  if( bcNames[i]=="bottom"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
  boundary[i].field = NULL;
  boundary[i].equation = "1.";
  }
    
  if( bcNames[i]=="lat1"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
  boundary[i].field = NULL;
  boundary[i].equation = "0.0";
 }
  if( bcNames[i]=="share"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = false;
//   boundary[i].field = NULL;
//   boundary[i].equation = "2.0";
  boundary[i].field = bdy;
    } 
  }     
        
    
    
    
  
  
  

  return;
 }
