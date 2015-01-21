#include <stdexcept>
#include "LIBMESH.h"
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
#include "Debug.h"
#include "Timer.h"
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
void getInput1(LIBMESH &,LIBMESH &, std::vector<sBC> & boundary);

int main(int argc, char **argv) {

  Timer T;         // time counter
  bool result;     // logic variable
  fDebug.open("convdiff", MPI_COMM_NULL); // file log -> test2.log

  // Reading data file ----------------------------------------------
  std::string dataFile =  "../MESH/Cube_down.med"; // Initial dataFile
   std::string dataFile1 =  "../MESH/Cube_up.med"; // Initial dataFile       
//   double good = -0.0530494;                    // Initial result check data
  if(argc > 1)  dataFile = argv[1];            // dataFile from argv[1]
//   if(argc > 2)  good = std::strtod(argv[2], NULL); // Check data  from argv[2]

  // only filename ("square2") ------------------------------
  std::string prefix;
  size_t p = dataFile.find_last_of('/');
  size_t q = dataFile.find_last_of('.');
  p = (p == std::string::npos) ? 0 : p+1;
  prefix = dataFile.substr(p, q-p);
   std::string prefix1;
  p = dataFile1.find_last_of('/');
  q = dataFile1.find_last_of('.');
  p = (p == std::string::npos) ? 0 : p+1;
  prefix1 = dataFile1.substr(p, q-p);
  // --------------------------------------------------------
  // coupling fields sol bdy  ------------------------------
  MEDCouplingFieldDouble *sol = NULL, *bdy = NULL;  
  MEDCouplingFieldDouble *sol1 = NULL, * bdy1 = NULL;

  
  
  // Preprocessing  ------------------------------
  // 1- initialize -> LIBMESH P
  T.start(); 
  LIBMESH P; 
  LIBMESH P1;  
  fDebug << "\n 1- initialize: t= " << T.elapsed() << std::endl;
             
  // 2. set problem  -> P.setType
   T.start();
   P.setType("convdiff");   
   P1.setType("convdiff"); 
   fDebug << "\n 2- set problem: t= " << T.elapsed() << std::endl;
//   // 3. set mesh  -> P.setType
   T.start();  
   P.setMesh(dataFile); 
               P1.setMesh(dataFile1); 
              fDebug << "\n 3- set mesh: t= " << T.elapsed() << std::endl;
//   // 4. set source  ->  P.setAnalyticSource("0")
  T.start();  P.setAnalyticSource("0.");
             P1.setAnalyticSource("0.");  
              fDebug << "\n 4- set source: t= " << T.elapsed() << std::endl;
//   // 5. boundary cds.  ->  P.setAnalyticBoundaryValues or P.setFieldBoundaryValues
// //  
 T.start();  std::vector<sBC> boundary;  
             std::vector<sBC> boundary1; 
  
  fDebug << "\n 5. set boundary cond. : t= " << T.elapsed() << std::endl;
 getInput(P, boundary);
 
  // time loop
  for(int itime=0;itime<60;itime++){
    T.start();
    // boundary
   
     for(int i=0; i<boundary.size(); i++) {
    sBC & b = boundary[i];
    if(b.isAnalytic)
      P.setAnalyticBoundaryValues(b.name, b.type, b.equation);
    else
      P.setFieldBoundaryValues(b.name, b.type, b.field);
  }
  // Solving --------------------------------------------------
  //   6. solving
    P.solve();  
     getInput1(P1, P,  boundary1);
    // boundary
    for(int i=0; i<boundary1.size(); i++) {
    sBC & b1 = boundary1[i];
    if(b1.isAnalytic)
       P1.setAnalyticBoundaryValues(b1.name, b1.type, b1.equation);
   else
     P1.setFieldBoundaryValues(b1.name, b1.type, b1.field);
  }
  // Solving --------------------------------------------------
  //   6. solving
    P1.solve();  
    
    fDebug << "\n 6. solve : t= " << T.elapsed() << std::endl;
  // Postprocessing --------------------------------------------------
  // 7. solution storage -----------------------------------------
   T.start();  
   sol = P.getOutputField("u");
   sol1 = P1.getOutputField("u");
   fDebug << "\n 7. get output : t= " << T.elapsed() << std::endl;
//   // 8. solution write
   T.start(); 
   std::string s = "RESU/convdiff_sol_bob_";  s += prefix;  s += ".med";
   std::string s1 = "RESU/convdiff_sol_bob_";  s1 += prefix1;  s1 += ".med";
   std::cout << " File "<< s <<"\n"; std::cout << " File "<< s1 <<"\n";
   std::ostringstream  s2; s2<<"RESU/convdiff_sol_bob" << prefix<<"."<<itime<<".med";
   std::ostringstream  s3; s3<<"RESU/convdiff_sol_bob" << prefix1<<"."<<itime<<".med";
    
   MEDLoader::WriteField(s2.str().c_str(), sol, true);
   MEDLoader::WriteField(s3.str().c_str(), sol1, true);
   
      std::vector<const MEDCouplingFieldDouble *> vsol(1);
      vsol[0] = sol;
      MEDCouplingFieldDouble::WriteVTK((s2.str()+".vtk").c_str(), vsol);
      vsol[0] = sol1;
      MEDCouplingFieldDouble::WriteVTK((s3.str()+".vtk").c_str(), vsol); 
  } // end time
   // boundary write
//    for(int i=0; i<boundary.size(); i++) {
//      bdy = P.getValuesOnBoundary(boundary[i].name, "u");
//      std::ostringstream s;  s << "RESU/testconvdiff_bdy_bob_" << prefix << "_" << i << ".med";
//      MEDLoader::WriteField(s.str().c_str(), bdy, true);
//    }

  P.terminate(); 
  P1.terminate();
  fDebug << "Elapsed (terminate)          : " << T.elapsed() << std::endl;

  if(sol) {
//     unsigned int n = sol->getNumberOfTuples()/2-1;
//     result = fabs(sol->getIJ(n,0) - good) < 1e-5;
//     std::ostringstream s;
//     s << "convdiff sol: result at node " << n
//       << " = " << sol->getIJ(n,0) << " expected " << good;
//     s << std::endl;
//     std::cerr << s.str() << std::endl;
    std::vector<const MEDCouplingFieldDouble *> vsol(1);
    vsol[0] = sol;
    MEDCouplingFieldDouble::WriteVTK("RESU/convdiff1.vtk", vsol);
    sol->decrRef();
  }
   if(sol1) {
//     unsigned int n = sol1->getNumberOfTuples()/2-1;
//     result = fabs(sol1->getIJ(n,0) - good) < 1e-5;
//     std::ostringstream s;
//     s << "convdiff sol1: result at node " << n
//       << " = " << sol1->getIJ(n,0) << " expected " << good;
//     s << std::endl;
//     std::cerr << s.str() << std::endl;
    std::vector<const MEDCouplingFieldDouble *> vsol(1);
    vsol[0] = sol1;
    MEDCouplingFieldDouble::WriteVTK("RESU/convdiff2.vtk", vsol);
    sol->decrRef();
  }
  

  return 0;
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
  
  //  2d================================
  
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
//   if( bcNames[i]=="Down_bottom"){
//   boundary[i].type = "Dirichlet";
//   boundary[i].isAnalytic = true;
//   boundary[i].field = NULL;
//   boundary[i].equation = "2";
//   }
//     
//   if( bcNames[i]=="Down_right"){
//   boundary[i].type = "Dirichlet";
//   boundary[i].isAnalytic = true;
// //   boundary[i].field = infield ;
// //   boundary[i].field = NULL;
//   boundary[i].equation = "0";
//  }
//   if( bcNames[i]=="Down_left"){
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
  
  
  
  /// 3D =========================================
  
  // Field
//  std::string inFieldFile = "MESH/in_square2_2.med";
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
    
    
  if( bcNames[i]=="Bottom_down"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
  boundary[i].field = NULL;
  boundary[i].equation = "1.";
  }
    
  if( bcNames[i]=="right_down"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
//   boundary[i].field = infield ;
//   boundary[i].field = NULL;
  boundary[i].equation = "0.";
 }
  if( bcNames[i]=="left_down"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
//   boundary[i].field = infield ;
//   boundary[i].field = NULL;
  boundary[i].equation = "0.";
 }
 if( bcNames[i]=="Up_down"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
//   boundary[i].field = infield ;
//   boundary[i].field = NULL;
  boundary[i].equation = "0.";
 }
  if( bcNames[i]=="Down_down"){
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
void getInput1(LIBMESH &P,LIBMESH &P_OLD, std::vector<sBC> & boundary)
{
  
  MEDCouplingFieldDouble *bdy=NULL;  
  bdy = P_OLD.getValuesOnBoundary("share", "u");
  


//   // boundary condition construction   
//   std::vector<std::string> bcNames = P.getBoundaryNames();
//   int nBC = bcNames.size();
//   boundary.resize(nBC);
  
//   for (int i=0; i<nBC; i++) {
//     boundary[i].name = bcNames[i];
//     boundary[i].type = "Neumann";
//     boundary[i].isAnalytic = true;
//     
//     
//   if( bcNames[i]=="share"){
//   boundary[i].type = "Dirichlet";
//   boundary[i].isAnalytic = false;
// //   boundary[i].field = NULL;
// //   boundary[i].equation = "2";
//   boundary[i].field = bdy ;
//   }
//     
//   if( bcNames[i]=="Up_left"){
//   boundary[i].type = "Dirichlet";
//   boundary[i].isAnalytic = true;
//   //   boundary[i].field = bdy ;
//   boundary[i].field = NULL;
//   boundary[i].equation = "0.0";
//  }
//    if( bcNames[i]=="Up_right"){
//   boundary[i].type = "Dirichlet";
//   boundary[i].isAnalytic = true;
//   //   boundary[i].field = bdy ;
//   boundary[i].field = NULL;
//   boundary[i].equation = "0.0";
//  }
// //   if( bcNames[i]=="share"){
// // //   boundary[i].type = "Dirichlet";
// // //   boundary[i].isAnalytic = true;
// // //   boundary[i].equation = "2.0";
// // //   boundary[i].field = NULL;
// //     } 
//   }
  
  
  ////// 3D ////////////////
 
 
  // boundary condition construction   
  std::vector<std::string> bcNames = P.getBoundaryNames();
  int nBC = bcNames.size();
  boundary.resize(nBC);
  
  for (int i=0; i<nBC; i++) {
    boundary[i].name = bcNames[i];
    boundary[i].type = "Neumann";
    boundary[i].isAnalytic = true;
    
    
  if( bcNames[i]=="Up_up"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
  boundary[i].field = NULL;
  boundary[i].equation = "0.0";
  }
    
  if( bcNames[i]=="right_up"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
  boundary[i].field = NULL;
  boundary[i].equation = "0.0";
 }
  if( bcNames[i]=="left_up"){
  boundary[i].type = "Dirichlet";
  boundary[i].isAnalytic = true;
  boundary[i].field = NULL;
  boundary[i].equation = "0.0";
  }
    
  if( bcNames[i]=="Down_up"){
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

/*
void getInput(LIBMESH &P, std::vector<sBC> & boundary) {

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
//     boundary[i].type = "Neumann";
     boundary[i].type = "Dirichlet";
    boundary[i].isAnalytic = true;
     boundary[i].equation = "0.";
    
  if( bcNames[i]=="top"){
   boundary[i].type = "Dirichlet";
   boundary[i].isAnalytic = true;
   boundary[i].field = NULL;
   boundary[i].equation = "1.";
   }
//     
//   if( bcNames[i]=="lat1"){
//   boundary[i].type = "Dirichlet";
//   boundary[i].isAnalytic = true;
// //   boundary[i].field = infield ;
// //   boundary[i].field = NULL;
//   boundary[i].equation = "0.";
//  }
//   if( bcNames[i]=="share"){
//  //   boundary[i].type = "Dirichlet";
// //   boundary[i].isAnalytic = true;
// //   boundary[i].equation = "2.0";
// //   boundary[i].field = NULL;
//     } 
  }
  
  
 return; 
}
*/
