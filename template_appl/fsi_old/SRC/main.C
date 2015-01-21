// libc+++ include
#include <iostream>
#include <cstdlib>
#include <sstream>

// configuration files -------------------------
#include   "Printinfo_conf.h"

// LibMesh library included ------------------------------
// #ifdef LM_INIT
// #include "libmesh.h" // for Libmesh library
// #endif

// solver library -------------------------------------
#include  "Solverlib_conf.h"  // Solver library options 
// Petsc
#ifdef HAVE_PETSCM
#include "petsc.h" // for Petsc solver
#endif
// Mpi
#ifdef HAVE_MPI
#include <mpi.h>   //For MPI_COMM_WORLD
#endif

// include local class
#include "MGFemusInit.h"
#include "MGUtils.h"
#include "MGSystem.h"
#include "MGGeomEl.h"
#include "MGMesh.h"
#include "MGFEMap.h"
#include "MGFE.h"
#include "MGEquationsSystem.h"
#include "MGTimeLoop.h"
#include "FEMUS.h"

#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"

// ==========================================================================
  ParaMEDMEM::MEDCouplingFieldDouble *build_field(
   const  std::string & path_mesh,
   const  std::string & expression,
    int id_interface
  );
#endif
/// Set up
// =======================================
// Main program
// =======================================

int main(int argc, char** argv) {

  argc = argc ; argv=argv;  // no unused warning
 
  std::cout<<" ============ MGUtils ===================================== \n";
  std::cout<<" =========================================================== \n";
  // setting MGUtils -> parameters and file name ------------------------
  std::vector<MGUtils*> mgutils;
  std::string mesh_nameP[NUM_MESH];
  std::ostringstream filenameP[2];  std::ostringstream osfilenameP[2];

  for(int i_mesh=0; i_mesh< NUM_MESH; i_mesh++) {
    // MGUtils constructor ----------------------------------------------------
    mgutils.push_back(new MGUtils(i_mesh+1));
    // mesh names -------------------------------------------------------------
    mesh_nameP[i_mesh]= mgutils[i_mesh]->get_file("F_MESH_READ"); // name mesh
    int posP = mesh_nameP[i_mesh].find(".");  // position of "live" in str
    filenameP[i_mesh] <<   mesh_nameP[i_mesh].substr(0,posP)  << "_gen.med" ;
    osfilenameP[i_mesh]<< mgutils[i_mesh]->_mesh_dir <<filenameP[i_mesh].str();
    std::cout<<" \n P mesh file "<< i_mesh+1 << "= "
                                  << osfilenameP[i_mesh].str().c_str() <<"\n "; 
  }
  std::cout<<" ============ end loop mesh ================================ \n";
  std::cout<<" =========================================================== \n";
  
// FEM class -----------------------------------------------------------
  MGFEMap *mgfemap; mgfemap=new MGFEMap();
  // MGFEMap mg_femap;
  MGFE *dfe_q;    dfe_q=new MGFE(2,ELTYPE); dfe_q->init_qua();
  mgfemap->set_FE(dfe_q); //// initialize quadratic fem
  MGFE *dfe_l;  dfe_l=new MGFE(1,ELTYPE); dfe_l->init_lin();
  mgfemap->set_FE(dfe_l); //initialize linear fem
  MGFE *dfe_k; dfe_k=new MGFE(0,ELTYPE);  dfe_k->init_pie();
  mgfemap->set_FE(dfe_k); //  initialize piecewise fem

  // MGGeomEl ----------------------------------------------------------
  MGGeomEl *mggeomel;  mggeomel=new  MGGeomEl();
  
  // system problem =========================================================
  std::vector<FIELDS> myproblemP; myproblemP.resize(1);
  // Problem to solve for each mesh
  myproblemP[0]=FS_F; 
//   myproblemP[1]=DS_F;
  
  // system 1
  // MGFemusInit --------------------------------------------------------------
  FEMUS P;                                        // constructor
  P.init_param(*mgutils[0]);                      // init parameter
  P.init_fem(*mggeomel,*mgfemap);                 // init fem      
  // setting mesh -------------------------------------------------------------
//   P.setMedMesh(osfilenameP[0].str().c_str());     // set med-mesh
  P.setMesh();                                    // set MGmesh   
  // setting system -----------------------------------------------------------
  P.setSystem(myproblemP);                         // set system
  P.init_interface(20,20, 2,filenameP[0].str().c_str());
//   P.init_interface(51,50, 2,filenameP[0].str().c_str());
//   P.init_interface(211,211,2,filenameP[0].str().c_str());
//   P.init_interface(210,211,2,filenameP[0].str().c_str());
//   P.init_interface(1,1,  2,filenameP[0].str().c_str(),false);
//   P.init_interface(2,2,  2,filenameP[0].str().c_str(),false);
//  P.init_interfaces(filenameP[0].str().c_str(),0);  // init interfaces
//   P.set_param_BCzone();                           // define param bc conditions 
  P.setAnalyticBoundaryValues(20,2, "IVec*(0.05)+JVec*0."); // write in the interfacefunction.field
  P.write_Boundary_value(20,"FSI0",2);             // write inside x_old the boundary values
  
//   P.setAnalyticBoundaryValues(51,1, "30");
//   P.write_Boundary_value(51,"T",1);               // write inside x_old the boundary values
// 
// //   P.setFieldBoundaryValues(111,3, NULL);
// //   P.setFieldBoundaryValues(111,1, NULL);
// //   P.setAnalyticSource(1,4, "IVec*0.+JVec*0.+KVec*100.+LVec*(0.)");
//   P.setAnalyticSource(2,4, "IVec*0.+JVec*0.+KVec*0.+LVec*(1000.)");

  
//   ParaMEDMEM::MEDCouplingFieldDouble *field=build_field(osfilenameP[0].str(),"10*x",2);
// //    P.setAnalyticSource(2,1, "10*x");
//   P.setFieldSource(2,1, field);
  
  
  // system 2
  // MGFemusInit --------------------------------------------------------------
//    FEMUS P1;                                        // constructor
//   P1.init_param(*mgutils[1]);                      // init parameter
//   P1.init_fem(*mggeomel,*mgfemap);                 // init fem      
   // setting mesh ----------------------------------------
//   P1.setMedMesh(osfilenameP[1].str().c_str());     // set med-mesh
//   P1.setMesh();                                    // set MGmesh   
   // setting system ----------------------------------------
//   P1.setSystem(myproblemP);                        // set system
  
//   P1.init_interface(211,211,2,filenameP[1].str().c_str(),filenameP[0].str().c_str(),P,1); // init interfaces
//   P1.init_interface(212,211,1,filenameP[1].str().c_str());
//   P.init_interface(212,211,1,filenameP[0].str().c_str(),filenameP[1].str().c_str(),P1,1); // init interfaces
//   P1.init_interface(1,1, 2,filenameP[1].str().c_str());
//   P1.setAnalyticSource(1,3, "IVec*0.+JVec*200.+KVec*0.");
  
//   P1.init_interface(210,211,2,filenameP[1].str().c_str(),filenameP[0].str().c_str(),P,1); // init interfaces
//   P1.init_interface(1,  2,filenameP[1].str().c_str(),false);
//   P1.setAnalyticSource(1,1, "1.*sin(3.14/0.4* x) +3*sin(3.14/0.4* y)+ 1*sin(3.14/0.4* z)");

  // solving
  int    n_steps = mgutils[0]->get_par("nsteps");
  double      dt = mgutils[0]->get_par("dt");
  int print_step = mgutils[0]->get_par("printstep");
  int    itime_0  = mgutils[0]->get_par("itime");
  double time    = 0.;
  
  P.solve_setup(itime_0,time);                    // initial time loop (t=0)
//   P1.solve_setup(itime_0,time);                   // initial time loop (t=0)

//    ParaMEDMEM::MEDCouplingFieldDouble *bdy=NULL; 
  
  // transient loop with  n_steps (i0time ->i0time+ n_steps)
  for(int itime=itime_0; itime< itime_0 + n_steps; itime++) {

    
    P.solve_onestep(itime_0,itime,print_step,time,dt);    // solving P
    
//     bdy = P.getValuesOnBoundary(211, "NS0",3);          // take field from P problem
//     P1.setFieldBoundaryValues(211,3, bdy);                // store field inside P1 problem
//     P1.write_Boundary_value(211,"NS0",3);               // write inside x_old the boundary values

//     bdy = P.getValuesOnBoundary(210, "T",1);            // take field from P problem
//     P1.setFieldBoundaryValues(210,1, bdy);              // store field inside P1 problem
//     P1.write_Boundary_value(210,"T",1);                 // write inside x_old the boundary values
    
    
              
//     P1.solve_onestep(itime_0,itime,print_step,time,dt);   // solving P1
        
//     bdy = P1.getValuesOnBoundary(212, "NS0",1,3);          // take field from P problem
//     P.setFieldBoundaryValues(212,1, bdy);                // store field inside P1 problem
//     P.write_Boundary_value(212,"NS0",1,3);               // write inside x_old the boundary values
  }   // end time loop

  // end ======================================================================
  // --------------------------------------------------------------------------
  P.terminate(); 
//   P1.terminate();

  // clean --------------------------------------------------------------------
  mgutils.clear();
  delete dfe_q;  delete dfe_l;   delete dfe_k;  // delete   fem
  delete mggeomel; delete mgfemap;
  return 0;
}

#ifdef HAVE_MED
// ==========================================================================
  ParaMEDMEM::MEDCouplingFieldDouble *build_field(
    const std::string & path_mesh,
    const std::string & expression,
    int id_interface
  ){
  std::vector<std::string> vG(1);
  std::ostringstream s;
  s<< id_interface;
  vG[0] =s.str();
  ParaMEDMEM::MEDCouplingUMesh * support = 
  MEDLoader::ReadUMeshFromGroups(path_mesh.c_str(),"Mesh_1", 0,vG);
   
  
    ParaMEDMEM::TypeOfField type = ParaMEDMEM::ON_NODES;
  int dim=support->getSpaceDimension();
  std::vector<std::string> vars(dim);
  if(dim > 0)    vars[0] = "x";
  if(dim > 1)    vars[1] = "y";
  if(dim > 2)    vars[2] = "z";

  const ParaMEDMEM::MEDCouplingFieldDouble * field_nodes =
    support->fillFromAnalytic3(type, 1, vars, expression.c_str());

  int n_elements_med= support->getNumberOfCells(); // n MED-elements
//   const ParaMEDMEM::DataArrayDouble * d =support->getCoords();// Med-mesh coordinates
  int n_nodes_med = field_nodes->getNumberOfTuples();             //  MED-nodes
//   int dim = field_nodes->getNumberOfComponents();             //  MED dimension

  int n_nod_el=NDOF_FEM;
  double sum;

  
 ParaMEDMEM::MEDCouplingFieldDouble *field;
 field= ParaMEDMEM::MEDCouplingFieldDouble::New(ParaMEDMEM::ON_CELLS,ParaMEDMEM::NO_TIME);
 field->setMesh(support);

  ParaMEDMEM::DataArrayDouble *array;
  array=ParaMEDMEM::DataArrayDouble::New();
  array->alloc(n_elements_med,1);

  double *field_av = new double [n_elements_med]; // MED el centers
//   std::vector<std::pair<double, int > > dist_med;       // MED el distance
//   std::vector<std::pair<double, int > >::iterator itdist;

//   _elemID.clear();  // clear map elements: FEMus-mesh -> MED-mesh
  std::vector<int> nodes1; // element nodes

  // Computing the MED center (xyz_med) ***************************************
  for(int ielem=0; ielem<  n_elements_med; ielem++) {
    sum=0.; // zeros
    support->getNodeIdsOfCell(ielem, nodes1);  // element nodes
    for(int inode=0; inode< n_nod_el; inode++) {
//       for(int idim=0; idim<dim; idim++)
      sum += field_nodes->getIJ(nodes1[inode],0);
    } // end inode
    nodes1.clear(); // clear element node vector
//     for(int idim=0; idim<dim; idim++)

    field_av[ielem]=sum/n_nod_el;
  }// *************************************************************************
  
   std::copy(field_av,field_av+n_elements_med,array->getPointer());

  field->setArray(array);
  field->setName("pippo");
  field->checkCoherency();
  
  return field;
  }
  
#endif

