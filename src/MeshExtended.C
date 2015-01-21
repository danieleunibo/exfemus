#include <iostream>
#include <iomanip>
#include <cstring>
#include <sstream>

#include "MeshExtended.h"
#include "MGGeomEl.h"
#include "MGUtils.h"


#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDFileMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#endif


// ===============================================================
// This is the constructor
MeshExtended::MeshExtended(
  const ParallelObjectM              & comm,
  MGUtils                            & mgutils,
  MGGeomEl                           & geomel
) : MGMesh(
    comm,
    mgutils,
    geomel
  )

  {

  read_bc_id(_NoLevels-1);
  read_mat_id(_NoLevels-1);

  return;
}
// =======================================================
MeshExtended::~MeshExtended() {
  _bc_id.clear(); _mat_id.clear(); 
}




void MeshExtended::read_bc_id(int Level) {

  std::string    input_dir = _mgutils._inout_dir;
  std::string    basemesh = _mgutils.get_file("BASEMESH");

  // Open an existing file. ---------------
  std::ostringstream meshname;  meshname << input_dir << basemesh << ".h5";
  std::cout << " Reading bc_id from= " <<  meshname.str() <<  std::endl;
  hid_t  file_id = H5Fopen(meshname.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

  // Reading  bc_id ---------------------------------------

  // Getting dataset
  std::ostringstream Name; Name << "NODES/COORD/BC";
  hid_t dtset = H5Dopen(file_id,Name.str().c_str()
#if HDF5_VERSIONM == 1812
                        ,H5P_DEFAULT
#endif
                       );
  hid_t filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  hsize_t dims[2];
  hid_t status  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status ==0) std::cerr << "SalomeIO::read dims not found";
//   _n_bc_id=dims[0];
  _bc_id.resize(dims[0]);

  // reading
  status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,&_bc_id[0]);
  
  
  //   int ielemn=0;
  for (int iproc=0; iproc<_n_subdom; iproc++) {
    int nel_b=_off_el[0][_NoLevels*iproc+Level];
    int nel_e=_off_el[0][Level+_NoLevels*iproc+1];
    for (int iel=nel_b; iel <nel_e; iel++) {
      
     
      for (int  inode=0; inode<NDOF_P; inode++)    {
        //get the global node number
        int el_conn = _el_map[0][iel*NDOF_FEM+inode];
        // get the element coordinades
        if (_bc_id[el_conn]>0)  _bc_id[el_conn] *=-1;
      }
  } // ---- end iel -------------------------
} // 2bB end interpolation over the fine mesh -
  
  return;
}

void MeshExtended::read_mat_id(int Level) {

  std::string    input_dir = _mgutils._inout_dir;
  std::string    basemesh = _mgutils.get_file("BASEMESH");

  // Open an existing file. ---------------
  std::ostringstream meshname;  meshname << input_dir << basemesh << ".h5";
  std::cout << " Reading mat_id from= " <<  meshname.str() <<  std::endl;
  hid_t  file_id = H5Fopen(meshname.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);

  // Reading  bc_id ---------------------------------------

  // Getting dataset
  std::ostringstream Name; Name << "ELEMS/SUB/MAT";
  hid_t dtset = H5Dopen(file_id,Name.str().c_str()
#if HDF5_VERSIONM == 1812
                        ,H5P_DEFAULT
#endif
                       );
  hid_t filespace = H5Dget_space(dtset);    /* Get filespace handle first. */
  hsize_t dims[2];
  hid_t status  = H5Sget_simple_extent_dims(filespace, dims, NULL);
  if(status <0) {std::cerr << " ******** SalomeIO::read dims not found"; abort();}
   _mat_id.resize(dims[0]);

  // reading
  status=H5Dread(dtset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,&_mat_id[0]);

  return;
}



#ifdef HAVE_MED
// ========================================================
// This function print the med file
void MeshExtended::print_med(int Level, std::string filename) {
  
  // med->libmesh map  second order only
#if ELTYPE==27
#if DIMENSION==3
  const unsigned int nodesinv[]   = {4,7,3,0,5,6,2,1,19,15,11,12,17,14,9,13,
                                     16,18,10,8,24,25,23,20,21,22,26};
  const unsigned int nodesinvbd[] = {3,0,1,2,7,4,5,6,8};
#else
  const unsigned int nodesinv[] = {3,0,1,2,7,4,5,6,8};
  const unsigned int nodesinvbd[] = {0,2,1};
#endif
#endif
#if ELTYPE==10
#if DIMENSION==3
  const unsigned int nodesinv[]   = {2,3,1,0,9,8,5,6,7,4};
  const unsigned int nodesinvbd[] = {1,2,0,4,5,3};
#else
  const unsigned int nodesinv[]   = {1,2,0,4,5,3};
  const unsigned int nodesinvbd[] = {0,2,1};
#endif
#endif

  int       el_neigh[NDOF_FEM];          // bd element connectivity
  int       sur_toply[NDOF_FEMB];
  int       el_conn[NDOF_FEM];
//   int       elb_conn[NDOF_FEMB];
  double    *xx_qnds=new double[_dim*NDOF_FEM];

  // setup  -------------------------------------------
//   Level=_NoLevels-1;
  int      n_nodes    =_NoNodes[Level];
  int      n_elements = _NoElements[0][Level];
  int       nodes_el   = _type_FEM[0];
  const int el_sides=  _GeomEl._n_sides[0];
  const int pt_sides_bd=NDOF_FEMB;
//
//   // init  -------------------------------------------
  std::cout << "n_nodes " << n_nodes << "n_elements "
            << n_elements << "nodes_el " << nodes_el << std::endl;
 int icount =0;  int ielemcount =0;
  int ibccount=0; int bd_n_elements =0;

//   // coords ------------------------------------------
  double * coord; coord = new double[n_nodes*_dim];
  for(int i=0; i<n_nodes; i++) {
    for(int idim=0; idim<_dim; idim++) {
      coord[i*_dim+idim]=_xyz[i+idim*n_nodes];
    }
  }

//
  int * conn_bd; conn_bd=new int [n_elements*el_sides*pt_sides_bd];
  int * conn; conn=new int [n_elements*nodes_el];
  std::vector<int> elem_bd_id2;
//   elem_id2[0]=0;
//   std::vector< int > asf; asf.resize(n_elements);
  for(int  iproc = 0; iproc <_n_subdom; iproc++) {
    const int nel_e = _off_el[0][Level+_NoLevels*iproc+1]; // start element
    const int nel_b = _off_el[0][Level+_NoLevels*iproc];   // stop element
    for(int el=0; el < (nel_e - nel_b); el++) {
//     for (int el =0 ;
//          el <_off_el[0][iproc*_NoLevels+Level+1]-
//          _off_el[0][iproc*_NoLevels+Level]; el++) {
      get_el_nod_conn(0,Level,el,el_conn,xx_qnds,iproc);
      get_el_neighbor(el_sides,0,Level,el,el_neigh,iproc);

      for(int  i = 0; i < nodes_el; i++) {
        conn[icount] = el_conn[nodesinv[i]];//el_map[0][el*nodes_el+i];
//          _proc_id[conn[icount]]= iproc;// map my-mesh->iproc
         icount++;
      }
      // element boundary
      for(int iside=0; iside< el_sides; iside++)  {
        if(el_neigh[iside] == -1){   // domain boundary
          int min=1000000;
          for(int idof=0; idof<pt_sides_bd; idof++) {                          //idof -> boundary node
            sur_toply[idof]=_GeomEl._surf_top[nodesinvbd[idof]+pt_sides_bd*iside]; //use map to find global node
            int idofb=sur_toply[idof];   //idofb -> element node
            conn_bd[ibccount]=el_conn[idofb];//_el_map[0][el*nodes_el+idofb];
            ibccount++;
//
            if(min>fabs(_bc_id[el_conn[idofb]]))  min=fabs(_bc_id[el_conn[idofb]]);
          }
          elem_bd_id2.push_back(min);
          bd_n_elements++;  // domain boundary found
        }
      }
//       _elem_id[el]=ielemcount+1; // map my-mesh->MEDcoupling

//       asf[ielemcount]=el+nel_b;

//       ielemcount++;
    }  // end for-el
  } // end for-iproc
  

  // MED mesh *************************************************
//
  // MEDCouplingUMesh mesh connectivity (volume)
  ParaMEDMEM::MEDCouplingUMesh *mesh1=ParaMEDMEM::MEDCouplingUMesh::New("Mesh_1",_dim);
  mesh1->allocateCells(n_elements);
  for(int  i = 0; i < n_elements; i++)
    mesh1->insertNextCell(MED_EL_TYPE,nodes_el,conn+i*nodes_el);
  mesh1->finishInsertingCells();

// MEDCouplingUMesh Mesh connectivity (boundary)
  ParaMEDMEM::MEDCouplingUMesh *mesh2=ParaMEDMEM::MEDCouplingUMesh::New("Mesh_1",DIMENSION-1);
  mesh2->allocateCells(bd_n_elements);
  for(int  i = 0; i < bd_n_elements; i++)
    mesh2->insertNextCell(MED_EL_BDTYPE,pt_sides_bd,conn_bd+i*pt_sides_bd);
  mesh2->finishInsertingCells();

  // coord (same node set for both meshes)
  ParaMEDMEM::DataArrayDouble *coordarr=ParaMEDMEM::DataArrayDouble::New();
  coordarr->alloc(n_nodes,_dim);
  std::copy(coord,coord+n_nodes*_dim,coordarr->getPointer());
  mesh1->setCoords(coordarr);
  mesh2->setCoords(coordarr);
// //
  // Setting MEDCouplingUMesh into MEDFileUMesh
  ParaMEDMEM::MEDFileUMesh *mm=ParaMEDMEM::MEDFileUMesh::New();
  mm->setName("Mesh_1");//name needed to be non empty
  mm->setDescription("Description Mesh_1");
  mm->setCoords(mesh1->getCoords());
  mm->setMeshAtLevel(0,mesh1,false);
  mm->setMeshAtLevel(-1,mesh2,false);

  //Volume Groups

  std::map <int,std::vector<int> > vol_group_elements;
//    int Level=_NoLevels-1;
  int count=0;
   for(int iproc=0; iproc<_n_subdom; iproc++) {
    int nel_b=_off_el[0][_NoLevels*iproc+Level];
    int nel_e=_off_el[0][Level+_NoLevels*iproc+1];
    for(int iel=nel_b; iel <nel_e; iel++) {
      vol_group_elements[_mat_id[iel]].push_back(count);
      count++;
    }// ---- end iel -------------------------
  } // ************************************************************************
  

  int n_vol_group=vol_group_elements.size();
  std::vector<const ParaMEDMEM::DataArrayInt *> gr_vol(n_vol_group);
  ParaMEDMEM::DataArrayInt ** g_vol=new ParaMEDMEM::DataArrayInt *[n_vol_group];

  int js=0;
  // defining the vol group data to store
  std::map<int,std::vector<int> >::iterator it_vol;
  for(it_vol=vol_group_elements.begin(); it_vol!=vol_group_elements.end(); ++it_vol) {
    int igroup=it_vol->first;  int is=it_vol->second.size();
    g_vol[js]=ParaMEDMEM::DataArrayInt::New();
    g_vol[js]->alloc(is,1);
    std::ostringstream name_p; name_p<< igroup;
    g_vol[js]->setName(name_p.str().c_str());
    int *val1=new int[is];
    for(int iv=0; iv<is; iv++)  val1[iv]=it_vol->second[iv];
    std::copy(val1,val1+is,g_vol[js]->getPointer());
    delete []val1;
    gr_vol[js]=g_vol[js] ;
    js++;
  }
  // inserting  the volume groups into the med-file
  mm->setGroupsAtLevel(0,gr_vol,false);


  // Boundary group ******************************************
  // Finding bd_group and  bd_group_elements map
  std::map <int,std::vector<int> > bd_group_elements;
//   std::map <int, int> bd_group;
  for(int i=0; i<bd_n_elements ; i++) {
//     bd_group [elem_bd_id2[i]]++;
    bd_group_elements[elem_bd_id2[i]].push_back(i);
  }
  int n_bd_group=bd_group_elements.size();
  // group vector
  std::vector<const ParaMEDMEM::DataArrayInt *> gr_bd(n_bd_group);
  ParaMEDMEM::DataArrayInt ** g_bd=new ParaMEDMEM::DataArrayInt *[n_bd_group];

  js=0;
  // defining the  group data to store
  std::map<int,std::vector<int> >::iterator it;
  for(it=bd_group_elements.begin(); it!=bd_group_elements.end(); ++it) {
    int igroup=it->first;  int is=it->second.size();
//     std::cout << igroup << "\n";
    g_bd[js]=ParaMEDMEM::DataArrayInt::New();
    g_bd[js]->alloc(is,1);
    std::ostringstream name_p; name_p<< igroup;
    g_bd[js]->setName(name_p.str().c_str());
    int *valb1=new int[is]; /*int icount=0;*/
    for(int iv=0; iv<is; iv++)  valb1[iv]=it->second[iv];
    std::copy(valb1,valb1+is,g_bd[js]->getPointer());
    delete []valb1;
    gr_bd[js]=g_bd[js] ;
    js++;
  }
  // insert the boundary groups into the med-file
  mm->setGroupsAtLevel(-1,gr_bd,false);

  // Printing Group and Family
  std::cout << "\n \n =================================== ";
  std::cout << "\n Group Names -> Families : \n";
  std::map<std::string, std::vector<std::string> > a=(mm->getGroupInfo());
  std::map<std::string, std::vector<std::string> >::iterator ita;
  for(ita=a.begin();  ita!=a.end(); ++ita) {
    std::string igroup=ita->first;  int is=ita->second.size();
    std::cout << "\n "<< igroup << "-> ";
    for(int i=0; i<is; i++) {
      std::string a2 =ita->second[i];
      std::cout << a2 << "  ";
    }
  }
  std::cout << "\n \n =================================== ";
  std::cout << " \n Families -> Groups id: \n";
  std::map<std::string,int> fama=mm->getFamilyInfo();
  std::map<std::string, int >::iterator itfa;
  for(itfa=fama.begin();  itfa!=fama.end(); ++itfa) {
    std::string igroup=itfa->first;
    std::cout << "\n "<< igroup << "-> ";
    int a2 =itfa->second;
    std::cout << a2 << "  ";
  }
  mm->write(filename.c_str(),2);

  // clean
  delete []coord; delete []xx_qnds; fama.clear(); a.clear();
   coordarr->decrRef(); mesh1->decrRef(); mesh2->decrRef();
   mm->decrRef();
   delete []conn_bd;   delete []conn;
   elem_bd_id2.clear();
   vol_group_elements.clear();bd_group_elements.clear();
//    for (int idel=0;idel<n_bd_group;idel++){
//      g_bd[idel]->decrRef(); gr_bd[idel]->decrRef();
//    }
//    for (int idel=0;idel<n_vol_group;idel++){
//      g_vol[idel]->decrRef(); gr_vol[idel]->decrRef();
//    }
return;
}
#endif //HAVE MED defined



