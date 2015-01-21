
//VOF stuff

//in EquationsMap.C
#ifdef TWO_PHASE
///readCC
void MGCase::readCC(const uint flag_print) {

  std::ostringstream filename;
  filename <<"output/ccf." << setw(_ndig) << setfill('0') << flag_print << EXT_H5;
  hid_t  file_id = H5Fopen ( filename.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT );

  hid_t dtset = H5Dopen ( file_id, "/cres", H5P_DEFAULT );
  hid_t filespace = H5Dget_space(dtset);
  hsize_t dims[2]; H5Sget_simple_extent_dims(filespace, dims, NULL);
  
  double *ccf = ( new double[dims[0]] ); double *xf = ( new double[4*dims[0]] );
  double *yf = ( new double[4*dims[0]] );double *zf = ( new double[4*dims[0]] );
  
  hid_t status=H5Dread ( dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ccf );
  H5Dclose ( dtset );
  dtset = H5Dopen ( file_id, "/X1", H5P_DEFAULT );
  status=H5Dread ( dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,xf );
  H5Dclose ( dtset );
  dtset = H5Dopen ( file_id, "/X2", H5P_DEFAULT );
  status=H5Dread ( dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,yf );
  H5Dclose ( dtset );
  dtset = H5Dopen ( file_id, "/X3", H5P_DEFAULT );
  status=H5Dread ( dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,zf );
  H5Dclose ( dtset );
  
  _mgsys.get_mgcc()->read_fine(dims[0],ccf,xf,yf,zf);
  
  delete []ccf; delete []xf; delete []yf; delete []zf; 
  return;
}


/// It prints the format Xdmf for one time step of CC


void MGCase::print_xmfCC ( std::ofstream &out, const uint t_step,
                           const uint n_lines, const uint n_cells ) {
  //  Mesh ---- ---------------------------------------------
  const MGMesh &mgmesh=_mgsys.get_mesh();
//  uint n_nodes=mgmesh._NoNodes[_cslevel];
  uint n_elements=mgmesh._NoElements[_cslevel];
  const uint nvrt=4* ( DIMENSION-1 );
#if DIMENSION==2
  std::string btype="Polyline";
  std::string mtype="Quadrilateral";
#else
  std::string btype="Triangle";
  std::string mtype="Hexahedron";
#endif
// time parameters
  const double dt=_mgsys.get_par ( "dt" );

  // ============
  // coarse color function
  // ============
  out << "<Attribute Name=\"CC\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_elements*nvrt << "  " << 1 <<
  "\" Format=\"HDF\">  \n";
  out << "sol." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":CC" << "\n";
  out << "</DataItem>\n";
  out << "</Attribute>\n";
  // Grid Collection end
  out << "</Grid> \n";

  // ==================
  //  Interface
  // =================
  out << "<Grid Name=\"Interf\"> \n";
  out << "<Time Value =\"" << t_step*dt << "\" /> \n";
  // +++++ Topology ++++++++++
  out << "<Topology Type=\""<< btype << "\"   Dimensions=\""<<  n_lines <<    "\"> \n";
  out << "<DataStructure DataType=\"Int\" Dimensions=\""<< n_lines <<"  "<< DIMENSION <<
  "\" Format=\"HDF\">  \n";
  out  << "sol." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":" << "intconn" << "\n";
  out << "</DataStructure> \n" << "</Topology> \n";
  // +++++++  Geometry +++++++++++++++++
  out << "<Geometry Type=\"X_Y_Z\"> \n";
  for ( uint ix=1;ix<4;ix++ )
    {
      out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_lines*DIMENSION << "  " << 1 <<
      "\" Format=\"HDF\">  \n";
      out <<  "sol." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":" << "intX" <<ix<< "\n";
      out << "</DataStructure> \n";
    }
  out << " </Geometry>\n";
  // Grid Collection end
  out << "</Grid> \n";
  // =========================
  // cc on fine grid
  // =========================
  out << "<Grid Name=\"ccf\"> \n";
  out << "<Time Value =\"" << t_step*dt << "\" /> \n";
  // +++++ Topology ++++++++++
  out << "<Topology Type=\""<< mtype << "\"   Dimensions=\""<<  n_cells <<    "\"> \n";
  out << "<DataStructure DataType=\"Int\" Dimensions=\""<< n_cells <<"  "<< nvrt <<
  "\" Format=\"HDF\">  \n";
  out  << "ccf." << setw(_ndig) << setfill('0') << t_step << EXT_H5 <<":" << "conn" << "\n";
  out << "</DataStructure> \n" << "</Topology> \n";
  // +++++++  Geometry +++++++++++++++++
  out << "<Geometry Type=\"X_Y_Z\"> \n";
  for ( uint ix=1;ix<4;ix++ )
    {
      out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_cells*nvrt << "  " << 1 <<
      "\" Format=\"HDF\">  \n";
      out <<  "ccf." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":X" <<ix<<"\n";
      out << "</DataStructure> \n";
    }
  out << " </Geometry>\n";
  // ccf
  out  << "<Attribute Name=\"ccf\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
  out  << "<DataItem Dimensions=\""<< n_cells << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\"> \n ";
  out << "ccf." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":" << "ccf" << "\n";
  out << "</DataItem>\n" << "</Attribute>\n";

  return;
}



void MGCase::print_h5CC ( hid_t file,const uint flag_print,
                          uint *n_l_out,uint *n_c_out )  {
  //  Mesh ---- ---------------------------------------------
  const MGMesh &mgmesh=_mgsys.get_mesh();
//  uint n_nodes=mgmesh._NoNodes[_cslevel];
  const uint n_elements=mgmesh._NoElements[_cslevel];

  int nvrt= ( DIMENSION-1 ) *4;
  // ==========
  //  function cc
  // ===========
  double *ucc;  ucc=new double[n_elements* ( DIMENSION-1 ) *4];
  hsize_t dimsf[2];  dimsf[0] = n_elements* ( DIMENSION-1 ) *4;  dimsf[1] = 1;
  MGSolCC *mgcc=_mgsys.get_mgcc();
  mgcc->Write_xmf ( ucc );
  std::string name="CC";  print_Dhdf5(file,name,dimsf,ucc);
  delete []ucc;
  // ==============
  // interface
  // =============
  const uint n_lines=mgcc->get_nel();
  // line coordinates
  double *xcoord;
  xcoord=new double[n_lines*DIMENSION];
  double *ycoord;
  ycoord=new double[n_lines*DIMENSION];
  double *zcoord;
  zcoord=new double[n_lines*DIMENSION];
  int *tcon;
  tcon=new int[n_lines*DIMENSION];
  mgcc->get_int ( n_lines,tcon,xcoord,ycoord,zcoord);
  // coordinate datasets --------------------
  dimsf[0] = n_lines*DIMENSION;  dimsf[1] = 1;
  std::string name="intX1";
    print_Dhdf5(file,name,dimsf,xcoord);
//   dataspace = H5Screate_simple ( 2, dimsf, NULL );
//   dataset = H5Dcreate ( file, "intX1",  H5T_NATIVE_DOUBLE,
//                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
//                       H5P_DEFAULT, xcoord );
//   H5Dclose ( dataset );
   name="intX2";
    print_Dhdf5(file,name,dimsf,ycoord);
//   dataset = H5Dcreate ( file, "intX2",  H5T_NATIVE_DOUBLE,
//                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
//                       H5P_DEFAULT, ycoord );
//   H5Dclose ( dataset );
 name="intX3";
    print_Dhdf5(file,name,dimsf,zcoord);
//   dataset = H5Dcreate ( file, "intX3",  H5T_NATIVE_DOUBLE,
//                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
//                       H5P_DEFAULT, zcoord );
//   H5Dclose ( dataset );
//   H5Sclose ( dataspace );
  // line connectivity ------------------------
    dimsf[0] = n_lines;  dimsf[1] = DIMENSION;
   name="intconn";    print_Ihdf5(file,name,dimsf,tcon);
//   dataspace = H5Screate_simple ( 2, dimsf, NULL );
//   dataset = H5Dcreate ( file, "intconn", H5T_NATIVE_INT, dataspace,
//                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//   status = H5Dwrite ( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
//                       H5P_DEFAULT,tcon );
//   H5Dclose ( dataset );
  // Close/release resources.
//   H5Sclose ( dataspace );
  delete []tcon;  delete []xcoord;
  delete []ycoord;  delete []zcoord;

  // ==========
  //  ccf
  // ==========
  std::ostringstream filename;
  filename << "output/ccf." << setw(_ndig) << setfill('0') << flag_print << EXT_H5;
  hid_t file1= H5Fcreate ( filename.str().c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT );

  // ccf
  const uint n_cells=mgcc->get_nef();
  // coordinates
  xcoord=new double[n_cells*nvrt];
  ycoord=new double[n_cells*nvrt];
  zcoord=new double[n_cells*nvrt];
  double *ccc = ( new double[n_cells] ) ;
  double *cres = ( new double[n_cells] ) ;
  tcon=new int[n_cells*nvrt];
  mgcc->Writefine_xmf ( ccc,cres,xcoord,ycoord,zcoord );
  // coordinate datasets --------------
  dimsf[0] = n_cells*nvrt;  dimsf[1] = 1;
  
  name="X1";  print_Dhdf5(file1,name,dimsf,xcoord);
//   dataspace = H5Screate_simple ( 2, dimsf, NULL );
//   dataset = H5Dcreate ( file1, "X1",  H5T_NATIVE_DOUBLE,
//                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
//                       H5P_DEFAULT, xcoord );
//   H5Dclose ( dataset );
 name="X2";  print_Dhdf5(file1,name,dimsf,ycoord);
//   dataset = H5Dcreate ( file1, "X2",  H5T_NATIVE_DOUBLE,
//                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
//                       H5P_DEFAULT, ycoord );
//   H5Dclose ( dataset );
  name="X3";  print_Dhdf5(file1,name,dimsf,zcoord);
//   dataset = H5Dcreate ( file1, "X3",  H5T_NATIVE_DOUBLE,
//                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
//                       H5P_DEFAULT, zcoord );
//   H5Dclose ( dataset );
//   H5Sclose ( dataspace );
  // ccf --------------------------------
  dimsf[0] = n_cells;  dimsf[1] = 1;
//   dataspace = H5Screate_simple ( 2, dimsf, NULL );
  name="ccf";  print_Dhdf5(file1,name,dimsf,ccc);
//   dataset = H5Dcreate ( file1, "ccf",  H5T_NATIVE_DOUBLE,
//                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
//                       H5P_DEFAULT, ccc );
//   H5Dclose ( dataset );
name="cres";  print_Dhdf5(file1,name,dimsf,cres);
//   dataset = H5Dcreate ( file1, "cres",  H5T_NATIVE_DOUBLE,
//                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
//                       H5P_DEFAULT, cres );
//   H5Dclose ( dataset );		      
//   H5Sclose ( dataspace );
  // line connectivity ----------------
  for ( uint it=0;it<n_cells*nvrt;it++ ) tcon[it]=it;
  dimsf[0] = n_cells;  dimsf[1] = nvrt;
  name="conn";  print_Dhdf5(file1,name,dimsf,tcon);
//   dataspace = H5Screate_simple ( 2, dimsf, NULL );
//   dataset = H5Dcreate ( file1, "conn", H5T_NATIVE_INT, dataspace,
// //                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
//   status = H5Dwrite ( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
//                       H5P_DEFAULT,tcon );
//   H5Dclose ( dataset );
//   H5Sclose ( dataspace );
  H5Fclose ( file1 );
  delete []tcon;  delete []ccc;
   delete []cres;  delete []xcoord;
  delete []ycoord;  delete []zcoord;

  *n_l_out = n_lines;  *n_c_out = n_cells;
  return;
}

#endif