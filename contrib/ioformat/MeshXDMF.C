#ifdef XDMF
// ==================================
// Print functions for Xdmf (return to file)
// functions:
// 	print
// 	print_mesh
// 	print_vol 
// 	print_bd
// =====================================

// ==============================================
/// It manages the printing in Xdmf format
// ==============================================
void MGMesh::print(const std::string & name,const unsigned int Level)
{ // ==============================================
  // print new connection mesh in hdf5 format
  print_mesh(name,Level,Level); 
  // print new connection boundary mesh in hdf5 format
  std::string namebd=name+"bd";
  print_mesh(namebd,Level,Level+_NoLevels);
  // print new connection mesh in xdmf format
  print_vol(name,Level);
  // print new connection boundary mesh in xdmf format
  print_bd(name,Level);
}

// ==============================================
/// It prints the volume Mesh (connectivity) in Xdmf format
// ==============================================
void MGMesh::print_vol(const std::string & name,const unsigned int Level)
{ // ==============================================
  const unsigned int n_nodes = _NoNodes[Level];
  const unsigned int n_el = _NoElements[Level];
  
  std::string namefile=name+".xmf";  std::ofstream out (namefile.c_str());
  int nvrt;std::string mtype;
#if DIMENSION==2
  nvrt = 4; mtype = "Quadrilateral";
#else
  nvrt = 8; mtype = "Hexahedron";
#endif
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" \n";
//   out << "[ <!ENTITY HeavyData \"mesh.h5 \"> ] ";  
  out << ">\n";
  out << " \n";
  out << "<Xdmf> \n" << "<Domain> \n"<< "<Grid Name=\"Mesh\"> \n";
  // Topology 
  out << "<Topology Type=\"" << mtype << "\"         Dimensions=\"" <<
    n_el*nvrt << "\"> \n";
  out << "<DataStructure DataType=\"Int\" Dimensions=\" " << n_el*
    nvrt << "  " << nvrt << "\" Format=\"HDF\">  \n";
  out << "meshxmf.h5:MSH"<<  Level << "CONN \n";  
  out << "</DataStructure> \n" << "</Topology> \n";
  // Geometry
  out << "<Geometry Type=\"X_Y_Z\"> \n";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" <<
  n_nodes << "  " << 1 << "\" Format=\"HDF\">  \n";
  out << "../data_in/mesh.h5:X1 \n" << "</DataStructure> \n";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" <<
  n_nodes << "  " << 1 << "\" Format=\"HDF\">  \n" << "../data_in/mesh.h5:X2 \n";
  out << "</DataStructure> \n";
  out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" <<
  n_nodes << "  " << 1 << "\" Format=\"HDF\">  \n";
  out << "../data_in/mesh.h5:X3 \n" << "</DataStructure> \n";

  out << " </Geometry>\n"<< "</Grid> \n" << "</Domain> \n" << "</Xdmf> \n";
  out.close ();
return;
}

// ==============================================
/// It prints the Xdmf fle for boundary visualization
// ==============================================
void MGMesh::print_bd(const std::string & name,const unsigned int Level)
{ // ==========================================
  const unsigned int n_nodes = _NoNodes[Level];
  const unsigned int n_el = _NoElements[Level];
  const unsigned int n_elb = _NoElements[_NoLevels + Level];

  std::string namefile=name+"bd"+".xmf";
  std::ofstream outb (namefile.c_str());
  int nvrtb;std::string top_type;
  
#if DIMENSION==2
  nvrtb = 2; top_type="Polyline";
#else
  nvrtb = 4; top_type="Quadrilateral";
#endif
    
  outb << "<?xml version=\"1.0\" ?> \n";
  outb << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" \n";
//   outb << "[ <!ENTITY HeavyData \"mesh.h5 \"> ] ";
  outb << ">\n";
  
  outb << "<Xdmf> \n" << "<Domain> \n" << "<Grid Name=\"Boundary\"> \n";
  
  // Topology 
  outb << "<Topology Type=\""<< top_type.c_str() <<"\"    Dimensions=\"" << n_elb*nvrtb 
	  << "\"> \n"<< "<DataStructure DataType=\"Int\" Dimensions=\" " 
	   << n_elb*nvrtb << "  " << nvrtb << "\" Format=\"HDF\">  \n";
  outb << "meshbdxmf.h5:MSH"<<_NoLevels + Level << "CONN \n"; 
  outb << "</DataStructure> \n" << "</Topology> \n";
  // Geometry
  outb << "<Geometry Type=\"X_Y_Z\"> \n";
  outb << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
	   << n_nodes << "  " << 1 << "\" Format=\"HDF\">  \n";
  outb << "../data_in/mesh.h5:X1 \n" << "</DataStructure> \n";
  outb << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
	   << n_nodes << "  " << 1 << "\" Format=\"HDF\">  \n";
  outb << "../data_in/mesh.h5:X2 \n"<< "</DataStructure> \n";
  outb << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\""
	   <<  n_nodes << "  " << 1 << "\" Format=\"HDF\">  \n";
  outb << "../data_in/mesh.h5:X3 \n"<< "</DataStructure> \n";
  outb << " </Geometry>\n";
  
  outb << "</Grid> \n" << "</Domain> \n" << "</Xdmf> \n";
  outb.close ();
  return;
}

// =========================================
/// It prints the connectivity in Xmdf format
/// The changes are only for visualization of quadratic FEM
// =========================================
 void MGMesh::print_mesh(const std::string & name,const unsigned int Level,
			 const unsigned int indx_mesh)
{  // ======================================
  unsigned int mode = _type_FEM[indx_mesh];
  unsigned int n_elements = _NoElements[indx_mesh];

  int conn[8][8];  unsigned int icount=0;
  unsigned int *gl_conn;
  
  switch(mode){
#if DIMENSION==2
    // -----------------------------------
    case 9:	// Quad 9  0-4-1-5-2-6-3-7-8
      gl_conn=new unsigned int[n_elements*4*4];
      conn[0][0] = 0;  conn[0][1] = 4;   conn[0][2] = 8; conn[0][3] = 7;// quad4  0-4-8-7
      conn[1][0] = 4;  conn[1][1] = 1;   conn[1][2] = 5; conn[1][3] = 8;// quad4  4-1-5-8
      conn[2][0] = 8;  conn[2][1] = 5;   conn[2][2] = 2; conn[2][3] = 6;// quad4  8-5-2-6
      conn[3][0] = 7;  conn[3][1] = 8;   conn[3][2] = 6; conn[3][3] = 3;// quad4  7-8-6-3
   
      for (unsigned int el = 0; el <n_elements; el++) {
	for (unsigned int se = 0; se < 4; se++) {
	  for (unsigned int i = 0; i < 4; i++) {
	    gl_conn[icount] = _conn_map[indx_mesh][el*mode+conn[se][i]]; icount++;
	  }
	}
      }
    break;
    // --------------------------------------
    case 3: // boundary edge 2 linear 0-2-1
      gl_conn=new unsigned int[n_elements*2*2];
      connb[0][0] = 0; connb[0][1] = 2;		// element 0-2
      connb[1][0] = 2; connb[1][1] = 1;		// element 1-2
      for (unsigned int el = 0; el < n_elements; el++) {
	for (unsigned int se = 0; se < 2; se++) {
	  for (unsigned int i = 0; i < 2; i++) {
	    gl_conn[icount]=_conn_map[indx_mesh][el*mode+conn[se][i]];
	    icount++;
	  }
	}
      } 
    break;
#else
  // ----------------------
  case  27: //  Hex 27 (8 Hex8)
    gl_conn=new unsigned int[n_elements*8*8];
    conn[0][0] = 0;  conn[0][1] = 8;  conn[0][2] = 20; conn[0][3] = 11; 
    conn[0][4] = 12; conn[0][5] = 21; conn[0][6] = 26; conn[0][7] = 24;  
    conn[1][0] = 8;  conn[1][1] = 1;  conn[1][2] = 9;  conn[1][3] = 20;
    conn[1][4] = 21; conn[1][5] = 13; conn[1][6] = 22; conn[1][7] = 26; 
    conn[2][0] = 11; conn[2][1] = 20; conn[2][2] = 10; conn[2][3] = 3;   
    conn[2][4] = 24; conn[2][5] = 26; conn[2][6] = 23; conn[2][7] = 15;
    conn[3][0] = 20; conn[3][1] = 9;  conn[3][2] = 2;  conn[3][3] = 10; 
    conn[3][4] = 26; conn[3][5] = 22; conn[3][6] = 14; conn[3][7] = 23;  
    conn[4][0] = 12; conn[4][1] = 21;  conn[4][2] = 26; conn[4][3] = 24;
    conn[4][4] = 4;  conn[4][5] = 16;  conn[4][6] = 25; conn[4][7] = 19; 
    conn[5][0] = 21;  conn[5][1] = 13; conn[5][2] = 22; conn[5][3] = 26;
    conn[5][4] = 16; conn[5][5] = 5;  conn[5][6] = 17; conn[5][7] = 25;
    conn[6][0] = 24; conn[6][1] = 26; conn[6][2] = 23; conn[6][3] = 15;
    conn[6][4] = 19; conn[6][5] = 25; conn[6][6] = 18; conn[6][7] = 7;
    conn[7][0] = 26; conn[7][1] = 22; conn[7][2] = 14; conn[7][3] = 23;
    conn[7][4] = 25; conn[7][5] = 17; conn[7][6] = 6;  conn[7][7] = 18;
    
    for (unsigned int el = 0; el < n_elements; el++) {
      for (unsigned int se = 0; se < 8; se++) {
	for (unsigned int i = 0; i < 8; i++) {
	  gl_conn[icount]=_conn_map[indx_mesh][el*mode+conn[se][i]];
	  icount++;
	}
      }
    } 
    break;
  // ---------------------------------------  
  case 9:  // Quad9 elements ( 4 Quad4)
    gl_conn=new unsigned int[n_elements*4*4];
    conn[0][0] = 0; conn[0][1] =4; conn[0][2] = 8; conn[0][3] = 7; 
    conn[1][0] = 4; conn[1][1] =1; conn[1][2] = 5; conn[1][3] = 8;  
    conn[2][0]= 8;  conn[2][1] =5; conn[2][2] = 2; conn[2][3] = 6;
    conn[3][0] = 7; conn[3][1] = 8;conn[3][2] = 6; conn[3][3] = 3;
    
    for (unsigned int el = 0; el <n_elements; el++) {
      for (unsigned int se = 0; se < 4; se++) {
	for (unsigned int i = 0; i < 4; i++) {
	  gl_conn[icount] = _conn_map[indx_mesh][el*mode+conn[se][i]]; icount++;
	}
      }
    }
  break;
#endif 
  // -----------------------------------------
  default:   // interior 3D
    gl_conn=new unsigned int[n_elements*mode];
    for (unsigned int el = 0; el <n_elements; el++) {
      for (unsigned int i = 0; i < mode; i++) {
	gl_conn[i+el*mode]=_conn_map[indx_mesh][i+el*mode];
      }
    }
    icount=n_elements*mode;
  break;
  }
  // storage in hf5 (Xdmf)
  std::string namefile=name+"xmf.h5"; std::cout << namefile << std::endl;
 
  hid_t file = H5Fcreate (namefile.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);	
  hsize_t dimsf[2]; dimsf[0] =icount;  dimsf[1] = 1;
   
  hid_t dtsp = H5Screate_simple (2, dimsf, NULL);
  std::ostringstream substor; substor << "MSH"<< indx_mesh << "CONN";
  hid_t dtset = H5Dcreate(file,substor.str().c_str(),H5T_NATIVE_INT,dtsp,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  H5Dwrite(dtset,H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT, gl_conn);	   
  H5Sclose(dtsp);H5Dclose (dtset);delete[]gl_conn;
  
  H5Fclose(file);
return;
}
#endif