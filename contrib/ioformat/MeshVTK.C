#ifdef Vtk
// ==================================
// Print functions for Xdmf (return to file)
// functions:
// 	print
// 	print_mesh
// 	print_vol 
// =====================================


// ==============================================
/// It manages the printing in Xdmf format
// ==============================================
void MGMesh::print(const std::string & name,const unsigned int Level)
{ // ==============================================
  // Print Mesh Volume -----------
  // FEM order for VTK
  print_mesh(name,Level,Level);
  // print new connection mesh in xdmf format
  print_vol(name,Level);
  // Print Mesh boundary -----------
  // FEM order for VTK
  std::string namefile=name+"bd";
  print_mesh(namefile,Level,Level+_NoLevels);
  // print new connection boundary mesh in xdmf format
  print_vol(namefile,Level+_NoLevels);
}


// ===================================================
/// It prints the boundary mesh in VTK format
// ==============================================
void  MGMesh::print_bd(const std::string & name,const unsigned int Level)
{ // ===================================================
return;
}

// ===================================================
void  MGMesh::print_mesh(const std::string & name,const unsigned int Level,
			 const unsigned int idx_map)
{ // ===================================================
   
  switch(_type_FEM[idx_map]){
    case 27:
      _ord_FEM[0] = 0;_ord_FEM[1] = 1;_ord_FEM[2] = 2;_ord_FEM[3] = 3;_ord_FEM[4] = 4;_ord_FEM[5] = 5;_ord_FEM[6] = 6;_ord_FEM[7]=7;
      _ord_FEM[8] = 8;_ord_FEM[9] = 9;_ord_FEM[10]=10;_ord_FEM[11]=11;
      _ord_FEM[12]=16;_ord_FEM[13]=17;_ord_FEM[14]=18;_ord_FEM[15]=19;
      _ord_FEM[16]=12;_ord_FEM[17]=13;_ord_FEM[18]=14;_ord_FEM[19]=15;
      _ord_FEM[20]=24;_ord_FEM[21]=22;_ord_FEM[22]=21;_ord_FEM[23]=23;_ord_FEM[24]=20;_ord_FEM[25]=25;_ord_FEM[26]=26;
      _ord_FEM[27]=29;
    break;
    case 9:
      _ord_FEM[0]=0;_ord_FEM[1]=1;_ord_FEM[2]=2;_ord_FEM[3]=3;_ord_FEM[4]=4;_ord_FEM[5]=5;_ord_FEM[6]=6;_ord_FEM[7]=7;_ord_FEM[8]=8;
      _ord_FEM[27]= 28;
    break;
     case 3:
      _ord_FEM[0]=0;_ord_FEM[1]=1;
      _ord_FEM[27]=3;
    break;
    default:
      std::cout << "Error FEM Mesh"; exit(0);
    break;  
  }
return;
}

// ===================================================
/// It prints the volume Mesh in VTK format
// ==============================================
void  MGMesh::print_vol(const std::string & name,const unsigned int Level)
{ // ===================================================

  unsigned int nodeFEMvtk=_ord_FEM[27]; 
  unsigned int n_nodes = _NoNodes[Level];
  unsigned int nodeFEM=_type_FEM[Level];
  unsigned int n_elements = _NoElements[Level];
  std::string namefile=name+".vtu";
  std::ofstream out (namefile.c_str());std::cout << namefile.c_str();
       
  out << "<?xml version=\"1.0\"?>\n";
  out << "<VTKFile type=\"UnstructuredGrid\"  byte_type=\"LittleEndian\">\n";
  out << "<UnstructuredGrid>\n";

  // --------------  Body --------------
  out << "<Piece NumberOfPoints=\""<<n_nodes << "\" NumberOfCells=\"" << n_elements<< "\">\n";   
  // write the nodes --------------------------------------
  out << "<Points> \n";
  out <<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(unsigned int i=0;i<n_nodes;i++) out <<_xyz[i]<<"  "<<_xyz[i+n_nodes]<<"  "<<  
#if DIMENSION==3	    
      _xyz[i+2*n_nodes]<< "  ";
#else
      0.<< "  ";
#endif	    
  out << "\n";
  out << "</DataArray>\n" "</Points>\n";
  // write the connectivity
  out << "<Cells>\n";
  out <<"<DataArray type=\"Int32\" Name=\"connectivity\"  format=\"ascii\">\n";
  for(unsigned int el=0;el<n_elements;el++) 
    for(unsigned int ik=0;ik<nodeFEM;ik++)  
      out<< _conn_map[Level][_ord_FEM[ik]+el*nodeFEM]<<" ";
  out << "\n</DataArray>\n";
  out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
  int offsets = 0;
  for (unsigned int el = 0; el < n_elements; el++){
    offsets += nodeFEM;    out << offsets << " ";
  } 
  out << "\n</DataArray>\n";
  out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
  for (unsigned int el = 0; el < n_elements; el++){out << nodeFEMvtk << " ";}
  out << "\n</DataArray>\n"<< "</Cells>\n";
  //  ------------  end body ----------------------
  out << "</Piece>\n" << "</UnstructuredGrid>\n" << "</VTKFile>\n";
  out.close();
  return;
}

#endif

