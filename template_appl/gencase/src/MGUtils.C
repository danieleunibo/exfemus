// std libraries ---------------------------------
#include <iomanip>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
// config includes --------------------------------
#include "Solverlib_conf.h"           // petsc conf

// local includes ---------------------------------
#include "MGUtils.h"


// ============================================================================
// ================ MGUtils Class functions ===================================
// ============================================================================
// ============================================================================
/// Constructor
MGUtils::MGUtils(
) {// =========================================================================
  read();
  _user_dir = _femus_dir +  "/USER_APPL/";
  _app_dir= _user_dir + _myapp_name;
  _inout_dir = _app_dir + "/RESU/";
  _mesh_dir = _user_dir + "/MESH/";
  _data_dir = _app_dir + "/DATA/";
  _fem_dir = _femus_dir + "/fem/";
  _contrib_dir=_femus_dir + "/contrib/";
  read_par();       // read parameters
  print();
  print_par();      // print in console parameters
}

// ============================================================================
// Constructor
MGUtils::MGUtils(
  int a  // mesh number <-
) {// =========================================================================

  std::ostringstream loc_par_file;
  loc_par_file << "/DATA/param_files_msh"<< a <<".in";
  const std::string & name=loc_par_file.str();
  read(name);
  _user_dir = _femus_dir +  "/USER_APPL/";
  _app_dir= _user_dir + _myapp_name;
  _inout_dir = _app_dir + "/RESU/";
  _mesh_dir = _user_dir + "/MESH/";
  _data_dir = _app_dir + "/DATA/";
  _fem_dir = _femus_dir + "/fem/";
  _contrib_dir=_femus_dir + "/contrib/";
  read_par();       // read parameters
  print();
  print_par();      // print in console parameters
}

// ============================================================================
/// This function reads all the  parameters
void MGUtils::read_par(
) {// =========================================================================

  // file names ------------------
  std::string    baseparutils = get_file("BASEPARUTILS");
  std::ostringstream filename;  filename << _data_dir << baseparutils;
  std::ifstream fin;  fin.open(filename.str().c_str());
  std::string buf="";
  double value;
  std::cout << "Init Reading" << std::endl;

  if(fin.is_open()) {
    while(buf != "/")  fin >> buf;
    fin >> buf;
    while(buf != "/") {
      if(buf == "#") getline(fin, buf);  // comment line
      else {
        fin >> value;    // set new parameter
        set_par(buf,value);
      }
      fin >> buf;// std::cerr <<buf.c_str() << "\n ";
    }
  } else {
    std::cerr << "MGUtils::read_par: no parameter file found" << std::endl;
    abort();
  }

  std::cout << "End Reading file " <<  filename.str() << std::endl;
  return;
}
// ============================================================================
/// This function reads the file names
void MGUtils::read(
  const std::string & name_file_in  // filename <-
) {//==========================================================================

  // read femus dir from shell -----------------------
  _femus_dir = getenv("FEMUS_DIR");
  if(_femus_dir == "") {
    std::cout << " Set FEMUS_DIR in your shell environment" << std::endl;
    abort();
  }
  _myapp_name = getenv("FM_MYAPP");
  if(_myapp_name == "") {
    std::cout << " Set MYAPP in your shell environment" << std::endl;
    abort();
  }

#ifdef PRINT_INFO
  std::cout << " femus_dir is " << _femus_dir << std::endl;
  std::cout << " myapp_name is " << _myapp_name << std::endl;
#endif

  // read param_files -----------------------------
  std::ostringstream filename;
  filename << _femus_dir <<  "/USER_APPL/"<< _myapp_name  <<"/" << name_file_in;

  std::cout << filename.str().c_str() <<" ";
  std::ifstream fin(filename.str().c_str());
  std::string buf="";
  std::string value;

  if(fin == NULL) {
    std::cerr << " MGFiles::read: no parameter file found" << std::endl;
    abort();
  }
  while(!fin.eof()) {
    fin >> buf;
    if(buf == "#")  fin.ignore(200,'\n');
    else { fin >> value;  set(buf,value);} // set new parameter
  }
  // cleaning and check
  fin.close();
return;
}

// ===========================================================================
/// This file prints the file name map
void MGUtils::print(
) {// =========================================================================
  std::cout << "\n ================================================ ";
  std::cout << "\n Class MGFiles: "<< (int)_mgfiles.size()
            << " file names: " << std::endl;
  std::map<std::string,std::string>::const_iterator pos=_mgfiles.begin();
  std::map<std::string,std::string>::const_iterator pos_e=_mgfiles.end();
  for(; pos!=pos_e; pos++) {
    std::string  name=pos->first; // get the string
    std::string value=pos->second;// get the name
    std::cout  <<" "<< std::left << std::setw(15)
               << name << " = " << value << std::endl;
  }
  return;
}

// ============================================================================
/// This function prints all the  parameters to stdout
void MGUtils::print_par(
) const {// ===================================================================
  std::cout << "\n ============================= ";
  std::cout << "\n   MGUtils: " << (int)_param_utils.size() << " parameters: \n";

  std::map<std::string,double>::const_iterator   pos=_param_utils.begin();
  std::map<std::string,double>::const_iterator pos_e=_param_utils.end();
  for(; pos!=pos_e; pos++) {
    std::string name=pos->first; // get name
    double value=pos->second;    // get value
    std::cout << "  " << std::left << std::setw(15) << name << " = "
              << std::setprecision(12) << value << std::endl;
  }
  return;
}

// ============================================================================
hid_t MGUtils::read_Dhdf5(
  hid_t file,                // file id        <-
  const std::string & name,  // dataset name   <-
  double data[]              // data to store  <-
) {// =========================================================================

  hid_t  dataset = H5Dopen(file,name.c_str()
#if HDF5_VERSIONM == 1812
                           , H5P_DEFAULT
#endif
                          );
  hid_t status=H5Dread(dataset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}

// ============================================================================
hid_t MGUtils::print_Dhdf5(
  hid_t file,const std::string & name,    // file name            <-
  hsize_t dimsf[],                        // dataspace dimension  <-
  double data[]                           // data to store        <-
) { // ========================================================================

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_DOUBLE,
                            dataspace, H5P_DEFAULT
#if HDF5_VERSIONM == 1812
                            , H5P_DEFAULT
#endif
                           );
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}

// ============================================================================
/// Print int data into dhdf5 file
hid_t MGUtils::print_Ihdf5(
  hid_t file,const std::string & name,     // file name            <-
  hsize_t dimsf[],                         // dataspace dimension  <-
  int data[]                               // data to store        <-
) { // ========================================================================

  hid_t dataspace = H5Screate_simple(2,dimsf, NULL);
  hid_t dataset = H5Dcreate(file,name.c_str(),H5T_NATIVE_INT,dataspace, H5P_DEFAULT
#if HDF5_VERSIONM == 1812
                            ,H5P_DEFAULT,H5P_DEFAULT
#endif
                           );
  hid_t  status = H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,H5P_DEFAULT,data);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  return status;
}

// ============================================================================
hid_t MGUtils::read_Ihdf5(
  hid_t file,                    // file id        <-
  const std::string & name,      // dataset name   <-
  int data[]                     // data to store  <-
) {// =========================================================================
  hid_t  dataset = H5Dopen(file,name.c_str()
#if HDF5_VERSIONM == 1812
                           ,H5P_DEFAULT
#endif
                          );
  hid_t status=H5Dread(dataset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,H5P_DEFAULT,data);
  H5Dclose(dataset);
  return status;
}
