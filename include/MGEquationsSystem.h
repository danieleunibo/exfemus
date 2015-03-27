#ifndef __mgequationsmap_h__
#define __mgequationsmap_h__

#include <map>
using namespace std;
#include <string>

#include "MGSolverBase.h"
#include "MGSystem.h"
class MGUtils;
// class MGSystem;
class MGMesh;
class MGFEMap;

///This class contains the system map 

class MGEquationsSystem: public MGSystem {  

  protected:
    // data --------------------------------
    map<string,MGSolBase*> _equations;   // system map
public:  
//     MGUtils& _mgutils;  ///< MGUtils class  pointer
//     MGSystem& _mgphys;  ///< MGSystem class pointer
//     MGMesh&   _mgmesh;  ///< MGMesh  class  pointer
    MGFEMap& _mgfemap;  ///< MGFEMap class  pointer
    
  // Constructor / Destructor -----------------------------------
  ///@{ \name CONSTRUCTOR-DESTRUCTOR
    MGEquationsSystem( MGUtils& mgutils_in,
// 				   MGSystem& mgphys_in,
				   MGMesh& mgmesh_in,
				   MGFEMap& mgfemap_in, 
                                   int np_data,
                                   int ncell_data
                     );
  ~MGEquationsSystem();
  ///@}
  
  ///@{ \name INIT AND CLEAN
  virtual void init(const std::vector<FIELDS> & pbName); ///< Initialize
  
  /// Clean all substructures
  void clean(); 
  ///@}
//---------------------------------------------------------------------------------------------
  ///@{ \name EQUATIONS GET/SET
  inline            void  set_eqs(MGSolBase* value)          {_equations.insert(make_pair(value->_eqname,value));}
  inline       MGSolBase* get_eqs(const string & name)       {return _equations.find(name)->second;}
  inline const MGSolBase* get_eqs(const string & name) const {return _equations.find(name)->second;}
  ///@}
//-----------------------------------------------------------------------------------------------
  ///@{ \name ITERATORS FOR EQUATION MAP
  typedef std::map<string, MGSolBase*>::iterator iterator;
  typedef std::map<string, MGSolBase*>::const_iterator const_iterator;

  inline iterator       begin()       { return _equations.begin();}
  inline iterator         end()       { return _equations.end();}
  inline const_iterator begin() const { return _equations.begin();}
  inline const_iterator   end() const { return _equations.end();}
  ///@}

  void setDofBcOpIc() ;
  void eqnmap_timestep_loop(const double time, const int delta_t_step_in);
  void eqnmap_timestep_loop_nonlinear(const double time, const int delta_t_step_in);
  void eqnmap_timestep_loop_control(const double time, const int delta_t_step_in); //for control problem
//------------------------------------------------------------------------------------------------
  ///@{ \name READ-PRINT FUNCTIONS 
  void print_soln(const int t_step); ///< Print solution
  void print_mesh_data(double vect_data[]);///< Print data from to mesh class (2-mesh code) 
  void print_case(const int t_init); ///< Print ic and bc
  
  /// Read solution
  void read_soln(const int t_step);  
  ///@}
#ifdef TWO_PHASE
  void  readCC(const int t_init);
#endif
//-------------------------------------------------------------------------------------------------
private:
  ///@{ \name PRINT XMF/H5 FUNCTIONS
 void print_soln_xmf(const int t_step, int n_l_out,int n_c_out);
 void print_soln_h5(const int t_flag);
 
 void print_case_xmf(const int t_init,const int n_lines,const int n_lin);
 void print_case_h5(const int t_init);
 ///@}
//  #ifdef TWO_PHASE
//  /// Print xmf file CC solution
//   void print_xmfCC(std::ofstream & out,const int t_init, const uint n_lines,const uint n_cells);  
//  /// Print CC solution in hdf5 format
//   void print_h5CC(hid_t file,const uint flag_print,int *n_l_out,uint *n_c_out); 
// #endif

};

#endif
