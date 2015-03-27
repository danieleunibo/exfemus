#ifndef __mgtimeloop_h__
#define __mgtimeloop_h__

// Forward class
class MGUtils;
class EquationSystemsExtendedM;


// ===============================================
//                  MGTimeLoop class
// ===============================================
/// Class doing generic time loop for the solution of the system
class MGTimeLoop  {

protected:
  // Data ---------------------------
  MGUtils &  _mgutils; ///< MGUtils vector of pointers
  EquationSystemsExtendedM &  _mgeqmap; ///< MGEquationsMap vector of pointers

public:
  ///@{ \name CONSTRUCTOR-DESTRUCTOR
  /// This function is the constructor
  MGTimeLoop(
    MGUtils & _mgutils,                             ///< MGUtils class (data)
    EquationSystemsExtendedM & _mgeqmap  ///< EquationSystemsExtendedM
  );
  /// This function is the class destructor
  ~MGTimeLoop() {};
  ///@}

//-----------------------------------------------------------------------------
  ///@{ \name FUNCTIONS FOR TIME LOOP
  /// This function setup a time loop (t=0)
  void transient_setup(
    const int & restart,         ///< restart iteration (0=no restart)  (in)
    int& t_in,                   ///< initial time iteration            (in)
    double&  time                ///< running time                      (in)
  );
/// This function runs a time loop
  void transient_loop(
    const int& t_in,               ///< initial time iteration          (in)
    double&  time                  ///< running time                    (in)
  );
/// This function runs  a time step
  void transient_onestep(
    const int  & t_in,                 ///< initial time iteration      (in)
    const int  & t_step,               ///< running time iteration      (in)
    const int  & print_step,           ///< print every                 (in)
    double     &  time,                ///< running time                (in)
    double     &  dt                   ///< step time                   (in)
  );
  /// This function runs  a time step
  void transient_non_linear_step(
    const int  & t_in,                 ///< initial time iteration      (in)
    const int  & t_step,               ///< running time iteration      (in)
    const int  & print_step,           ///< print every                 (in)
    double     &  time,                ///< running time                (in)
    double     &  dt                   ///< step time                   (in)
  );
  /// This function  prints a file time.****.xmf in RESU
  void transient_print_xmf(
    int t_init                          ///< initial time iteration     (in)
  );
///@}
};

#endif //  ---------  end header -------------------------
