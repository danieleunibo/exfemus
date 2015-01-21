// ==================================================================
//                  Class MGTimeLoop
// ==================================================================
// lib include
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>

// config include
// #include "FilePaths_conf.h"
#include "Printinfo_conf.h"
#include "Equations_conf.h"

// local include
#include "MGTimeLoop.h"

// class include
#include "MGUtils.h"
#include "EquationSystemsExtendedM.h"
#include "MGMesh.h"

// =============================================================================
/// Constructor
MGTimeLoop::MGTimeLoop(
  MGUtils & mgutils_in,  
  EquationSystemsExtendedM & mgeqmap_in
):
  _mgutils(mgutils_in),_mgeqmap(mgeqmap_in)
{// ===========================================================================
  
}

// ============================================================================
///This function does the  start or the restart
void MGTimeLoop::transient_setup(
    const int & restart,         ///< restart iteration (0=no restart)  (in)
    int& t_in,             ///< initial time iteration                  (in)
    double&  time                ///< running time                      (in)
)  {// ========================================================================

  if(restart != 0)  { _mgeqmap.read_soln(restart);}
  else   {
    _mgeqmap.print_soln(restart);  //print initial step sol.0.h5 and sol.0.xmf
    if(_mgeqmap._mgmesh._iproc == 0) _mgeqmap._mgmesh.write_c(restart);
  }
  // print --------------------------------------
  _mgeqmap.print_case(restart);   //print caseN.xmf&h5 = IC + BC flags
  this->transient_print_xmf(restart); //print timeN.xmf

  // restart -------------------------------------------
  if(restart) {
    double restartime=_mgutils.get_par("restartime");     // read  from file
    t_in=restart; time +=restartime; // set initial time
    std::cout << "\n *+*+* MGTimeLoop::transient_setup: RESTART  "
              << restart << " time= " << restartime << std::endl;
  }
  return;
}
// ============================================================================
/// This function controls the transient loop
void MGTimeLoop::transient_loop(
    const int& t_in,               ///< initial time iteration          (in)
    double&  time                  ///< running time                    (in)
)  { //========================================================================

  //  parameters
  int    n_steps = _mgutils.get_par("nsteps");
  double      dt = _mgutils.get_par("dt");
  int print_step = _mgutils.get_par("printstep");
  // transient loop
  for(int t_step= t_in; t_step< t_in + n_steps; t_step++) {
    transient_onestep(t_in,t_step,print_step,time,dt);  ///< step time  
  }   // end time loop
  return;
}


// ===========================================================================
/// This function controls the transient loop
void MGTimeLoop::transient_onestep(
    const int  & t_in,                 ///< initial time iteration      (in)
    const int  & t_step,               ///< running time iteration      (in)
    const int  & print_step,           ///< print every                 (in)
    double     &  time,                ///< running time                (in)
    double     &  dt                   ///< step time                   (in)  
)  { // =====================================================================

  // A Soving the system ****************************************************
#if PRINT_TIME==1 // only for cpu time check ----------------------------------
  std::clock_t  start_time=std::clock();
#endif   // ------------------------------------------------------------------- 
  std::cout << "\n  ** Solving time step " << t_step+1
            << ", time = " << time + dt << " ***" << std::endl;
  _mgeqmap.eqnmap_timestep_loop(time, t_step-t_in); // solve one step
#if PRINT_TIME==1 // only for cpu time check ----------------------------------
  std::clock_t    end_time=std::clock();
#endif            // ----------------------------------------------------------
  // B) print iteration t_step every print_step ****************************
  // print solution in xdmf format
  if((t_step+1-t_in)%print_step == 0) {
    _mgeqmap.print_soln(t_step+1);      // print sol.N.h5 and sol.N.xmf
    _mgeqmap._mgmesh.write_c(t_step+1); // print mesh.N.h5
  }
  time += dt;
#if PRINT_TIME==1 // only for cpu time check -----------------------------------
  std::clock_t    end_time2=std::clock();
  std::cout <<" Time solver ----->= " << double(end_time- start_time) / CLOCKS_PER_SEC
            <<" Time printing ----->= " << double(end_time2- end_time) / CLOCKS_PER_SEC <<
            std::endl;
#endif  // ----------------------------------------------------------------------

  return;
}



// ============================================================================
/// Xdmf transient  print
 /// This function  prints a file time.****.xmf in RESU
  void  MGTimeLoop::transient_print_xmf(
    int t_init                          ///< initial time iteration     (in)
)  {// ========================================================================

  // time parameters
  const int ndigits     = _mgutils.get_par("ndigits");
  const int print_step   = (int)(_mgutils.get_par("printstep") +0.5);
  const int nsteps       = (int)(_mgutils.get_par("nsteps") +0.5);
// 	for(int imesh=0;imesh<NUM_MESH_MAIN;imesh++){
  // dir names
  std::string basetime = _mgutils.get_file("BASETIME");
  std::string  basesol = _mgutils.get_file("BASESOL");
  std::string aux_xdmf = _mgutils.get_file("AUX_XDMF");
  // file

  std::ostringstream Name;
  Name << _mgutils._inout_dir << basetime << "."
       << setw(ndigits) << setfill('0') <<  t_init << ".xmf";
  std::ofstream out(Name.str().c_str());

  int nprt=1;
  std::string gname[3];  gname[0]=basesol;

#ifdef TWO_PHASE
  nprt=3;        gname[1]="int";        gname[2]="ccf";
#endif
  //   Mesh -----------------------------------
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM "
      <<  "\"" << _mgutils._contrib_dir << aux_xdmf << "\"" << "[]>\n";
  out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
  out << "<Domain> \n";
  for(int kp=0; kp< nprt; kp++)    {
    out << "<Grid Name=\""<< gname[kp].c_str() <<"\"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
    // time loop for grid sequence
    for(int it=t_init+1; it<=t_init+nsteps; it++) if(it%print_step ==0)   {
        out << "<xi:include href=\""
            << basesol << "."
            << setw(ndigits) << setfill('0') << it <<  ".xmf" << "\""
            << " xpointer=\"xpointer(//Xdmf/Domain/Grid["<< kp+1 <<"])\" >\n";
        out << "<xi:fallback />\n";
        out << " </xi:include>\n";
      }
    out << "</Grid> \n";
  }
  // Grid Collection end
  out << "</Domain> \n";
  out << "</Xdmf> \n";
  out.close();
//       }
//       else {
// 	        std::ostringstream Name;
//         Name << _mgutils._inout_dir << basetime << "."
//              << setw(ndigits) << setfill('0') <<  t_init << ".xmf";
//         std::ofstream out(Name.str().c_str());
//
//         int nprt=1;
//         std::string gname[3];  gname[0]=basesol;
//
// #ifdef TWO_PHASE
//         nprt=3;        gname[1]="int";        gname[2]="ccf";
// #endif
//         //   Mesh -----------------------------------
//         out << "<?xml version=\"1.0\" ?> \n";
//         out << "<!DOCTYPE Xdmf SYSTEM "
//         <<  "\"" << _mgutils._contrib_dir << aux_xdmf << "\"" << "[]>\n";
//         out << "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.2\"> \n";
//         out << "<Domain> \n";
//         for (int kp=0;kp< nprt; kp++)    {
//           out << "<Grid Name=\""<< gname[kp].c_str() <<"\"  GridType=\"Collection\" CollectionType=\"Temporal\"> \n";
//           // time loop for grid sequence
//           for (int it=t_init+1;it<=t_init+nsteps;it++) if (it%print_step ==0)   {
//               out << "<xi:include href=\""
//               << basesol<< "."
//               << setw(ndigits) << setfill('0') << it <<  ".xmf" << "\""
//               << " xpointer=\"xpointer(//Xdmf/Domain/Grid["<< kp+1 <<"])\" >\n";
//               out << "<xi:fallback />\n";
//               out << " </xi:include>\n";
//             }
//           out << "</Grid> \n";
//         }
//         // Grid Collection end
//         out << "</Domain> \n";
//         out << "</Xdmf> \n";
//         out.close();
//       }
  return;
}
