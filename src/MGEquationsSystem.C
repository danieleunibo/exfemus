// std libraries -----------------
#include <iomanip>
#include <sstream>

// class
#include "MGEquationsSystem.h"

// conf files ----------------------------
#include "Equations_conf.h"   //choose the EQNS to solve
#include "Printinfo_conf.h"
#include "MGFE_conf.h"
#include "Domain_conf.h"

// local inlcudes ----------------------
#include "MGUtils.h"
#include "MGSystem.h"
#include "MGMesh.h"
#include "MGGeomEl.h"
#include "MGFEMap.h"
#include "numeric_vectorM.h"

// classes included in the map ---------


#ifdef  NS_EQUATIONS
#include "MGSolverNS.h"  // NS ================
#endif // -----------------------------------------
#ifdef  NSA_EQUATIONS
#include "MGSolverNS-A.h"  // NS ================
#endif // -----------------------------------------
#ifdef  TBK_EQUATIONS  // Turbulence -------------
#include "MGSolverTBK.h"
#endif // -----------------------------------------
#ifdef  TBKA_EQUATIONS  // Turbulence -------------
#include "MGSolverTBK-A.h"
#endif // -----------------------------------------
#ifdef  T_EQUATIONS // Temperature ================
#include "MGSolverT.h"
#endif // -----------------------------------------
#ifdef  T_ADJ_EQUATIONS // Temperature adjoint ================
#include "MGSolverT_ADJ.h"
#endif // -----------------------------------------
#ifdef  T_G_EQUATIONS // Temperature control ================
#include "MGSolverT_G.h"
#endif // -----------------------------------------
#ifdef  T_COUP_EQUATIONS // Temperature control ================
#include "MGSolverT_COUP.h"
#endif // -----------------------------------------
#ifdef  ALFA_EQUATIONS
#include "MGSolverALFA.h"  // two fluids ===========
#endif // -----------------------------------------
#ifdef TTBK_EQUATIONS  // Turbulence -------------
#include "MGSolverTTBK.h"
#endif
#ifdef SM_EQUATIONS  // Structural Mechanics ======
#include "MGSolverSM.h"
#include "MGSolverDS.h"
#endif
#ifdef FSI_EQUATIONS  // FSI =======================
#include "MGSolverFSI.h"
#include "MGSolverDS.h"
#ifdef FSIC_EQUATIONS
#include "MGSolverFSIC.h"
#endif
#endif
#ifdef COLOR_EQUATIONS
#include "MGSolverCOL.h"
#endif

#include "MGSolverDA.h"


// ====================================================
/// This function constructs all the MGSystems
MGEquationsSystem::MGEquationsSystem(
  MGUtils& mgutils_in,// MGUtils pointer
//   MGSystem& mgphys_in,// MGSystem pointer
  MGMesh& mgmesh_in,  // MGMesh pointer
  MGFEMap& mgfemap_in, // MGFEMap pointer
  int np_data,
  int ncell_data
):MGSystem(mgutils_in,mgmesh_in,np_data,ncell_data),
//   _mgutils(mgutils_in),
//   _mgphys(mgphys_in),
//   _mgmesh(mgmesh_in),
  _mgfemap(mgfemap_in)   {// ====================================
}


// ====================================================
MGEquationsSystem::~MGEquationsSystem() {
  clean(); //deallocates the map of equations

}

// ====================================================
/// This function destroys all the MGSystems
void MGEquationsSystem::clean() {
  /*
   for(MGEquationsSystem::iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
     delete eqn->second;
   } */
  _equations.clear();
}

// =====================================================
// Table external fields
// ----------------------------------------------------


// ============================================
/// This function set the various MGSystems
void MGEquationsSystem::init(const std::vector<FIELDS> & pbName)  {

  int nvars_in[3]; nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=0;
  int n_equations=pbName.size();
//  for(int iname=0;iname<n_equations;iname++){


#ifdef DA_EQUATIONS // 
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== DA_F) {
      //number of variables of [0] piecewise [1] linear and [2] quadratic order
      nvars_in[0]=1;  nvars_in[1]=3; nvars_in[2]=1;
      MGSolDA   *mgsDA=new MGSolDA(*this,nvars_in);  // class def
      set_eqs(mgsDA);                                // set class -> equation_map
    }
#endif   // end DA_EQUATION  (Test system)


// ====================================================================
#ifdef NS_EQUATIONS // ---------  Navier-Stokes -----------
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== NS_F || pbName[iname]==NSX_F ||
        pbName[iname]==NSY_F || pbName[iname]==NSZ_F) {

      nvars_in[1]=NS_EQUATIONS%2;
      nvars_in[2]=((NS_EQUATIONS==2)?1:DIMENSION);// quadratic(2) Approximation
// if (_mgutils.get_file("MESHNUMBER")=="mesh1") {
#if NS_EQUATIONS==2     // - NS_EQUATIONS==2 -
      MGSolNS* mgs= new MGSolNS(*this,nvars_in,"NS0X","u");  set_eqs(mgs);
      MGSolNS* mgs2=new MGSolNS(*this,nvars_in,"NS0Y","v");  set_eqs(mgs2);
#if DIMENSION==3
      MGSolNS* mgs3=new MGSolNS(*this,nvars_in,"NS0Z","w");  set_eqs(mgs3);
#endif
#endif
#if (NS_EQUATIONS==0 || NS_EQUATIONS==1)      // - NS_EQUATIONS==0,1 -
      MGSolNS* mgs=new MGSolNS(*this,nvars_in);  set_eqs(mgs);
#endif

#if (NS_EQUATIONS%2==0) // - NS_EQUATIONS==0,2  projection -
      nvars_in[0]=0;  nvars_in[1]=1; nvars_in[2]=0;// only  Linear(1) approx
      MGSolNSP   *mgsP=new MGSolNSP(*this,nvars_in);  set_eqs(mgsP);
#endif
    }
#endif // -------------  end  Navier-Stokes ---------------

// ====================================================================
#ifdef NSA_EQUATIONS // ---------  Navier-Stokes Adjoint -----------
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== FS_F) {

      nvars_in[1]=1; nvars_in[2]=((NS_EQUATIONS==2)?1:DIMENSION); // quadratic(2) Approximation
// if (_mgutils.get_file("MESHNUMBER")=="mesh1") {
#if NSA_EQUATIONS==2     // - NS_EQUATIONS==2 -
      MGSolNSA* mgs= new MGSolNSA(*this,nvars_in,"NSA0X","ua");   set_eqs(mgs);
      MGSolNSA* mgs2=new MGSolNSA(*this,nvars_in,"NSA0Y","va");   set_eqs(mgs2);
#if DIMENSION==3
      MGSolNSA* mgs3=new MGSolNSA(*this,nvars_in,"NSA0Z","wa");   set_eqs(mgs3);
#endif
#endif
#if (NSA_EQUATIONS==0 || NSA_EQUATIONS==1)      // - NS_EQUATIONS==0,1 -
      MGSolNSA* mgs=new MGSolNSA(*this,nvars_in);  set_eqs(mgs);
#endif

    }
#endif // -------------  end  Navier-Stokes Adjoint ---------------



#ifdef ALFA_EQUATIONS
  nvars_in[0]=0; nvars_in[1]=1; nvars_in[2]=0;
  MGSolALFA   *mgsALFA=new MGSolALFA(*this,nvars_in, "ALFA", "a");     // class def
  set_eqs(mgsALFA);                                 // set class -> equation_map
#endif

// ====================================================================
// Energy-Equation
// ====================================================================
#ifdef T_EQUATIONS
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== T_F) {
      nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=1;  // only Quadratic[2] approx
      MGSolT   *mgsT=new MGSolT(*this,nvars_in); set_eqs(mgsT);     // class def
    }
#endif

// ====================================================================
// Energy-Equation Boundary control
// ====================================================================

#ifdef T_ADJ_EQUATIONS
  for(int iname=0; iname<n_equations; iname++) if(pbName[iname]== TA_F) {

      nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=1;
      MGSolT_ADJ   *mgsT_adj=new MGSolT_ADJ(*this,nvars_in); set_eqs(mgsT_adj);     // class def
    }
#endif
#ifdef T_G_EQUATIONS
  for(int iname=0; iname<n_equations; iname++) if(pbName[iname]== FS_F) {
      nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=1;
      MGSolT_G   *mgsT_g=new MGSolT_G(*this,nvars_in); set_eqs(mgsT_g);     // class def
    }
#endif

#ifdef T_COUP_EQUATIONS
  for(int iname=0; iname<n_equations; iname++) if(pbName[iname]== TA_F) {
#if T_COUP_EQUATIONS == 1
      nvars_in[0]=0;  nvars_in[1]=0; nvars_in[2]=2;
#else
      nvars_in[0]=0;  nvars_in[1]=1; nvars_in[2]=1;
#endif
      MGSolT_COUP   *mgsT_coup=new MGSolT_COUP(*this,nvars_in); set_eqs(mgsT_coup);     // class def
    }
#endif

// ====================================================================
#ifdef SM_EQUATIONS // ---------  STRUCTURAL MECHANICS -----------
// ====================================================================
  for(int iname=0; iname<n_equations; iname++)   if(pbName[iname]== SM_F) {
      nvars_in[0]=0; nvars_in[1]=SM_EQUATIONS%2;   // Costant(1)  Linear(0)
      nvars_in[2]=((SM_EQUATIONS==2)?1:DIMENSION); // quadratic(2) Approximation


#if (SM_EQUATIONS==2)     // - SM_EQUATIONS==2 --------------------------
      MGSolSM* mgsm= new MGSolSM(*this,nvars_in,"SM0X","u");   set_eqs(mgms);
      MGSolSM* mgsm2=new MGSolSM(*this,nvars_in,"SM0Y","v");   set_eqs(mgsm2);
#if (DIMENSION==3)
      MGSolSM* mgsm3=new MGSolSM(*this,nvars_in,"SM0Z","w"); set_eqs(mgsm3);
#endif
#else                   // - SM_EQUATIONS==0,1 -----------------------
      MGSolSM* mgsm=new MGSolSM(*this,nvars_in);  set_eqs(mgsm);
#endif
    }
#endif // =================  end SM ===================================
// ====================================================================


// ====================================================================
#ifdef FSI_EQUATIONS //    FLUID-STRUCTURE ---  FSI
// ====================================================================
  for(int iname=0; iname<n_equations; iname++) if(pbName[iname]== FS_F) {

      nvars_in[0]=0; nvars_in[1]=FSI_EQUATIONS%2;   // Costant(1)  Linear(0)
      nvars_in[2]=((FSI_EQUATIONS==2)?1:DIMENSION); // quadratic(2) Approximation
#if FSI_EQUATIONS!=2     // - NS_EQUATIONS==2 -
      MGSolFSI* mgs=new MGSolFSI(*this,nvars_in,"FSI0"); set_eqs(mgs);
#else
      MGSolFSI* mgs= new MGSolFSI(*this,nvars_in,"FSI0X","u");  set_eqs(mgs);
      MGSolFSI* mgs2=new MGSolFSI(*this,nvars_in,"FSI0Y","v");  set_eqs(mgs2);
#if DIMENSION==3
      MGSolFSI* mgs3=new MGSolFSI(*this,nvars_in,"FSI0Z","w");  set_eqs(mgs3);
#endif
#endif
    }
#ifdef DS_EQUATIONS // -------------  Displacement --------------------
  nvars_in[0]=0; nvars_in[1]=0;   // Costant(1)  Linear(0)
  nvars_in[2]=1;                  // quadratic(2) Approximation

  MGSolDS* mgsdsdx=new MGSolDS(*this,nvars_in,"SDSX","dx");  set_eqs(mgsdsdx);
  MGSolDS* mgsdsdy=new MGSolDS(*this,nvars_in,"SDSY","dy");  set_eqs(mgsdsdy);
#if (DIMENSION==3)
  MGSolDS* mgsdsdz=new MGSolDS(*this,nvars_in,"SDSZ","dz");  set_eqs(mgsdsdz);
#endif
#if (FSI_EQUATIONS%2==0) // - NS_EQUATIONS==0,2  projection -
  nvars_in[0]=0;  nvars_in[1]=1; nvars_in[2]=0;// only  Linear(1) approx
  MGSolFSIP   *mgsP=new MGSolFSIP(*this,nvars_in);  set_eqs(mgsP);
#endif
#endif // ----------------  end  Disp ---------------------------------  

#ifdef COLOR_EQUATIONS //   COLORS
// ====================================================================
//    if (_mgutils.get_file("MESHNUMBER")=="mesh1") {
  nvars_in[0]=1; nvars_in[1]=0;   // Costant(1)  Linear(0)
  nvars_in[2]=1; // quadratic(2) Approximation

  MGSolCOL* mgscol= new MGSolCOL(*this,nvars_in);   set_eqs(mgscol);
#endif // =================  end COLORS ==================================
// ====================================================================

#endif // =================  end FSI ==================================



// ========================================================================================
// ==================  Turbulence NS   ====================================================
#ifdef TBK_EQUATIONS
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== K_F ||  pbName[iname]==EW_F) {
      nvars_in[0]=0;                     // Costant(0)
      nvars_in[1]=0;                     // Linear(1)
      nvars_in[2]=(TBK_EQUATIONS%2)+1;   // Quadratic(2)
// ---------------  Turbulence K model -> 1 ----------------------------------------------
#if ((TBK_EQUATIONS/2)==0)
      MGSolTBK   *mgsTBK=new MGSolTBK(*this,nvars_in); // class def
      set_eqs(mgsTBK);                                 // set class -> equation_map
#endif    // end K-e model ++++++++++++++++++++++

#if ((TBK_EQUATIONS/2)==1)    // Turbulence K-E model -> 1 +++++++++++++++++++++++++++++++
#if (TBK_EQUATIONS%2==0)                        // splitting [2]
      MGSolTBK   *mgsTBK=new MGSolTBK(*this,nvars_in,"K","kt");  // class def (kt)
      set_eqs(mgsTBK);                                           // set class -> equation_map
      MGSolTBK   *mgsTBW=new MGSolTBK(*this,nvars_in,"K2","et"); // class def (et)
      set_eqs(mgsTBW);                                           // set class -> equation_map
#else                                           // coupled [3]  (kt,et)
      MGSolTBK   *mgsTBK=new MGSolTBK(*this,nvars_in);           // class def  (kt,et)
      set_eqs(mgsTBK);                                           // set class -> equation_map
#endif
#endif               // end K-e model ++++++++++++++++++++++

#if ((TBK_EQUATIONS/2)==2)    // Turbulence KW model -> 1 +++++++++++++++++++++++++++++++ 
#if (TBK_EQUATIONS%2==0)                         // splitting  [4]
      MGSolTBK   *mgsTBK=new MGSolTBK(*this,nvars_in,"K","kt");  // class def
      set_eqs(mgsTBK);                                           // set class -> equation_map
      MGSolTBK   *mgsTBW=new MGSolTBK(*this,nvars_in,"K2","wt"); // class def
      set_eqs(mgsTBW);                                           // set class -> equation_map
#else                                            // coupled [5]
      MGSolTBK   *mgsTBK=new MGSolTBK(*this,nvars_in);         // class def   (kt,et)
      set_eqs(mgsTBK);                                         // set class -> equation_map
#endif                                                     // end splitting-coupled 
#endif     // end K-w model ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
#endif    // end TBK model =============================================================

// ========================================================================================
// ==================  Turbulence adjoint NS   ====================================================
#ifdef TBKA_EQUATIONS
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== KTT_F) {
      nvars_in[0]=0;                     // Costant(0)
      nvars_in[1]=0;                     // Linear(1)
      nvars_in[2]=(TBKA_EQUATIONS%2)+1;   // Quadratic(2)
// ---------------  Turbulence K model -> 1 ----------------------------------------------
#if ((TBKA_EQUATIONS/2)==0)
      MGSolTBKA   *mgsTBKA=new MGSolTBKA(*this,nvars_in); // class def
      set_eqs(mgsTBKA);                                 // set class -> equation_map
#endif    // end K-e model ++++++++++++++++++++++
#if ((TBKA_EQUATIONS/2)==2)    // Turbulence KW model -> 1// coupled [5] +++++++++++++++++++++++++++++++ 
      MGSolTBKA   *mgsTBKA=new MGSolTBKA(*this,nvars_in);         // class def   (kt,et)
      set_eqs(mgsTBKA);                                         // set class -> equation_map
#endif     // end K-w model ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
#endif    // end TBK adjoint model =============================================================



// ============================================================================================
// ================     Turbulence energy  ====================================================
#ifdef TTBK_EQUATIONS
  for(int iname=0; iname<n_equations; iname++)
    if(pbName[iname]== KTT_F || pbName[iname]== EWTT_F) {
      nvars_in[0]=0;                     // Costant   -> 0
      nvars_in[1]=0;                     // Linear    -> 1
      nvars_in[2]=(TTBK_EQUATIONS%2)+1;  // Quadratic -> 2
// --------------------- Turbulence K equation (kh)  ------------------------------------------
#if ((TTBK_EQUATIONS/2)==0)    // Turbulence K model -> 1 +++++++++++++++++++++++++++++++++++++
      MGSolTTBK   *mgsTTBK=new MGSolTTBK(*this,nvars_in);          // class def   (kt)
      set_eqs(mgsTTBK);                                            // set class -> equation_map
#endif    // end K-e model ++++++++++++++++++++++
// --------------------- Turbulence K-e equation (kh,eh)---------------------------------------
#if ((TTBK_EQUATIONS/2)==1)          // Turbulence K-E model -> 1 +++++++++++++++++++++++++++++
#if (TTBK_EQUATIONS%2==0)                             // splitting [2]
      MGSolTTBK   *mgsTTBK=new MGSolTTBK(*this,nvars_in,"TK","kh");  // class def   (kt,et)
      set_eqs(mgsTTBK);                                              // set class -> equation_map
      MGSolTTBK   *mgsTTBW=new MGSolTTBK(*this,nvars_in,"TK2","eh"); // class def   (kt,et)
      set_eqs(mgsTTBW);                                              // set class -> equation_map
#else                                                  // coupled [3]
      MGSolTTBK   *mgsTTBK=new MGSolTTBK(*this,nvars_in);            // class def   (kt,et)
      set_eqs(mgsTTBK);                                              // set class -> equation_map
#endif
#endif                             // end K-e model +++++++++++++++++++++++++++++++++++++++++++
// -----------------------Turbulence K-w equation (kh,wh)--------------------------------------
#if ((TTBK_EQUATIONS/2)==2)    // Turbulence KW model -> 1 ++++++++++++++++++++++++++++++++++++
#if (TTBK_EQUATIONS%2==0)          // splitting [4]
      MGSolTTBK   *mgsTTBK=new MGSolTTBK(*this,nvars_in,"TK","kh");  // class def   (kh)
      set_eqs(mgsTTBK);                                              // set class -> equation_map
      MGSolTTBK   *mgsTTBW=new MGSolTTBK(*this,nvars_in,"TK2","wh"); // class def   (wh)
      set_eqs(mgsTTBW);// set class -> equation_map
#else                                                 // coupled [5]
      MGSolTTBK   *mgsTTBK=new MGSolTTBK(*this,nvars_in);             // class def   (kh,wh)
      set_eqs(mgsTTBK);                                               // set class -> equation_map
#endif                   // end splitting-coupled 
#endif     // end K-w model ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    }
#endif  // ========================end Turbulence energy======================================


  // ===================================  external ==============================================
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
    MGSolBase* mgsol = eqn->second;// get the pointer
    mgsol -> set_ext_fields(pbName);          // init ext fields
  }
  // ================================================================================
  return;
}

// ====================================================
/// This sets dof initial and boundary conditions and sets the operators
void  MGEquationsSystem::setDofBcOpIc() {

  // Reading operators
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
    MGSolBase* mgsol = eqn->second;// get the pointer
    mgsol -> MGDofBcOp();          // init dof, GenBc, ReadOperators
    mgsol -> GenIc();              // initial solution
  }
  return;
}

// ==========================================================================================
/// This function performes all the MGSystem time step routines
void MGEquationsSystem::eqnmap_timestep_loop(
  const double time,             // real time
  const int delta_t_step_in     // integer time
) {
  // loop for time steps
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
    MGSolBase* mgsol = eqn->second;
    mgsol -> MGTimeStep(time,delta_t_step_in);
  }
  return;
}
// ==========================================================================================
/// This function performes all the MGSystem time step routines for control problem
void MGEquationsSystem::eqnmap_timestep_loop_control(
  const double time,             // real time
  const int delta_t_step_in     // integer time
) {

  MGSolBase* mgsolns  =  get_eqs("NS0");
#if NS_EQUATIONS%2==0
  MGSolBase* mgsolnsp =  get_eqs("NSP");
#endif
  MGSolBase* mgsolk   =  get_eqs("K");
  MGSolBase* mgsolnsa =  get_eqs("NSA0");
  MGSolBase* mgsolka  =  get_eqs("KA");
  
  const int NoLevels=mgsolns->_NoLevels;
  
  mgsolnsa->x_old[NoLevels-1]->localize(*mgsolnsa ->x_oold[NoLevels-1]);
  
  const int max_iter=3000; const int max_iter_adj=8000; const int max_func_iter=50; const double toll_conv=1.e-6;
//   mgsolns ->_dumping = 2000.;
  
//   std::cout << "\nDumping now: " << mgsolns ->_dumping << endl; 
  for(int isolns=0;isolns<max_iter;isolns++)  {
    (mgsolns->x_old[NoLevels-1])->localize(*mgsolns ->disp[NoLevels-1]);
    (mgsolk->x_old[NoLevels-1])->localize(*mgsolk ->disp[NoLevels-1]);
    
//     std::cout << " Non linear iteration: ";
//     for(int inonlinear=0;inonlinear<5;inonlinear++)  {
      mgsolns -> MGTimeStep(time,delta_t_step_in);
#if NS_EQUATIONS%2==0
      mgsolnsp -> MGTimeStep(time,delta_t_step_in);
#endif
      mgsolk -> MGTimeStep(time,delta_t_step_in);
//       std::cout<< " " << isolns+1 << endl;
//     }
    const double err_norm_ns=fabs(((mgsolns->disp[NoLevels-1]->l2_norm())-(mgsolns ->x_old[NoLevels-1]->l2_norm()))/
	(mgsolns ->x_old[NoLevels-1]->l2_norm()));
    const double err_norm_k=fabs(((mgsolk->disp[NoLevels-1]->l2_norm())-(mgsolk ->x_old[NoLevels-1]->l2_norm()))/
	(mgsolk ->x_old[NoLevels-1]->l2_norm()));
    if (err_norm_ns< toll_conv && err_norm_k< toll_conv) {
      std::cout << "\n Steady state NS-K found, iteration " << isolns+1 << '\n';
      break;
    }
    if ( isolns == max_iter-1) 
      std::cout << "\n Maximum iteration reached in NS!!! Error norm is " << err_norm_ns<< " and " << err_norm_k << '\n';
  }
//   mgsolns ->_dumping = mgsolns ->_dumping/10.;
//   
//   
//   
  double func0= mgsolns->MGFunctional(0,1.); //compute functional with _eta=1
//   std::cout << "\n \n Divergence is " << func0 << endl; 
//   
  for(int isolad=0;isolad<max_iter_adj;isolad++)  {
    (mgsolnsa->x_old[NoLevels-1])->localize(*mgsolnsa ->disp[NoLevels-1]);
    (mgsolka->x_old[NoLevels-1])->localize(*mgsolka ->disp[NoLevels-1]);
    mgsolnsa -> MGTimeStep(time,delta_t_step_in);
    mgsolka -> MGTimeStep(time,delta_t_step_in);
    const double err_norm_nsa=fabs(((mgsolnsa ->disp[NoLevels-1]->l2_norm())-(mgsolnsa ->x_old[NoLevels-1]->l2_norm()))/
      (mgsolnsa ->x_old[NoLevels-1]->l2_norm()));
    const double err_norm_ka=fabs(((mgsolka ->disp[NoLevels-1]->l2_norm())-(mgsolka ->x_old[NoLevels-1]->l2_norm()))/
      (mgsolka ->x_old[NoLevels-1]->l2_norm()));
    if (err_norm_nsa < toll_conv && err_norm_ka < toll_conv) {
      std::cout << "\n Steady state NSA-KA found , iteration " << isolad+1 << '\n';
      break;
    }
    if ( isolad == max_iter_adj-1) 
      std::cout << "\n Maximum iteration reached in NSA!!! Error norm is " 
      << err_norm_nsa << " " << err_norm_ka << '\n';
  }

  (mgsolns->x_old[NoLevels-1])->localize(*mgsolns ->x_user[NoLevels-1]); //store the good value of NS
	 
  for(int isolfunc=0;isolfunc<max_func_iter;isolfunc++)  {
    for(int isolns=0;isolns<max_iter;isolns++)  {
      (mgsolns->x_old[NoLevels-1])->localize(*mgsolns->disp[NoLevels-1]);
      mgsolns -> MGTimeStep(time,delta_t_step_in);
#if NS_EQUATIONS%2==0
      mgsolnsp -> MGTimeStep(time,delta_t_step_in);
#endif
      mgsolk -> MGTimeStep(time,delta_t_step_in);
    const double err_norm_ns=fabs(((mgsolns->disp[NoLevels-1]->l2_norm())-(mgsolns ->x_old[NoLevels-1]->l2_norm()))/
	(mgsolns ->x_old[NoLevels-1]->l2_norm()));
    const double err_norm_k=fabs(((mgsolk->disp[NoLevels-1]->l2_norm())-(mgsolk ->x_old[NoLevels-1]->l2_norm()))/
	(mgsolk ->x_old[NoLevels-1]->l2_norm()));
    if (err_norm_ns< toll_conv && err_norm_k< toll_conv) {
      std::cout << "\n Steady state NS-K found, iteration " << isolns+1 << '\n';
      break;
    }
      if ( isolns == max_iter-1) 
	std::cout << "\n Maximum iteration reached in NS!!! Error norm is " << err_norm_ns << " and " << err_norm_k << '\n';
    }
    double func1= mgsolns->MGFunctional(0,0);
    std::cout <<"\nIter : " << isolfunc+1<< " old functional is " << func0 << " and new " << func1; 
    if ((fabs(func0-func1)/func0) < toll_conv) {
      std::cout << "\n Convergence of the optimal control problem reached for equal functionals!!";
      mgsolns->MGFunctional(3,0);  // call compute adjoint to save the good control
      (mgsolns->x_user[NoLevels-1])->localize(*mgsolnsa ->x_old[NoLevels-1]);
      mgsolns->MGFunctional(1,0.);   //set _eta=0;
      (mgsolns->x_old[NoLevels-1])->localize(*mgsolns ->x_user[NoLevels-1]);
      break; 
    }
    else if (func0<func1) {
      mgsolns->MGFunctional(1,0.66);  // set _eta=0.5*_eta
      (mgsolns->x_user[NoLevels-1])->localize(*mgsolns ->x_old[NoLevels-1]); //reset the good value of NS
    }
    else if (func0>func1) {
      mgsolns->MGFunctional(3,0);  // call compute adjoint to save the good control
      (mgsolns->x_user[NoLevels-1])->localize(*mgsolnsa ->x_old[NoLevels-1]);
      mgsolns->MGFunctional(1,1.5);   //set _eta=2.*_eta;
      (mgsolns->x_old[NoLevels-1])->localize(*mgsolns ->x_user[NoLevels-1]);
      break;
    }
    if(isolfunc==max_func_iter-1) {
      std::cout << "\n Convergence of the optimal control problem reached for maximum iteration!! \n";
// // //       mgsolns->MGFunctional(1,0.);   //set _eta=0;
      (mgsolns->x_user[NoLevels-1])->localize(*mgsolns ->x_old[NoLevels-1]); //reset the good value of NS
      mgsolns->MGFunctional(3,0);  // call compute adjoint to save the good control
      (mgsolns->x_user[NoLevels-1])->localize(*mgsolnsa ->x_old[NoLevels-1]);
      break;
    }
  }

  std::cout << endl;

 
  return;
}

// =================================================================


// ==========================================================================================
/// This function performes all the MGSystem time step routines
void MGEquationsSystem::eqnmap_timestep_loop_nonlinear(
  const double time,             // real time
  const int delta_t_step_in     // integer time
) {
  // loop setup
  std::cout << " ============= NON LINEAR ITERATOR :setup ===========  "<<std::endl;
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
    MGSolBase* mgsol = eqn->second;  mgsol -> MGTimeStep_nl_setup(time,delta_t_step_in);
  }
  //  nonlinear loop
  /// D) Update of the old solution at the top Level  (MGSolFSI::OldSol_update),

//     disp[_NoLevels-1]->zero();
//   x_oold[_NoLevels-1]->close();
//   x_oold[_NoLevels-1]->zero();
  int iter=0; int itest_count=0;int itest=0;
  while(iter<3 &&  itest_count==0 ) { 
    std::cout << " ============= NON LINEAR ITERATOR============  "<< iter<<std::endl;
    itest_count=0;itest=0;
    for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {  
      MGSolBase* mgsol = eqn->second;
      itest=mgsol -> MGTimeStep_nl_iter(time,delta_t_step_in);
      itest_count+=itest;
    }
    iter++;
  }
  std::cout << " ============= NON LINEAR ITERATOR :solution up ===========  "<<std::endl;
  for(iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++)  {
    MGSolBase* mgsol = eqn->second;
    mgsol -> MGTimeStep_nl_sol_up(time,delta_t_step_in);
  }

return;
}
// =================================================================


/// This function prints xdmf and hdf5 file
void MGEquationsSystem::print_soln(const int t_step // time step
                                  ) {

  const int    iproc =_mgmesh._iproc;
  if(iproc==0) {  // print only one processor

    print_soln_h5(t_step);                   // print sol h5
    int n_lines=0,n_cells=0;
// #ifdef TWO_PHASE  // --------- cc ----------------
//     print_h5CC(file,t_flag,n_lines,n_lines); // print CC
// #endif // ----------------- cc --------------------
    print_soln_xmf(t_step,n_lines,n_cells);  // print xdmf file
  }

  return;
}

// =================================================================
/// This function prints the attributes into the corresponding hdf5 file
void MGEquationsSystem::print_soln_h5(const int t_flag // time flag
                                     ) {

  const int NoLevels = _mgutils.get_par("nolevels");
  const int ndigits  = _mgutils.get_par("ndigits");

  // file  ---------------------------------------------
  // file name
  std::ostringstream filename;
  filename << _mgutils._inout_dir << _mgutils.get_file("BASESOL") << "." << setw(ndigits) << setfill('0') << t_flag << ".h5";
  // open file for hf5 storage
  hid_t   file= H5Fcreate(filename.str().c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  H5Fclose(file);

  // print all systems ---------------------------
  MGEquationsSystem::const_iterator pos=_equations.begin();
  MGEquationsSystem::const_iterator pos_e=_equations.end();
  for(; pos!=pos_e; pos++)    {
    MGSolBase *mgsol=pos->second;
    mgsol->print_u(filename.str(),NoLevels-1);
//     mgsol->print_u_med("mesh_sol.med",NoLevels-1);
  }
// printing cell system data (MGSystem-> printdata)
  for(int idata=0; idata<_n_data[0]+_n_data[1]; idata++) {
    std::ostringstream dir_name;
    dir_name << "DATA" << idata;
    print_data_view(filename.str(),idata,dir_name.str());
  }
  return;
}

// ===================================================================
/// It prints the attributes in  Xdmf format for one time step
void MGEquationsSystem::print_soln_xmf(const int t_step, int /*n_lines*/,int /*n_cells*/) {

  const int NoLevels = _mgutils.get_par("nolevels");
  const int ndigits  = _mgutils.get_par("ndigits");

  //  Mesh ----------------------
  const MGMesh& mgmesh=_mgmesh;
  int n_nodes    = mgmesh._NoNodes[NoLevels-1];
  int n_elements    = mgmesh._NoElements[0][NoLevels-1];

  // get parameters
  std::string inout_dir  = _mgutils._inout_dir;
  std::string basesol    = _mgutils.get_file("BASESOL");
  std::string basemesh   = _mgutils.get_file("BASEMESH");
  std::string contrib_dir= _mgutils.get_file("CONTRIB_DIR");
//   std::string aux_xdmf   = _mgutils.get_file("AUX_XDMF");
  std::string connlin    = _mgutils.get_file("CONNLIN");
  const double dt = _mgutils.get_par("dt");

  // files
//   std::ostringstream conn_file; // connectivity file (mesh_conn_lin.h5)
//   conn_file /* <<femus_dir  <<  "/" << appl_dir << "/"  << myapp << "/" << input_dir */ << basemesh;
  std::ostringstream topol_file; // topology file file (mesh_conn_lin.h5)
  topol_file  << basemesh   << connlin <<".h5";
//   conn_file << ".h5";
  std::ostringstream coord_time_file; // connectivity file (mesh_conn_lin.h5)
  coord_time_file /*<<femus_dir  << "/"<< output_dir << "/"*/ << basemesh << "." << std::setw(ndigits) << std::setfill('0') << t_step << ".h5";

  std::ostringstream filename; //  solution file xmf
  filename << inout_dir  << basesol << "." << setw(ndigits) << setfill('0') << t_step;
  std::ostringstream attr_file;//  solution file h5
  /* attr_file << filename.str()<< ".h5"; */
  attr_file << basesol << "." << setw(ndigits) << setfill('0') << t_step<< ".h5";
  filename << ".xmf";

  // solution file  xmf
  std::ofstream out(filename.str().c_str());
  std::cout << " " << filename.str().c_str() << std::endl;
//  ++++++++++++ Header ++++++++++++++
  out << "<?xml version=\"1.0\" ?> \n";
//   out << "<!DOCTYPE Xdmf SYSTEM ";
//   out <<  "\"" << aux_xdmf << "\" > \n";
  out << "<Xdmf> \n" << "<Domain> \n"<< "<Grid Name=\"Mesh\"> \n";
  // time
  out << "<Time Value =\"" << t_step*dt << "\" /> \n";
  // +++++ Topology ++++++++++
  _mgmesh.print_xmf_topology(out,topol_file.str(),NoLevels-1,0);
  // +++++++  Geometry +++++++++++++++++
  _mgmesh.print_xmf_geometry(out,coord_time_file/*conn_file*/.str(),NoLevels-1,0);
  // ++++  Attributes ++++++++++++
  MGEquationsSystem::const_iterator pos1   = _equations.begin();
  MGEquationsSystem::const_iterator pos1_e = _equations.end();
  for(; pos1!=pos1_e; pos1++)   {
    MGSolBase *mgsol=pos1->second;
    mgsol->print_xml_attrib(out,n_nodes,n_elements*4,attr_file.str());
  }
//   printining cell attributes
  if(_n_data[0]+_n_data[1]>0) print_xml_attrib(out,n_elements,n_nodes,attr_file.str());

// #ifdef TWO_PHASE
//   // print of CC
//   print_xmfCC(out,t_step,n_lines,n_cells);
// #endif
  out << "</Grid>\n" << "</Domain> \n" << "</Xdmf> \n";
  out.close();
  return;
}

// ========================================================================
/// This function read the solution for all the system (restart)
void MGEquationsSystem::read_soln(const int t_step) {
// --------------------------------------------------------------

  const int ndigits  = _mgutils.get_par("ndigits");
  const  int restart_lev_flag = _mgutils.get_par("restart_lev");

  // open file -----------------------------
  std::ostringstream namefile;
  namefile << _mgutils._inout_dir << _mgutils.get_file("BASESOL") << "."
           << setw(ndigits) << setfill('0') << t_step << ".xmf";

#ifdef PRINT_INFO // --------  info ------------------ 
  std::cout << "\n MGEquationsSystem::read_soln: Reading time  from "
            << namefile.str().c_str();
#endif  // -------------------------------------------
  std::ifstream in ; in.open(namefile.str().c_str());  //associate the file stream with the name of the file
  if(!in.is_open()) { std::cout << " MGCase: restart .xmf file not found "  << std::endl; abort(); }

  // reading time from xmf file --------------
  std::string buf="";  while(buf != "<Time") in >> buf;
  in >> buf >> buf;  buf=buf.substr(2,buf.size()-3);
//create an istringstream from a string
  std::istringstream buffer(buf); double restart_time;  buffer >> restart_time;

  //add parameter to system
  _mgutils.set_par("restartime",restart_time);

  // ---------------------------------------------------
  // reading data from  sol.N.h5
  // ---------------------------------------------------
  // file name -----------------------------------------
  namefile.str("");  //empty string
  namefile << _mgutils._inout_dir << _mgutils.get_file("BASESOL")
           << "." << setw(ndigits) << setfill('0') << t_step << ".h5";

#ifdef PRINT_INFO  // --------------- info ---------------
  std::cout << "\n MGEquationsSystem::read_soln: Reading from file "
            << namefile.str().c_str() << std::endl;
#endif // ---------------------------------------------
  // loop reading over the variables ---------------------
  for(MGEquationsSystem::const_iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
    MGSolBase *mgsol=eqn->second;
    mgsol->read_u(namefile.str(),restart_lev_flag);
  } //  loop --------------------------------------------------------

// #ifdef TWO_PHASE
//   readCC(t_step);
// //   _mgsys.get_mgcc()->read_fine ( t_init,_mgsys.get_mgcc()->_Nlev-1 );
// //   readCC ( t_step );
// #endif

  return;
}


// ========================================================================
/// This function print data from single class equation to the mesh system  on vect_data
void MGEquationsSystem::print_mesh_data(double vect_data[]) {
// --------------------------------------------------------------

  // loop reading/printing over the equation ---------------------
  int count=0;
//   eqn=_equations.begin();
  for(MGEquationsSystem::const_iterator eqn=_equations.begin(); eqn != _equations.end(); eqn++) {
    if(_mgutils.get_file("MESHNUMBER")=="mesh1") {
      if(count==1) {
        MGSolBase *mgsol=eqn->second;
        mgsol->print_ext_data(vect_data); // compute from a system/mesh to vect_data
      }
      count++;
    }
    if(_mgutils.get_file("MESHNUMBER")=="mesh2") {
      if(count==0) {
        MGSolBase *mgsol=eqn->second;
        mgsol->print_ext_data(vect_data); // compute from a system/mesh to vect_data
      }
      count++;
    }
  } //  loop --------------------------------------------------------

  return;
}
// ====================================================

// ====================================================
/// This function prints initial and boundary data in xdmf+hdf5 format
void MGEquationsSystem::print_case(const int t_init) {

  const int    iproc =_mgmesh._iproc;

  if(iproc==0) {  // print only one processor
    print_case_h5(t_init); // ic+bc print format h5
    int n_lines=0,n_cells=0; // for VOF
// #ifdef TWO_PHASE
//     print_h5CC(hid_t file,t_init,&n_lines,&n_cells);
// #endif
    print_case_xmf(t_init,n_lines,n_cells); // xml format
  }
  return;
}





// =============================================================================
/// This function prints initial and boundary data in hdf5 fromat
/// in the file case.h5
void MGEquationsSystem::print_case_h5(const int t_init) {

  const int NoLevels  = _mgutils.get_par("nolevels");
  const int ndigits   = _mgutils.get_par("ndigits");
  std::string output_dir = _mgutils._inout_dir;
  std::string   basecase = _mgutils.get_file("BASECASE");
  //  Mesh ---- ---------------------------------------------
  const MGMesh& mgmesh = _mgmesh;

  // file ---------------------------------------
  std::ostringstream filename; // file name
  filename << output_dir << basecase << "." << setw(ndigits) << setfill('0') << t_init << ".h5";

  hid_t   file = H5Fcreate(filename.str().c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
  mgmesh.print_subdom_hf5(filename.str())  ; // PID (n processor)
#ifdef TBK_EQUATIONS
  mgmesh.print_dist_hf5(filename.str(),mgmesh._dist,"DIST")  ; // PID (n processor)
#endif
  H5Fclose(file);

  // loop over all systems ---------------------------
  MGEquationsSystem::const_iterator pos   = _equations.begin();
  MGEquationsSystem::const_iterator pos_e = _equations.end();
  for(; pos!=pos_e; pos++) {
    MGSolBase *mgsol=pos->second;
    mgsol->print_u(filename.str(),NoLevels-1);   // initial solution
    mgsol->print_bc(filename.str(),NoLevels-1); // boundary condition
  }

  return;
}



// ====================================================================
/// It prints the Xdmf file to read the initial and boundary conditions
void MGEquationsSystem::print_case_xmf(const int t_init,const int /*n_lines*/,
                                       const int /*n_cells*/) {

  // file ----------------------------------------------
  std::string inout_dir  = _mgutils._inout_dir;
  std::string   basecase = _mgutils.get_file("BASECASE");
  std::string   basemesh = _mgutils.get_file("BASEMESH");
  std::string  contrib_dir = _mgutils.get_file("CONTRIB_DIR");
  std::string  aux_xdmf = _mgutils.get_file("AUX_XDMF");
  std::string  connlin = _mgutils.get_file("CONNLIN");
  const int NoLevels = _mgutils.get_par("nolevels");
  const int ndigits  = _mgutils.get_par("ndigits");

  // files
  std::ostringstream conn_file; // connectivity file (mesh_conn_lin.h5)
  conn_file << basemesh;
  std::ostringstream topol_file; // topology file file (mesh_conn_lin.h5)
  topol_file << conn_file.str() << connlin <<".h5";
  conn_file << ".h5";
  std::ostringstream filename; //  solution file xmf
  filename <<  inout_dir  << basecase << "." << setw(ndigits) << setfill('0') <<  t_init<< ".xmf";

  //  Mesh ---- ---------------------------------------------
  const MGMesh&    mgmesh = _mgmesh;
  int n_nodes    = mgmesh._NoNodes[NoLevels-1];
  int n_elements = mgmesh._NoElements[0][NoLevels -1];
  std::string var_name[3];    std::string var_type[3];

  //   File xdmf -------------------------------
  std::ofstream out(filename.str().c_str());
  out << "<?xml version=\"1.0\" ?> \n";
  out << "<!DOCTYPE Xdmf SYSTEM ";
  out <<  "\"" <<  aux_xdmf << "\" \n";
//    out << " [ <!ENTITY HeavyData \"\"> ] ";
  out << ">\n";
  out << "<Xdmf> \n" << "<Domain> \n" << "<Grid Name=\"Mesh\"> \n";
  // +++++ Topology ++++++++++
  _mgmesh.print_xmf_topology(out,topol_file.str(),NoLevels-1,0);
  // +++++++  Geometry +++++++++++++++++
  _mgmesh.print_xmf_geometry(out,conn_file.str(),NoLevels-1,0);

  // ++++  Attributes for each system +++++++++++++++++++
  MGEquationsSystem::const_iterator pos1=_equations.begin();
  MGEquationsSystem::const_iterator pos1_e=_equations.end();
  for(; pos1!=pos1_e; pos1++)   {
    MGSolBase *mgsol=pos1->second;
    for(int ivar=0; ivar<mgsol->_nvars[1]+mgsol->_nvars[2]; ivar++)     {

      // Volume and boundary conditions
      var_name[0] = mgsol->_var_names[ivar];
      var_name[1] =var_name[0]+"bd";
      var_name[2] =var_name[0]+"vl";
      var_type[0] ="Float";
      var_type[1] ="Float";
      var_type[2] ="Float";
      for(int ibvar=0; ibvar<3; ibvar++) {
        out << "<Attribute Name=\""<< var_name[ibvar] <<"\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        out << "<DataItem  DataType=\""<< var_type[ibvar].c_str()
            << "\" Precision=\"8\" Dimensions=\"" << n_nodes << "  "
            << 1 << "\" Format=\"HDF\">  \n";
        out <</* femus_dir << "/" << output_dir <<*/ basecase << "."
            << setw(ndigits) << setfill('0') << t_init << ".h5"
            << ":" << var_name[ibvar].c_str() << "\n";
        out << "</DataItem>\n" << "</Attribute>\n";
      }
    }
    for(int ivar=0; ivar<mgsol->_nvars[0]; ivar++)     {
      // Volume and boundary conditions
      var_name[0] = mgsol->_var_names[mgsol->_nvars[1]+mgsol->_nvars[2]+ivar];
      var_name[1] =var_name[0]+"bd";
      var_name[2] =var_name[0]+"vl";
      var_type[0] ="Float";
      var_type[1] ="Int";
      var_type[2] ="Int";
      for(int ibvar=0; ibvar<1; ibvar++) { // only initial condition ->var_name[0]
        out << "<Attribute Name=\""<< var_name[ibvar] <<"\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
        out << "<DataItem  DataType=\""<< var_type[ibvar].c_str()
            << "\" Precision=\"8\" Dimensions=\"" << n_elements*_mgmesh._GeomEl.n_se[0] << "  "
            << 1 << "\" Format=\"HDF\">  \n";
        out << basecase << "."
            << setw(ndigits) << setfill('0') << t_init << ".h5"
            << ":" << var_name[ibvar].c_str() << "\n";
        out << "</DataItem>\n" << "</Attribute>\n";
      }
    }
  }

  // ----------------------------------------------------------
  out << "<Attribute Name=\"PID\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_elements*_mgmesh._GeomEl.n_se[0] << "  "
      << 1 << "\" Format=\"HDF\">  \n";
  out << basecase << "."
      << setw(ndigits) << setfill('0')<< t_init << ".h5" << ":PID\n";
  out << "</DataItem>\n" << "</Attribute>\n";
#ifdef TBK_EQUATIONS
  // ----------------------------------------------------------
  out << "<Attribute Name=\"DIST\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
  out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\""
      << n_elements*_mgmesh._GeomEl.n_se[0] << "  "
      << 1 << "\" Format=\"HDF\">  \n";
  out << basecase << "."
      << setw(ndigits) << setfill('0')<< t_init << ".h5" << ":DIST\n";
  out << "</DataItem>\n" << "</Attribute>\n";
#endif
// #ifdef TWO_PHASE
//   // print of CC
//   print_xmfCC(out, t_init,n_lines,n_cells);
// #endif
  out << "</Grid>\n"<< "</Domain> \n" << "</Xdmf> \n";
  out.close();

  return;
}


// #ifdef TWO_PHASE
// // ==================================================
// ///readCC
// void MGCase::readCC(const int flag_print) {
//
//   std::ostringstream filename;
//   filename <<"output/ccf." << setw(_ndig) << setfill('0') << flag_print << EXT_H5;
//   hid_t  file_id = H5Fopen(filename.str().c_str(),H5F_ACC_RDWR, H5P_DEFAULT);
//
//   hid_t dtset = H5Dopen(file_id, "/cres", H5P_DEFAULT);
//   hid_t filespace = H5Dget_space(dtset);
//   hsize_t dims[2]; H5Sget_simple_extent_dims(filespace, dims, NULL);
//
//   double *ccf = (new double[dims[0]]); double *xf = (new double[4*dims[0]]);
//   double *yf = (new double[4*dims[0]]); double *zf = (new double[4*dims[0]]);
//
//   hid_t status=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,ccf);
//   H5Dclose(dtset);
//   dtset = H5Dopen(file_id, "/X1", H5P_DEFAULT);
//   status=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,xf);
//   H5Dclose(dtset);
//   dtset = H5Dopen(file_id, "/X2", H5P_DEFAULT);
//   status=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,yf);
//   H5Dclose(dtset);
//   dtset = H5Dopen(file_id, "/X3", H5P_DEFAULT);
//   status=H5Dread(dtset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,H5P_DEFAULT,zf);
//   H5Dclose(dtset);
//
//   _mgsys.get_mgcc()->read_fine(dims[0],ccf,xf,yf,zf);
//
//   delete []ccf; delete []xf; delete []yf; delete []zf;
//   return;
// }
//
//
// /// It prints the format Xdmf for one time step of CC
// void MGCase::print_xmfCC(std::ofstream &out, const int t_step,
//                          const int n_lines, const int n_cells) {
//   //  Mesh ---- ---------------------------------------------
//   const MGMesh &mgmesh=_mgsys.get_mesh();
// //  int n_nodes=mgmesh._NoNodes[_cslevel];
//   int n_elements=mgmesh._NoElements[_cslevel];
//   const int nvrt=4* (DIMENSION-1);
// #if DIMENSION==2
//   std::string btype="Polyline";
//   std::string mtype="Quadrilateral";
// #else
//   std::string btype="Triangle";
//   std::string mtype="Hexahedron";
// #endif
// // time parameters
//   const double dt=_mgsys.get_par("dt");
//
//   // ============
//   // coarse color function
//   // ============
//   out << "<Attribute Name=\"CC\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
//   out << "<DataItem  DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_elements*nvrt << "  " << 1 <<
//       "\" Format=\"HDF\">  \n";
//   out << "sol." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":CC" << "\n";
//   out << "</DataItem>\n";
//   out << "</Attribute>\n";
//   // Grid Collection end
//   out << "</Grid> \n";
//
//   // ==================
//   //  Interface
//   // =================
//   out << "<Grid Name=\"Interf\"> \n";
//   out << "<Time Value =\"" << t_step*dt << "\" /> \n";
//   // +++++ Topology ++++++++++
//   out << "<Topology Type=\""<< btype << "\"   Dimensions=\""<<  n_lines <<    "\"> \n";
//   out << "<DataStructure DataType=\"Int\" Dimensions=\""<< n_lines <<"  "<< DIMENSION <<
//       "\" Format=\"HDF\">  \n";
//   out  << "sol." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":" << "intconn" << "\n";
//   out << "</DataStructure> \n" << "</Topology> \n";
//   // +++++++  Geometry +++++++++++++++++
//   out << "<Geometry Type=\"X_Y_Z\"> \n";
//   for(int ix=1; ix<4; ix++) {
//     out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_lines*DIMENSION << "  " << 1 <<
//         "\" Format=\"HDF\">  \n";
//     out <<  "sol." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":" << "intX" <<ix<< "\n";
//     out << "</DataStructure> \n";
//   }
//   out << " </Geometry>\n";
//   // Grid Collection end
//   out << "</Grid> \n";
//   // =========================
//   // cc on fine grid
//   // =========================
//   out << "<Grid Name=\"ccf\"> \n";
//   out << "<Time Value =\"" << t_step*dt << "\" /> \n";
//   // +++++ Topology ++++++++++
//   out << "<Topology Type=\""<< mtype << "\"   Dimensions=\""<<  n_cells <<    "\"> \n";
//   out << "<DataStructure DataType=\"Int\" Dimensions=\""<< n_cells <<"  "<< nvrt <<
//       "\" Format=\"HDF\">  \n";
//   out  << "ccf." << setw(_ndig) << setfill('0') << t_step << EXT_H5 <<":" << "conn" << "\n";
//   out << "</DataStructure> \n" << "</Topology> \n";
//   // +++++++  Geometry +++++++++++++++++
//   out << "<Geometry Type=\"X_Y_Z\"> \n";
//   for(int ix=1; ix<4; ix++) {
//     out << "<DataStructure DataType=\"Float\" Precision=\"8\" Dimensions=\"" << n_cells*nvrt << "  " << 1 <<
//         "\" Format=\"HDF\">  \n";
//     out <<  "ccf." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":X" <<ix<<"\n";
//     out << "</DataStructure> \n";
//   }
//   out << " </Geometry>\n";
//   // ccf
//   out  << "<Attribute Name=\"ccf\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
//   out  << "<DataItem Dimensions=\""<< n_cells << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\"> \n ";
//   out << "ccf." << setw(_ndig) << setfill('0') << t_step << EXT_H5 << ":" << "ccf" << "\n";
//   out << "</DataItem>\n" << "</Attribute>\n";
//
//   return;
// }
//
//
// // =========================================================
// void MGCase::print_h5CC(hid_t file,const int flag_print,
//                         int *n_l_out,int *n_c_out)  {
//   //  Mesh ---- ---------------------------------------------
//   const MGMesh &mgmesh=_mgsys.get_mesh();
// //  int n_nodes=mgmesh._NoNodes[_cslevel];
//   const int n_elements=mgmesh._NoElements[_cslevel];
//
//   int nvrt= (DIMENSION-1) *4;
//   // ==========
//   //  function cc
//   // ===========
//   double *ucc;  ucc=new double[n_elements*(DIMENSION-1) *4];
//   hsize_t dimsf[2];  dimsf[0] = n_elements* (DIMENSION-1) *4;  dimsf[1] = 1;
//   MGSolCC *mgcc=_mgsys.get_mgcc();
//   mgcc->Write_xmf(ucc);
//   std::string name="CC";  print_Dhdf5(file,name,dimsf,ucc);
//   delete []ucc;
//   // ==============
//   // interface
//   // =============
//   const int n_lines=mgcc->get_nel();
//   // line coordinates
//   double *xcoord;
//   xcoord=new double[n_lines*DIMENSION];
//   double *ycoord;
//   ycoord=new double[n_lines*DIMENSION];
//   double *zcoord;
//   zcoord=new double[n_lines*DIMENSION];
//   int *tcon;
//   tcon=new int[n_lines*DIMENSION];
//   mgcc->get_int(n_lines,tcon,xcoord,ycoord,zcoord);
//   // coordinate datasets --------------------
//   dimsf[0] = n_lines*DIMENSION;  dimsf[1] = 1;
//   std::string name="intX1";
//   print_Dhdf5(file,name,dimsf,xcoord);
// //   dataspace = H5Screate_simple ( 2, dimsf, NULL );
// //   dataset = H5Dcreate ( file, "intX1",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, xcoord );
// //   H5Dclose ( dataset );
//   name="intX2";
//   print_Dhdf5(file,name,dimsf,ycoord);
// //   dataset = H5Dcreate ( file, "intX2",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, ycoord );
// //   H5Dclose ( dataset );
//   name="intX3";
//   print_Dhdf5(file,name,dimsf,zcoord);
// //   dataset = H5Dcreate ( file, "intX3",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, zcoord );
// //   H5Dclose ( dataset );
// //   H5Sclose ( dataspace );
//   // line connectivity ------------------------
//   dimsf[0] = n_lines;  dimsf[1] = DIMENSION;
//   name="intconn";    print_Ihdf5(file,name,dimsf,tcon);
// //   dataspace = H5Screate_simple ( 2, dimsf, NULL );
// //   dataset = H5Dcreate ( file, "intconn", H5T_NATIVE_INT, dataspace,
// //                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT,tcon );
// //   H5Dclose ( dataset );
//   // Close/release resources.
// //   H5Sclose ( dataspace );
//   delete []tcon;  delete []xcoord;
//   delete []ycoord;  delete []zcoord;
//
//   // ==========
//   //  ccf
//   // ==========
//   std::ostringstream filename;
//   filename << "output/ccf." << setw(_ndig) << setfill('0') << flag_print << EXT_H5;
//   hid_t file1= H5Fcreate(filename.str().c_str(),H5F_ACC_TRUNC, H5P_DEFAULT,H5P_DEFAULT);
//
//   // ccf
//   const int n_cells=mgcc->get_nef();
//   // coordinates
//   xcoord=new double[n_cells*nvrt];
//   ycoord=new double[n_cells*nvrt];
//   zcoord=new double[n_cells*nvrt];
//   double *ccc = (new double[n_cells]) ;
//   double *cres = (new double[n_cells]) ;
//   tcon=new int[n_cells*nvrt];
//   mgcc->Writefine_xmf(ccc,cres,xcoord,ycoord,zcoord);
//   // coordinate datasets --------------
//   dimsf[0] = n_cells*nvrt;  dimsf[1] = 1;
//
//   name="X1";  print_Dhdf5(file1,name,dimsf,xcoord);
// //   dataspace = H5Screate_simple ( 2, dimsf, NULL );
// //   dataset = H5Dcreate ( file1, "X1",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, xcoord );
// //   H5Dclose ( dataset );
//   name="X2";  print_Dhdf5(file1,name,dimsf,ycoord);
// //   dataset = H5Dcreate ( file1, "X2",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, ycoord );
// //   H5Dclose ( dataset );
//   name="X3";  print_Dhdf5(file1,name,dimsf,zcoord);
// //   dataset = H5Dcreate ( file1, "X3",  H5T_NATIVE_DOUBLE,
// //                         dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
// //   status = H5Dwrite ( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
// //                       H5P_DEFAULT, zcoord );
// //   H5Dclose ( dataset );
// //   H5Sclose ( dataspace );
//   // ccf --------------------------------
//   dimsf[0] = n_cells;  dimsf[1] = 1;
//   name="ccf";  print_Dhdf5(file1,name,dimsf,ccc);
//   name="cres";  print_Dhdf5(file1,name,dimsf,cres);
//   // line connectivity ----------------
//   for(int it=0; it<n_cells*nvrt; it++) tcon[it]=it;
//   dimsf[0] = n_cells;  dimsf[1] = nvrt;
//   name="conn";  print_Dhdf5(file1,name,dimsf,tcon);
//
//   H5Fclose(file1);
//   delete []tcon;  delete []ccc;
//   delete []cres;  delete []xcoord;
//   delete []ycoord;  delete []zcoord;
//
//   *n_l_out = n_lines;  *n_c_out = n_cells;
//   return;
// }
//
// #endif
