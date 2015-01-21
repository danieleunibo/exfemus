
//c++ includes
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

//*************************
// level
// coarse element axial length
// number of elements in the core section
#define NZC 10

#define D_EQ 0.0129     // equivalent diameter fo the fuel bundle subchannel
#define R_CORE 0.5408   // R_CORE = v/v* (A_channel/A_core_section)

// Geometry

// half fuel assembly length on section
#define HR 0.147
//#define HR   0.5

////////////// ELSY PARAMETERS /////////////////
#define HIN  (2.25)  // active-core inlet height
#define HOUT (3.15)  // active-core inlet height

#define C_IN (1.3)      // core inlet height
#define C_OUT (3.245)    // core outlet height

// #define HEIGHT_MAX_GV (9)       // vapor generator height (m)
// #define HEIGHT_MIN_GV (8.2)     // inlet/outlet vapor generator height (m)
// #define Q_GV          (10.)     // gv heat absorbed per unit length
// #define T_OUT_GV      (673.15)  // T out vapour generator (K)

// #define H_GV (12.485)  // 13.485-0.4-0.6        gv outflow height: case h0
// #define H_GV (11.245)  // 13.485-0.4-1.24-0.6    gv outflow height: case h1
#define H_GV (8.245)  // 13.485-0.4-1.24-3-0.6  gv outflow height: case h2

#define LONG_GV (0.6)      // gv dimension (height)
#define LONG_PUMP (0.3)      // pump dimension (height)
#define Q_GV (10.)      // gv heat absorbed per unit length
#define T_OUT_GV (673.15)   // T out vapour generator (K)

#define Z_LOW (0.01) // lower point of the tube mesh

////////////////////////////////////////////////
// freezing simulation constants
#define T_FREEZE_SOLID (598.52)
#define DELTA_T_FREEZE (1.63)

////////////// CORE PARAMETERS /////////////////
// #define HIN  0.95  // active-core inlet height
// #define HOUT 1.85  // active-core inlet height
// 
// #define C_IN   0.      // core inlet height
// #define C_OUT  1.945   // core outlet height
// 
// #define TOP  1.85     // 2mesh data passing (MGSolverNS.C --> search: "write file for other mesh")
// #define S_X  0.8085   //        "
// #define D_X  1.3965   //        "
////////////////////////////////////////////////

//////////// CORE ALTO PARAMETERS //////////////
// #define HIN  2.25  // active-core inlet height
// #define HOUT 3.15  // active-core inlet height
// 
// #define C_IN   1.3      // core inlet height
// #define C_OUT  3.245    // core outlet height
// 
// #define TOP  4.36      // 2mesh data passing (MGSolverNS.C --> search: "write file for other mesh")
// #define S_X  0.8085    //        "
// #define D_X  1.3965    //        "
////////////////////////////////////////////////

////////////// BOX PARAMETERS /////////////////
// #define HIN  2.25  // active-core inlet height
// #define HOUT 3.15  // active-core inlet height
// 
// #define C_IN  1.3       // core inlet height
// #define C_OUT 3.245     // core outlet height
// 
// #define TOP  4.3   // 2mesh data passing (MGSolverNS.C --> search: "write file for other mesh")
// #define S_X  0.4   //         "
// #define D_X  0.6   //         "
////////////////////////////////////////////////

// PRESSURE LOSS
// #define IN_OUT_LOSS  // to activate in & out core abrupt area change pressure drop
// #define GRIDS_LOSS_ON   // to activate grids spacer pressure drop
#define GRIDS_LOSS_OFF   // to deactivate grids spacer pressure drop

// ---------------------------
// zone
// --------------------------

#define INNER 0
#define INTER 1
#define OUTER 2
#define CTRLR 3
#define DUMMY 4


#define CTL_R 0.
//*******************************

#include "Domain_conf.h"
#include "MGFE_conf.h"
#include "Printinfo_conf.h"
class MGMesh;
class MGSystem;
class MGUtils;


class ReactData {


public:
    ReactData() {};
    ~ReactData() {};

    static const int mat_zone[8][8];
    static const double mat_pf[8][8];// Power factor
    static const double axpf[10][3];
    static const double axlpf[10+4][3];// axial pressure loss
    static const double axllf[10+4];



    void  datagen(const MGMesh & mesh_in,   const MGSystem & phys_in, const MGUtils & _mgutils);
    int zone(  double xm, double ym ) ;
    double axial_lf( double z);
    double lossf(  double x,  double y, double z) ;
    double axial_pf(  double x, double y,  double z);
    double asby_power(  double xm,  double ym);


};




