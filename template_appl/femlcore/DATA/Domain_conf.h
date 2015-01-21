#ifndef _domainconf
#define _domainconf

//-------------------------------------------------------------
//Dimension of the problem
// DIMENSION = 1 or 2 or 3 (D)
#define DIMENSION  (3)
// #define AXISYM
//-----------------------------------------------------------
#define NUM_MESH (2)
#define  BDRY_TOLL  1.e-12//tolerance for setting the BCs
//-----------------------------------------------------------

#define COUPLED_MESH // if defined make the meshes pass data.
#define  BDRY_TOLL  1.e-12//tolerance for setting the BCs


#define GAMBIT_INTERFACE   // boundary conditions and material from Gambit



//Physical dimensions of the mesh:

#if DIMENSION==2  // =============================

//beam turek
#define LXB (0)
#define LXE (1)
#define LYB (0)
#define LYE (1)
#define LZB (0.)
#define LZE (0.)

//     #define LXB (0.0041)
//   #define LXE (.0302)
//    #define LYB (0.)
//    #define LYE (2.)
//   #define LZB (0.)
//   #define LZE (.05)

//cylinder axialsymmetric
//   #define LXB (0)
//   #define LXE (0.005)
//   #define LYB (0.)
//   #define LYE (0.06)
//   #define LZB (0.)
//   #define LZE (0.)

#endif

#if DIMENSION==3 // ========================= 

//   #define LXB -0.1
//   #define LXE .1
//   #define LYB -.1
//   #define LYE .1
//   #define LZB 0.
//   #define LZE .8

//   #define LXB (0.0041)
//   #define LXE (.0302)
//   #define LYB (0.)
//   #define LYE (.2)
//   #define LZB (0.)
//   #define LZE (.05)

#define LXB (0.)
#define LXE (1.)
#define LYB (0.)
#define LYE (1.)
#define LZB (0.)
#define LZE (4.485)

#endif









#endif
