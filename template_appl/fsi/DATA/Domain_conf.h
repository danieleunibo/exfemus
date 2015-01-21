#ifndef _domainconf
#define _domainconf

//-------------------------------------------------------------
//Dimension of the problem
// DIMENSION = 2 or 3 (D) 
#define DIMENSION  (3)
// #define AXISYM 
//-----------------------------------------------------------
#define  BDRY_TOLL  1.e-12//tolerance for setting the BCs
#define MATBC_INTERFACE
 
#define NUM_MESH (1)
 
#define HAVE_GROUP (1) // see in the gencase file
 
//Physical dimensions of the mesh: 

#if DIMENSION==2  // =========================
  #define LXB (0.)
  #define LXE (1.)
  #define LYB (0.)
  #define LYE (1.)
  #define LZB (0.)
  #define LZE (0.)
#endif

#if DIMENSION==3 // ========================= 
  #define LXB (0.0)
  #define LXE (2.)
  #define LYB (0.)
  #define LYE (1.)
  #define LZB (0.)
  #define LZE (1.)
#endif

#endif