#ifndef _domainconf
#define _domainconf

//-------------------------------------------------------------
//Dimension of the problem
// DIMENSION = 2 or 3 (D) 
#define DIMENSION  (2)
// #define AXISYM 
//-----------------------------------------------------------
#define  BDRY_TOLL  1.e-12//tolerance for setting the BCs
// #define MATBC_INTERFACE       // boundary conditions and material
//Physical dimensions of the mesh: 

 
#define NUM_MESH (2)
 
#define HAVE_GROUP (1) // see in the gencase file


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
  #define LXE (1.)
  #define LYB (0.)
  #define LYE (1.)
  #define LZB (0.)
  #define LZE (1.)
#endif

#endif