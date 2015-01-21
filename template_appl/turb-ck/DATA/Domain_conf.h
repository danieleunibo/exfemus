#ifndef _domainconf
#define _domainconf

//-------------------------------------------------------------
//Dimension of the problem
// DIMENSION = 2 or 3 (D) 
#define DIMENSION  (3)
//  #define AXISYM
#define NUM_MESH (1)
//-----------------------------------------------------------

#define  BDRY_TOLL  1.e-12 //tolerance for setting the BCs
// #define MATBC_INTERFACE       // boundary conditions and material

#define HAVE_GROUP (1) // see in the gencase file



//Physical dimensions of the mesh: 

#if DIMENSION==2  // =============================
  
  // annulus
// #define LXB (0.00605)
// #define LXE (0.03025)
// #define LXB (0.0041)
// #define LXE (0.03025)

// boundary layer
  #define LXE (0.1)
  #define LXB (0.)
  #define LYB (0.)
  #define LYE (5.)
  #define LZB (0.)
  #define LZE (0.)

  
#endif

#if DIMENSION==3 // ========================= 
  
/*  #define LXB -0.1
  #define LXE .1
  #define LYB -.1
  #define LYE .1
  #define LZB 0.
  #define LZE .8 */ 
  
//     #define LXB (0.0041)
//   #define LXE (.0302)
//    #define LYB (0.)
//    #define LYE (.2)
//   #define LZB (0.)
//   #define LZE (.05)

  #define LXB (0.)
  #define LXE (0.1)
  #define LYB (0.)
  #define LYE (0.1)
  #define LZB (0.)
  #define LZE (5.)
  
#endif









#endif
