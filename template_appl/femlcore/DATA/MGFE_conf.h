// *****************************************************************
//      MGGeomElconf
// *****************************************************************

#ifndef _quadtri_cfg
#define _quadtri_cfg

// ======= GEOMETRIC ELEMENT TYPE ====
//Geometric element,
// only quadratic,
// volume and boundary

//here no structure, but the numbers

#include "Domain_conf.h"


//***********
#define ELTYPE 27        //quadrilateral
// #define ELTYPE 10     //triangular 
//***********

#define NUMGEOMELS 1//=========number of geometric element types


//***********
#if ELTYPE == 27

#if DIMENSION==2 
  #define  NNDS        (9)  //Lagrange Quad9
  #define  NNDS_B      (3)  //Lagrange Edge3
#endif
#if DIMENSION==3 
  #define  NNDS       (27)
  #define  NNDS_B      (9)  //Lagrange 
#endif

#endif

#if ELTYPE == 10

#if DIMENSION==2 // ----------------------
#define  NNDS        (6)  //Lagrange Tri6
#define  NNDS_B      (3)  //Lagrange Edge3
#endif

#if DIMENSION==3  // ------------------------ 
#define  NNDS       (10)
#define  NNDS_B      (6)
#endif

#endif

// #endif


// *****************************************************************
//      MGFEconf
// *****************************************************************


// #ifndef _feconf
// #define _feconf

// enum FECouple {
//  Q2Q1=0, Q2Q0, Q2S0, Q2P1 
// }; 
//    static const FECouple myfemc = Q2P1;  
//    static const FECouple myfemc = Q2Q1;      
// #define  Q2Q (myfemc) 

#define  Q2Q1   
// #define  NDOF_K      (0)
// #define  Q2Q0 
// #define  NDOF_K      (1)
// #define  Q2P1
// #define  NDOF_K      (3)
// #define  Q2S0
// #define  NDOF_K      (1)   


#if DIMENSION==2 // ===============================

// Quad9-Quad8-Quad4
#if ELTYPE == 27

#define  NDOF_FEM    (9)
#define  NDOF_FEMB   (3)
#define  NDOF_P  (4)
#define  NDOF_PB (2)



#define  NSUBDOM    (4)  
#endif
// Tri6-Tri3
#if ELTYPE == 10

#define  NDOF_FEM    (6)
#define  NDOF_FEMB   (3)
#define  NDOF_P  (3)
#define  NDOF_PB (2)
#define  NSUBDOM    (3)   
#endif

// 2d axisymmetric case -----------------
//#define AXISYMX
#endif


#if DIMENSION==3 // ========================

// Hex27-Hex20-Hex8
#if ELTYPE==27

#define  NDOF_FEM   (27)
#define  NDOF_FEMB   (9)
#define  NDOF_P      (8)
#define  NDOF_PB     (4)
#define  NDOF_K      (1)
#define  NSUBDOM    (8)  
#define  NDOF_BK      (0)
#endif
// Tetra10-Tetra4
#if ELTYPE==10

#define  NDOF_FEM   (10)
#define  NDOF_FEMB   (6)
#define  NDOF_P  (4)
#define  NDOF_PB (3)
#define  NSUBDOM    (4)   
#endif


#endif



#endif


