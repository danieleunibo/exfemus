#ifndef __solverlibconf_h__
#define __solverlibconf_h__


// ======== SOLVER LIBRARY ======
//   #define HAVE_LASPACKM
    #define HAVE_PETSCM 
    #define HAVE_CUSP
    #define HAVE_MPI   //what if I want to use Petsc without MPI?

#define LSOLVER LASPACK_SOLVERSM  
#ifdef HAVE_PETSCM
#undef LSOLVER
 #define LSOLVER  PETSC_SOLVERSM
#endif



//******PETSC VERSION ***88
#ifdef HAVE_PETSCM

//was: #include "libmesh_config.h"   //TODO REMOVE THIS, as it depends on libmesh
//Instead of libmesh_config we put only the variables 
// assuming that the rest of the code picks the LIBMESH variables
//DIRECTLY from somewhere else


/* PETSc's major version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_MAJOR 
#define LIBMESH_DETECTED_PETSC_VERSION_MAJOR  3 
#endif

/* PETSc's minor version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_MINOR 
#define LIBMESH_DETECTED_PETSC_VERSION_MINOR  2
#endif

/* PETSc's subminor version number, as detected by LibMesh */
#ifndef LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR 
#define LIBMESH_DETECTED_PETSC_VERSION_SUBMINOR  0 
#endif

#endif



#endif