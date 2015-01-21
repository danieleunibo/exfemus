#ifndef __mggauss_h__
#define __mggauss_h__

#include        "Domain_conf.h"

#include "Printinfo_conf.h"

#include <fstream>
#include <vector>
#include "MGUtils.h" 


/// Class containing mathematical information about the finite element.
/** 
 * 
 * */
class MGFE
#ifdef LM_REFCOUNT
      : public ReferenceCountedObject<MGFE>
#endif
{

public:

//   const MGUtils& _mgutils; ///< MGUtils pointer

  ///@{ \name FINITE ELEMENT DEFINITION 
  
  const int _dim; ///< Dimension (3 or 2) used in the 3D FEM
  
  /// Order of the shape functions (2=phi=QUADRATIC, 1=psi=LINEAR)
  const int _order;
  ///@}

  ///@{ \name FEM INTEGRATION 
  int _NoGauss1[3]; ///< Number of Gaussian points in 1-2-3D (for example, for HEX27 ngauss[3]=(3,9,27))
  
  /// Number of shape functions in 1-2-3D
  int _NoShape[3]; 
  /// Polinomial degree
  int _deg[3];
  ///@}
  
  ///@{ \name  SHAPES [0]=1D,[1]=2D,[2]=3D in Gaussian points
  double* _weight1[3];      ///< Weight
  double* _phi_map1[3];     ///< Shape functions
  double* _dphidxez_map1[3]; ///< Shape derivative functions in gaussian points
  
  /// Shape derivative functions in nodal points
  double* _dphidxez_map1_nodes[3];
  ///@}

//------------------------------------------------------------
  ///@{ \name CONSTRUCTOR-DESTRUCTOR 
  MGFE(
//    const MGUtils& mgutils_in,  ///< MGUtils pointer
   const int order,          ///< Order of the shape functions
   
   /// Fem type: 27(Quad,Hex) or 10(Tri,Tet)
   const int fem_type         
  );  
  ~MGFE();
  ///@}
  
  ///@{ \name INITIALIZING AND CLEANING
  void init_qua();    ///< Generates the Lagrangian quad shape functions
  void init_lin();    ///< Generates the Lagrangian linear shape functions
  void init_pie();    ///< Generates the Lagrangian piecewise shape functions
  
  /// Clear data substructures
  void clear();  
  ///@}

  
  ///@{ \NAME RETURN SHAPE DERIVATIVE FUNCTIONS IN GAUSS POINTS
   ///< \param[in]  <dim>    Dimension
   ///< \param[in]  <qp>     Gaussian point
   ///< \param[in]  <InvJac> Jacobean
   ///< \param[out] <dphi>   Derivative
  
  void get_dphi_gl_g(const int dim,const int qp,const  double InvJac[], double dphi[]);
  void get_dphi_gl_g(const int dim,const int qp,const  double InvJac[], double dphi[],int sdim);
  void get_dphi_gl_g(const int dim,const int qp, const double InvJac[],std::vector<double> &dphi);
  void get_dphi_node(const int dim,const int qp,const  double InvJac[], double dphi[]);
  ///@}
  
  ///@{ \NAME RETURN SHAPE FUNCTIONS IN GAUSS POINTS
   ///< \param[in]  <dim>    Dimension
   ///< \param[in]  <qp>     Gaussian point
   ///< \param[out] <phi>    Shape function
  
  void get_phi_gl_g(const int dim, const int qp,double phi[]);
  void get_phi_gl_g(const int dim, const int qp, std::vector<double>& phi);
  ///@}
      
  //  Jacobian functions ---------------------------------
  // Jacobian at gaussian points
  double Jac(const int ng, double x[], double InvJac[]) ;  ///< Jacobian (3D-2D)
  double Jac1D(const int ng,double x[], double InvJac[]);///< Jacobian (1D)
  double Jac2D(const int ng,double x[], double InvJac[]);  ///< Jacobian (2D)
  double Jac3D(const int ng,double x[], double InvJac[]);///< Jacobian (3D)
  // Jacobian at nodal points
  double Jac_nodes(const int ng, double x[], double InvJac[]) ;  ///< Jacobian (3D-2D)
  double Jac1D_nodes(const int ng,double x[], double InvJac[]);///< Jacobian (1D)
  double Jac2D_nodes(const int ng,double x[], double InvJac[]);  ///< Jacobian (2D)
  double Jac3D_nodes(const int ng,double x[], double InvJac[]);///< Jacobian (3D)

  // Jacobian dim -1
  double JacSur(const int ng, double x[], double InvJac[]) const;      ///< Boundary Jacobian (2D-1D)
  double JacSur2D(const int ng, double x[], double InvJac[]) const;
  double JacSur3D(const int n_gauss,double x[], double InvJac[]) const;///< Surface Jacobian (2D)
//   double JacSur2D(const int ng,double x[]) const;     ///<  Line Jacobian    (1D)
  double JacSur1D(const int ng,double x[], double InvJac[]) const;     ///<  P Jacobian    (1D)
  // functions -----------------------------------------
  /// Compute normal
  void normal_g(const double* xx, double* normal_g) const;
 
  // Reading - writing -------------------------------------
  /// Write
  void write(const std::string& name,const int kdim);
private:
  /// Reading
  void read_c(std::istream& infile, const int kdim);
  /// Writing
  void write_c(std::ostream& infile, const int kdim);
};

#endif
