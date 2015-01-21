// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"

#ifdef NS_EQUATIONS

// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

// class files --------------------------------------------------
#include "MGSclass_conf.h"   // Navier-Stokes class conf file
#include "MGSolverNS.h"      // Navier-Stokes class header file
#include "ReactData.h"       // Reactor data

// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include -----------------------------------
#include "MGMesh.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsMap.h"  // Equation map class

// standard lib -----------------------------------------------
#include <string.h>          // string library

// local alg lib -----------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ==============================================================

// ==================================================================
/// This routine constructs the FSI class:
MGSolNS::MGSolNS(
  MGEquationsMap& mg_equations_map_in,
  int             nvars_in[],
  std::string     eqname_in,
  std::string     varname_in
):  MGSolDA(mg_equations_map_in,nvars_in,eqname_in,varname_in),
  /// A) reading parameters
  _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes (top level)
  _dt(_mgutils.get_par("dt")),            // parameter  dt
  _uref(_mgphys.get_par("Uref")),         // parameter  u reference
  _lref(_mgphys.get_par("Lref")),         // parameter  l reference
  _rhof(_mgphys.get_par("rho0")),         // parameter density
  _muf(_mgphys.get_par("mu0")) {          // parameter viscosity
  // ================================================================

  /// B) setting class variables
  // class variable names
  _dir=0;
  _mgphys._flag_probe=0;  // to define the surface node to probe
#if NS_EQUATIONS==2       //   segregated ( P in NSP_EQUATIONS)
  if(!varname_in.compare("u")) _dir=0;   // u-equation
  if(!varname_in.compare("v")) _dir=1;   // v-equation
  if(!varname_in.compare("w")) _dir=2;   // w-equation
  _var_names[0]=varname_in;
  _refvalue[0]=_uref;
#else
  _var_names[0]="u";
  _refvalue[0]=_uref; // velocity 2D
  _var_names[1]="v";
  _refvalue[1]=_uref;
#if NS_EQUATIONS==1       //  coupled  (+P)
  _var_names[DIMENSION]="p";
  _refvalue[DIMENSION]=_rhof*_uref*_uref;  // pressure
#endif
#if DIMENSION==3
  _var_names[2]="w";
  _refvalue[2]=_uref;// velocity 3D
#endif
#endif

  /// C ) setting solver type
  for(int l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVERNS);

  /// D) setting nondimensional parameters
  _IRe=_muf/(_rhof*_lref*_uref);          // Reynolds number
  _IFr=9.81*_lref/(_uref*_uref);          // Froud number
  _dirg[0] = _mgphys.get_par("dirgx");    // x-gravity
  _dirg[1] = _mgphys.get_par("dirgy");    // y-gravity
  _dirg[2] = _mgphys.get_par("dirgz");    // z-gravity

  return;
}//  =================================================================


// ====================================================================
/// This function assembles the matrix and the rhs:
//  ===================================================================

void  MGSolNS::GenMatRhs(const double time,
                         const int Level,
                         const int mode) {
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

  /// a) Set up
  // geometry -----------------------------------------------------------------------------------
  const int  ndim = DIMENSION;                                           //dimension
  const int  offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
  const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];                     // element connectivity
  int        el_neigh[NDOF_FEM];                                         // element connectivity
  int        sur_toply[NDOF_FEMB];                                       // boundary topology
  double     xx_qnds[DIMENSION*NDOF_FEM], xxb_qnds[DIMENSION*NDOF_FEMB]; // element node coords
  double     normal[DIMENSION];                                          // normal to the boundary
  double     xm[DIMENSION];

  // Gauss integration ---------------------------------------------------------------------------
  const int  el_ngauss = _fe[2]->_NoGauss1[ndim-1];                //elem gauss points
  const int  elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];          //elem gauss points
  double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];           // Jac, Jac*w Jacobean
  double dphijdx_g[3][DIMENSION];
  double dphiidx_g[3][DIMENSION]; // global derivatives at g point

  // Number of  element dof: constant[0]-linear[1]-quadratic[2] -----------------------------------
  int el_ndof[3];  el_ndof[0]=NDOF_K;              // number of element  dofs   [const(0)]
  for(int ideg=1; ideg<3; ideg++)  el_ndof[ideg]=((_nvars[ideg]>0)? _fe[ideg]->_NoShape[ndim-1]:0);
  int elb_ndof[3];      elb_ndof[0]=NDOF_BK;       // number of element boundary dofs [const(0)]
  elb_ndof[1]=NDOF_PB;  elb_ndof[2]=NDOF_FEMB;            //  [lin(1)+quad(2)]
  int el_mat_nrows =0;                             // Global dof -> Matrix
  for(int ideg=0; ideg<3; ideg++) el_mat_nrows +=_nvars[ideg]*_el_dof[ideg];
  int el_mat_ncols = el_mat_nrows;                        // square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);          // dof index vector


  // fields -> Navier-Stokes + two-meshes -----------------------------------------------------------
  const int ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NS_F]];       // NS equation [NS_F]]
  double vel_g[DIMENSION];  double vel_gdx[DIMENSION*DIMENSION];   // velocity field
  double val_tbg[10];  double h_eff[DIMENSION];                    // turbulence, h_eff
  double Pe_h[DIMENSION];  double f_upwind[DIMENSION];             // local Peclet, upwind

  const int t_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[T_F]];  //Temperature
  // two-meshes
  double mid_t=0;        int mid_t_count=0;                        // temperature
  double mid_p=0;                              // temperature
  double mid_vel=0;                         // velocity
  double mid_pressure=0; int mid_pressure_count=0;                 //  pressure
  double p_old;
  int probe_count=0;

  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) --------------------------
  A[Level]->zero();  if(mode ==1) b[Level]->zero();                 // global matrix+rhs
  DenseMatrixM KeM;  DenseVectorM FeM;                              // local  matrix+rhs

  KeM.resize(el_mat_nrows,el_mat_ncols);  FeM.resize(el_mat_nrows); // resize  local  matrix+rhs


  int nel_lev=0;
  for(int isubd=0; isubd<_iproc; isubd++) {
    nel_lev += _mgmesh._off_el[0][isubd*_NoLevels+Level+1]-_mgmesh._off_el[0][isubd*_NoLevels+Level];//const
  }
  int ndof_lev=0;
  for(int pr=0; pr <_mgmesh._iproc; pr++) {
    int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
    ndof_lev +=delta;
  }
  double u1_in=0; double w1_in=0;
  double p1_in=82000.; // editing pressure to impose the pressure value at inlet (NOT coupled mesh!)
  int und_step=0; double T1_in=608.15; double v1_in=0;
  // ============================================================================
  // ============================================================================
#ifdef COUPLED_MESH // (see Domain_conf.h!)

  if((int)(time/_dt)%4==3) und_step=1;

  // ===============---- Pressure  (primary) reading from file ==================
//   char buf[100];
//   FILE *fp;    fp = fopen("cmsh2.txt", "r");      // open the file
//   fscanf(fp,"%lf,%lf,%lf", &T1_in, &v1_in, &p1_in);              // read from file
//   fclose(fp);                                     // close the file
  u1_in=_mgmesh._mesh2_data[5+0]; v1_in=_mgmesh._mesh2_data[5+1];
  w1_in=_mgmesh._mesh2_data[5+2]; p1_in=_mgmesh._mesh2_data[5+DIMENSION];
  T1_in=_mgmesh._mesh2_data[DIMENSION+1];


//   if(pres<100) pres_mesh2=100000.;
//   pres_mesh2=280000.;
  // dp=7.735062*10500 v\approx 0.8 no-gravity
  // ===============---- Pressure integral (reactor) ===========================
  mid_pressure=0; mid_pressure_count=0;
  const int  nel_e1 = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1];
  const int  nel_b1 = _mgmesh._off_el[0][Level+_NoLevels*_iproc];
  for(int iel=0; iel < (nel_e1 - nel_b1); iel++) {

    // geometry element quantities ---------------------------
    // Element connectivity  and coordinates (xx_qnds)
    _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);  // get connectivity and coord
    _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);// get neighboring element

    // local system fields
    get_el_dof_bc(Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

//     // element nodes coordinates
//     for(int idim=0; idim<DIMENSION; idim++) {
//       xm[idim]=0.;
//       for(int d=0; d< NDOF_FEM; d++) {
//         _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];
//         xm[idim] +=xx_qnds[idim*NDOF_FEM+d];
//       }
//     }
    // pressure linear field ------------------------------------
    for(int id=0; id<_el_dof[1]; id++)    {
      const int  kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+(_nvars[2])*offset]; // dof from top level
      _data_eq[1].ub[id]= ((*x_old[_NoLevels-1])(kdof_top));     // element sol
    } // end linear  ------------------------------------------------

    // Pressure on Boundary
    for(int  iside=0; iside< el_sides; iside++) {

      if(el_neigh[iside] == -1) {

        // setup boundary element -> connectivity+coordinates
        for(int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) {
          int lnode=_mgmesh._GeomEl._surf_top[lbnode+NDOF_FEMB*iside];// local nodes
          sur_toply[lbnode]=lnode;          // lbnode -> lnode
          elb_conn[lbnode]=el_conn[lnode];  // connctivity el_conn->elb_conn
//                  if (_bc_vol[sur_toply[lbnode]+3*NDOF_FEM] ==10) flag_write=1;
//                     for (int idim=0; idim<ndim; idim++) { // coordinates
//                         double xyz=xx_qnds[idim*NDOF_FEM+lnode];
//                     }
        }

        // Pressure average
        if(_bc_vol[sur_toply[NDOF_FEMB-1]+2*NDOF_FEM]==0) { // Dirichlet
          if(_bc_bd[sur_toply[NDOF_FEMB-1]+2*NDOF_FEM]==4) { // Homogeneous
            double sum_node=(_data_eq[1].ub[sur_toply[0]]+_data_eq[1].ub[sur_toply[1]]
                             +_data_eq[1].ub[sur_toply[2]]+_data_eq[1].ub[sur_toply[3]])*0.25;
            mid_pressure+=sum_node;
            mid_pressure_count++;
          }
        }
      }   //end if side
    } //end for iside
  }
  if(mid_pressure_count >0) p_old=mid_pressure/mid_pressure_count;
  printf(" iproc %d from  %d to %d :   p_old  %f = %f / %d press-prim %f \n",
         _iproc, nel_b1,nel_e1,p_old, mid_pressure,mid_pressure_count,p1_in);
  //===============---- end pressure integral =================================

#endif // ---------  end COUPLED_MESH ----------------------
  // ============================================================================





  //element type loop
  const int  nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1];
  const int  nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];


  for(int iel=0; iel < (nel_e - nel_b); iel++) {
    /// b) Element  Loop over the volume (n_elem)
    // set to zero matrix and rhs
    KeM.zero();    FeM.zero();

    // geometry element quantities ------------------------------------------------------------
    // Element connectivity  and coordinates (xx_qnds)
    _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);  // get connectivity and coord
    _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);// get neighboring element

    // local system fields
    get_el_dof_bc(Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

    for(int ivar0=0; ivar0<_nvars[0]; ivar0++) {
      for(int idof0=0; idof0<el_ndof[0]; idof0++) {
        el_dof_indices[el_ndof[2]*_nvars[2]+idof0] =
          _node_dof[Level][idof0+(iel+nel_lev)*el_ndof[0]+(ivar0+_nvars[2])*offset];
      }
    }
    // element nodes coordinates
    for(int idim=0; idim<DIMENSION; idim++) {// dimension loop
      xm[idim]=0.;                               // element middle point xm[idim]
      for(int d=0; d< NDOF_FEM; d++) {           // element point loop
        _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];
        xm[idim] +=xx_qnds[idim*NDOF_FEM+d];     // element middle point
      }
      xm[idim] /= NDOF_FEM;                      // element middle point
      // element nodes coordinates
      // element grid distance
      h_eff[idim]=fabs(xx_qnds[idim*NDOF_FEM+2+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+0]);
      double h_1=fabs(xx_qnds[idim*NDOF_FEM+3+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+1]);
      if(h_eff[idim] <  h_1) h_eff[idim]=h_1;  // Max dx diagonal term
    }
    // external fields (from constant 0 to quadratic 2)
//     for (int deg=0; deg<3; deg++) {
//       for (int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
//         _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-
//                                              _data_eq[deg].indx_ub[eq],el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],
//                                              _data_eq[deg].ub);
//       }
//     }

    ////-----------------------------------------------------------------------------------------------------------------------------------------------------
// #if NS_EQUATIONS%2==0 // pressure as external field (projection or splitting)
//     _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol(0,1,el_ndof[1],el_conn, offset,1,_data_eq[1].ub);
// #else
//     _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_sol(DIMENSION,1,el_ndof[1],el_conn, offset,0,_data_eq[1].ub);
// #endif
    ////-----------------------------------------------------------------------------------------------------------------------------------------------------
    for(int eq=0; eq<_data_eq[2].n_eqs; eq++) {
      _data_eq[2].mg_eqs[eq]->get_el_sol(0,_data_eq[2].indx_ub[eq+1]-_data_eq[2].indx_ub[eq],
                                         _el_dof[2],el_conn,offset,_data_eq[2].indx_ub[eq],_data_eq[2].ub);
    }
//         }
    //linear field
//     _data_eq[1].mg_eqs[0]->get_el_sol(_nvars[2],_data_eq[1].indx_ub[0],_el_dof[1],el_conn,offset,0,_data_eq[1].ub);
    for(int  id=0; id<_el_dof[1]; id++)    {
      const int  kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+(_nvars[2])*offset]; // dof from top level
      _data_eq[1].ub[id]= ((*x_old[_NoLevels-1])(kdof_top));     // element sol
    } // end quadratic ------------------------------------------------

    //  external node quantities -------------------------------------
    _data_eq[0].mg_eqs[0]->get_el_sol_piece(
      _nvars[2]+_nvars[1],_data_eq[0].indx_ub[0],_el_dof[0],iel+ndof_lev,offset,0,_data_eq[0].ub);

    // element fields ---------------------------------------------------------------------
    double phase=1.;    double rho=1.;    double mu=1.;      // non-dimensional properties
    double fpp=0.;   if(xm[2]>=C_IN && xm[2]<=C_OUT) fpp=1.; // distributed pressure losses
#ifdef TBK_EQUATIONS
    _y_dist=_mgmesh._dist[ iel+nel_b];                       // distance from the wall
#endif
    double bbeta=_mgphys._sys_data[1][iel+nel_b];            // accidental pressure losses

    // ------------------------------------------------------------------------------------
    /// c) gaussian integration loop (n_gauss)
    // ------------------------------------------------------------------------------------
    for(int qp=0; qp<  el_ngauss; qp++) {

      // shape functions at gaussian points -----------------------------------------------
      // quadratic continuous
//       det[2]      = _fe[2]->Jac(qp,xx_qnds,InvJac[2]);     // Jacobian
//       JxW_g[2] =det[2]*_fe[2]->_weight1[ndim-1][qp];       // weight
//       _fe[2]->get_phi_gl_g(ndim,qp,_phi_g[2]);
//       _fe[2]->get_dphi_gl_g(ndim,qp,InvJac[2],_dphi_g[2]); // global coord deriv
      // discontinuous and continuous linear
      for(int ideg=0; ideg<2; ideg++) if(_nvars[ideg]>0)
          _fe[ideg]->get_phi_gl_g(ndim,qp,_phi_g[ideg]);  // shape funct

      for(int ideg=1; ideg<3; ideg++) {  // linear-quadratic  [1]=linear [2]=quad
        det[ideg]      = _fe[ideg]->Jac(qp,xx_qnds,InvJac[ideg]);     // Jacobian
        JxW_g[ideg] =det[ideg]*_fe[ideg]->_weight1[ndim-1][qp];       // weight
        _fe[ideg]->get_phi_gl_g(ndim,qp,_phi_g[ideg]);               // shape funct
        _fe[ideg]->get_dphi_gl_g(ndim,qp,InvJac[ideg],_dphi_g[ideg]); // global coord deriv
      }


      // Interpolations with linear and quadratic fields ----------------------------------
      // linear fields -> pressure
      // ------------------------------------------------------------------------------------------------------------
// #if (NS_EQUATIONS%2==0)  // projection and segregated
//       interp_el_sol(_data_eq[1].ub,0,2,_phi_g[1],el_ndof[1],_ub_g[1]);
// #else                   // coupled
//       interp_el_sol(_data_eq[1].ub,0,1,_phi_g[1],el_ndof[1],_ub_g[1]);
//       _ub_g[1][1]=0.;
// #endif
      // -----------------------------------------------------------------------------------------------------------

      // quadratic fields
      interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],
                    _phi_g[2],_el_dof[2],_ub_g[2]);
      // pressure linear field
      // interp_el_sol(_data_eq[1].ub,0,_data_eq[1].indx_ub[0],_phi_g[1],_el_dof[1],_ub_g[1]); // linear
      for(int id=0; id<_el_dof[1]; id++)    {
        const int  kdof_top = _node_dof[_NoLevels-1][ el_conn[id]+(_nvars[2])*offset]; // dof from top level
        _data_eq[1].ub[id]= ((*x_old[_NoLevels-1])(kdof_top));     // element sol
      } // end linear  ------------------------------------------------

      interp_el_gdx(_data_eq[2].ub,ns_idx,DIMENSION,_phi_g[2],el_ndof[2],vel_gdx); // derivatives
#ifdef AXISYM
      JxW_g[2]  *=_ub_g[2][0];      JxW_g[1]  *=_ub_g[2][0];
#endif
      double IRe_eff=_IRe;

#ifdef T_EQUATIONS  // Temperature (quad,_indx_eqs[T_F]) ----------------------------------
      double dens = _rhof;                                // reference density   (673.15K)
      double visc = _muf;                                 // reference viscosity (673.15K)
#ifdef TEMPERATURE_DEPENDENCE_NS // TEMPERATURE DEPENDENCY (see MGSclass_conf.h)
      double Temp_g=_ub_g[2][t_idx]; // local temp
      dens = densityNS(Temp_g);       // density (densityNS(x) in MGSclass_conf.h)
      visc = viscosityNS(Temp_g);    // viscosity (viscosityNS(x) in MGSclass_conf.h)
      rho =  dens/_rhof;              // normlized  local density
      mu =   visc/_muf;               // normlized local viscosity
      IRe_eff=_IRe*rho/mu;
      // printf(" Temp %f; rho %f %f ;  visc %f %f",Temp_g, rho, dens, mu, visc);
      if (Temp_g > T_FREEZE_SOLID && Temp_g < T_FREEZE_SOLID+DELTA_T_FREEZE){ 
	if(fpp==0) fpp=1.e+5*(1+(T_FREEZE_SOLID-Temp_g)/(DELTA_T_FREEZE));
	else if (fpp==1) fpp=1.e+5*(1+(T_FREEZE_SOLID-Temp_g)/(DELTA_T_FREEZE))+1.;
      }
	#endif

      
#endif // ---------------------------------------------------------------------------------

      // velocity[NS_F] -> (quad,_indx_eqs[NS_F+idim]) ------------------------------------
      double mod2_vel=1.e-20;
      _sP=1.e-20;

      for(int idim=0; idim< DIMENSION; idim++) {
        // velocity field  -> gaussian points <- _ub_g[2][ns_idx+idim]
        vel_g[idim]=_ub_g[2][ns_idx+idim];
        mod2_vel +=vel_g[idim]*vel_g[idim];
        // turbulence production term  -> _sP
        for(int jdim=0; jdim< DIMENSION; jdim++) {
          const double sss=vel_gdx[idim+jdim*DIMENSION]+vel_gdx[jdim+idim*DIMENSION];
          _sP += sss*sss;
        }
        // upwind term -> f_upwind
        Pe_h[idim]=0.5*mod2_vel*h_eff[idim]/_IRe;             // local  Peclet
        f_upwind[idim]=UP_WIND_NS*rho*0.5*(1./tanh(Pe_h[idim])-1./Pe_h[idim])*h_eff[idim]/sqrt(mod2_vel);
      }
      mod2_vel =sqrt(mod2_vel);


      // -------------------- Temperature[T_F] -> (quad,_indx_eqs[T_F]) -------------------------------------

      // ----------------- Turbulent viscosity [K_F] -> (quad,_indx_eqs[K_F]) ---------------------------
      _mu_turb=0.;// eff visc turb at g pt
#ifdef TBK_EQUATIONS
      const int k_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[K_F]];  // kappa-equations
      _kappa_g[0]= _ub_g[2][k_idx];
      _kappa_g[1]= _ub_g[2][k_idx+1];
      f_mu(val_tbg);
#endif
      IRe_eff = _IRe*(1.+_mu_turb);//  visc_eff  at g pt


      // ---------------------------------------------------------------------------------------------------
      /// d) Assemblying NS equation
      // ---------------------------------------------------------------------------------------------------
      for(int i=0; i< el_ndof[2]; i++) {  // i-node loop quad ++++++++++++++++++++++++++++++++++++++++++++++
        // set up row i
        for(int idim=0; idim<ndim; idim++) {
          dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
        }

        const double phii_g=_phi_g[2][i];
        // loop u,v,w
        for(int  ivar=0; ivar< _nvars[2]; ivar++)    { // i-var loop quad ++++++++++++++++++++++++++++++++++

          int    indx=i+ivar*el_ndof[2];                       // index
          double dtxJxW_g=JxW_g[2]*_bc_vol[i+ivar*el_ndof[2]]; // jac+bc

          // ---------------------------------------------------------------
          // Only vertical velocity at the core inlet
          if(ivar !=2 &&   xm[2]>1.3   &&   xm[2]<1.4) {
            dtxJxW_g=0.;
            KeM(indx,i+ivar*el_ndof[2])=JxW_g[2];
            FeM(indx)=0.;
          }
          // ---------------------------------------------------------------

          // Assemblying rhs --------------------------------------------------------------
          if(mode == 1)                { 
            FeM(indx)  +=
              dtxJxW_g*(                                           // J rho
                +rho*vel_g[ivar+_dir]*phii_g/_dt                             // time  (v_t)
                +0.981*(rho-1.)*_dirg[ivar+_dir]*phii_g  //*((_ub_g[2][2]>HIN)? 0:1) // gravity (rho.g)
#if (NS_EQUATIONS%2==0)      //       projection segregated          //old pressure
                +rho*(_ub_g[1][0] + _ub_g[1][1])*dphiidx_g[2][ivar+_dir]*(1-_nvars[1])
#endif
              );
          }
          // Assemblying Matrix quad ---------------------------------
          // ---------------------------------------------------------
          for(int j=0; j<el_ndof[2]; j++) {                 // j-loop quad shape functions
            const double phij_g= _phi_g[2][j];              // j-shape (quad) functions

            // set up
            double Lap_g=0.,Adv_g=0.,Div_g= 0;
            // Grids pressure drops using Rheme approach
// 			double aa = 3.5, bb = 73.14, cc = 2.79e+10, dd = 2.6;
// 			double ep = 0.27237;
#ifdef  AXISYM
            Lap_g =(IRe_eff+f_upwind*vel[ivar+_dir]*vel[ivar+_dir])*
                   (1-ivar)*phii_g*phij_g/(_ub_g[2][0]*_ub_g[2][0]);// axysimmetry
#endif
            for(int kdim=0; kdim<ndim; kdim++) {
              dphijdx_g[2][kdim] =_dphi_g[2][j+kdim*el_ndof[2]];
              Adv_g +=vel_g[kdim]*dphijdx_g[2][kdim]*phii_g/R_CORE;   // R_CORE is defined in ReactData.h
              Lap_g +=(IRe_eff+f_upwind[kdim]*vel_g[kdim]*vel_g[kdim])*dphijdx_g[2][kdim]*dphiidx_g[2][kdim];
              Div_g +=_ub_dxg[kdim+kdim*ndim];
            }
            // i-variable diagonal blocks := Fem[ivar],Kem[ivar][ivar]
            KeM(indx,j+ivar*el_ndof[2]) +=dtxJxW_g*rho*(
                                            phij_g*phii_g/_dt  // time
                                            + Adv_g            // advection
                                            + Lap_g            // viscous Laplacian
                                            + (IRe_eff+f_upwind[ivar+_dir]*vel_g[ivar+_dir]*vel_g[ivar+_dir])// upwind
                                            *(dphijdx_g[2][ivar+_dir]*dphiidx_g[2][ivar+_dir])   // viscous tensor

                                            // distributed pressure drop
                                            // R_CORE & D_EQ are defined in ReactData.h
                                            +0.5*fpp*mod2_vel*phij_g*4*0.0791/(R_CORE*R_CORE*D_EQ*sqrt(sqrt(dens
                                                *(mod2_vel/R_CORE)*D_EQ/visc)))*phii_g
                                            // accidental pressure drop
                                            +0.5*bbeta*vel_g[ivar+_dir]*phij_g*phii_g/(R_CORE*R_CORE)

// 					    +0.5*bbeta/0.52*(ep*ep)*fmin(aa+(bb/pow(dens*(mod2_vel/R_CORE)*D_EQ/visc,0.264))+(cc/pow(dens*(mod2_vel/R_CORE)*D_EQ/visc,2.79)),(dd/(ep*ep)))
// 					    *vel_g[ivar+_dir]*phij_g*phii_g/(R_CORE*R_CORE) //  grid spacer pressure loss with Rheme correlation
                                          );

            // i-variable off-diagonal block:= Fem[ivar],Kem[ivar][ivar+1] ++++++++++++++++++
            int idimp1=(ivar+1+_dir)%ndim;
#if NS_EQUATIONS==2                        // splitting  (rhs term) 
            FeM(indx) += -_data_eq[2].ub[ns_idx*NDOF_FEM+j+idimp1*el_ndof[2]]* // sum_j
#else                                      // no splitting  (matrix term)
            KeM(indx,j+idimp1*el_ndof[2]) +=
#endif                                     //  splitting (rhs) or not-splitting (matrix)
                         dtxJxW_g* rho*(IRe_eff+f_upwind[ivar+_dir]*vel_g[ivar+_dir]*vel_g[idimp1])*
                         dphijdx_g[2][ivar+_dir]*dphiidx_g[2][idimp1];

#if DIMENSION==3   // I-VARIABLE OFF-DIAGONAL BLOCK:= Fem[ivar],Kem[ivar][ivar+1] (only 3d)
            int idimp2=(ivar+2+_dir)%ndim; //Fem[ivar],Kem[ivar][ivar+2] ------------
#if NS_EQUATIONS==2                        // splitting    (rhs)      
            FeM(indx)+=  -_data_eq[2].ub[ns_idx*NDOF_FEM+j+idimp2*el_ndof[2]]* // sum j
#else                                     // no splitting (matrix)
            KeM(indx,j+idimp2*el_ndof[2])+=
#endif                                    //  splitting (rhs) or not-splitting (matrix)
                         dtxJxW_g* rho*(IRe_eff+f_upwind[ivar+_dir]*vel_g[ivar+_dir]*vel_g[idimp2])*
                         dphijdx_g[2][ivar+_dir]*dphiidx_g[2][idimp2];

#endif     // -------END I-VARIABLE OFF-DIAGONAL BLOCK (3D)  
          } // end A element matrix quad -quad (end loop on j)-------------------------------


#if NS_EQUATIONS==1      // B^T element matrix ( p*div(v) )---------------------------------
          for(int  jkl=1; jkl<2; jkl++) {   // jkl=0 discontinuous ikl=1 continuous pressure
            for(int  j=0; j<el_ndof[jkl]; j++) {   // j-loop shape linear -------------------
              for(int  jvarl=0; jvarl<_nvars[jkl]; jvarl++) { // j-var-loop ------------------
                double psij_g=_phi_g[1][j];          // j-shape linear function
                KeM(indx,j+ndim*el_ndof[2]+jvarl*el_ndof[1]) +=dtxJxW_g*(
                      -psij_g*dphiidx_g[2][ivar]   //  p*div(v)
#ifdef AXISYM
                      -(1-ivar)*psij_g*phii_g/_ub_g[0] // v/r
#endif   // end-if  AXISYM
                    );
              }       // end loop jvarl-loop quad -----------------------------------------
            }         // end-loop j shape linear ------------------------------------------
          }           // end loop jkl-loop ------------------------------------------------
#endif   // end-if  NS_EQUATIONS==1  B^T element matrix     



        }              // end loop ivar-loop quad +++++++++++++++++++++++++++++++++++++++++
      }                // end loop i-loop quad shape ++++++++++++++++++++++++++++++++++++++

//------------------    QL    ------------------------------------------------------------
#if NS_EQUATIONS==1   // only coupled Assemblying Matrix linear ---------------------------
      for(int  ikl=1; ikl<2; ikl++) {   // ikl=0 discontinuous ikl=1 continuous pressure
        for(int   ivarl=0; ivarl< _nvars[ikl]; ivarl++) { // i-var lin-test-function ******
          for(int   i=0; i< el_ndof[ikl]; i++) {  // i-loop: test functions +++++++++++++++
            // set up row i
            int  indx=i+el_ndof[2]*_nvars[2];                   // row index
            const double psii_g=_phi_g[ikl][i];                 // test function
            if(mode == 1)   FeM(indx)  +=
                JxW_g[2]*psii_g*_ub_g[1][0]*KOMP_NS/_dt;// p_old/dt

            // linear-quadratic Assemblying Matrix ----------------------------------------
            for(int j=0; j<el_ndof[2]; j++)     { // j-loop: quad shape functions
              const double phij_g= _phi_g[2][j];           // shape function
              for(int  jvar=0; jvar< _nvars[2]; jvar++) {  // j-shape variable
                // p-equation
                KeM(indx,j+jvar*el_ndof[2]) +=
                  JxW_g[2]*rho*(                         //  Jac x rho (
                    psii_g*_dphi_g[2][j+jvar*el_ndof[2]] //      div u
#ifdef AXISYM
                    +(1-jvar)*psii_g*phij_g/_ub_g[2][0]  //     + v/r)
#endif
                  );
              }// jvar
            }  // j-shape (quad)  ----------------------------------------------------------

            // linear-linear Assemblying Matrix --------------------------------------------
            for(int j=0; j<el_ndof[1]; j++)  {          // j-loop: linear shape functions
              const double psij_g=_phi_g[1][i];         // j-shape function
              KeM(i+el_ndof[2]*_nvars[2],j+ndim*el_ndof[2])  +=
                JxW_g[2]*psii_g*psij_g*KOMP_NS/_dt; // p/dt
            }       // end j-linear shape linear liner -------------------------------------


          }         // end i-linear test functions +++++++++++++++++++++++++++++++++++++++++
        }           // end ivar end linear *************************************************
      }             // end const-linear-quad variable loop
#endif  // ------------   coupled NS equation ----------------------------------------------


    } // end of the quadrature point qp-loop


    // ====================== end volume (element) =======================================


    // ================================================================================
    //         boundary surface
    // ================================================================================
//     double xx=xx_qnds[0*NDOF_FEM+26];
//     double yy=xx_qnds[1*NDOF_FEM+26];
//     double zz=xx_qnds[2*NDOF_FEM+26];


    for(int  iside=0; iside< el_sides; iside++) {
      if(el_neigh[iside] == -1) {
        // local matrix and rhs
        int flag_write=0;
        // setup boundary element -> connectivity+coordinates
        for(int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) {
          int lnode=_mgmesh._GeomEl._surf_top[lbnode+NDOF_FEMB*iside];// local nodes
          sur_toply[lbnode]=lnode;          // lbnode -> lnode
          elb_conn[lbnode]=el_conn[lnode];  // connctivity el_conn->elb_conn
          if(_bc_vol[sur_toply[lbnode]+3*NDOF_FEM] ==10) flag_write=1;
          for(int idim=0; idim<ndim; idim++) {  // coordinates
            double xyz=xx_qnds[idim*NDOF_FEM+lnode];
            _data_eq[2].ub[idim*NDOF_FEM+lbnode]=xxb_qnds[idim*NDOF_FEMB+lbnode]=xyz;
          }
        }

        // normal
        _fe[2]->normal_g(xxb_qnds,normal);
// #ifdef COUPLED_MESH // (press <- mesh2)	===================================================
//         if(flag_write==1) {
//           double loc_vel=0.;   double loc_temp=0.; double loc_pres=0.;
//           for(int i=0; i< elb_ndof[2]; i++) {
//             for(int idim =0; idim<DIMENSION; idim++) {
//               loc_vel += _data_eq[2].ub[ns_idx*NDOF_FEM+sur_toply[i]+idim*el_ndof[2]]*normal[idim];
//             }
// //             _node_glob_probe[probe_count]=el_conn[sur_toply[i]];
// // 	    probe_count++;
//             loc_temp +=_data_eq[2].ub[t_idx*NDOF_FEM+sur_toply[i]];// printf ("  % d     ", elb_conn[i ] );
//             loc_pres +=_data_eq[1].ub[sur_toply[i]];// printf ("  % d     ", elb_conn[i ] );
//           }
//           mid_vel+=(loc_vel)/ elb_ndof[2];
//           mid_t+=(loc_temp)/ elb_ndof[2];
//           mid_p+=(loc_pres)/ elb_ndof[2];
//           mid_t_count++ ;
// //        printf ("  %f   %d      iel   %d       \n",  sqrt(loc_vel/ elb_ndof[2]), mid_vel_count, iel);
//         }
// #endif // =================================================================================


        for(int ivar=0; ivar< _nvars[2]; ivar++)    {

          // Dirichlet boundary conditions  ***********************
          if(_bc_vol[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] ==0) {

            double det= _fe[2]->JacSur(elb_ngauss-1,xxb_qnds);              // jacobian
            double Ipenalty=det/_dt;                        // Dirichlet bc flag

            for(int i=0; i< elb_ndof[2]; i++) {  // +++++++++++++++++++++++

              const int  indx_sol=sur_toply[i]+(ivar+_dir)*el_ndof[2];//ivar=0;
              const int  indx_row=sur_toply[i]+(ivar)*el_ndof[2];//ivar=0;

              // flag setup (\int_[S bc_normal*(u \dot n)+bc_tg*(u-(u.n)n)+bc_rhs*val] dS)
              int bc_s=(int)_bc_bd[indx_row]; //int bc_v=(int)_bc_vol[indx_row];
              int bc_rhs=((bc_s&4)>>2);  // (1??) nonhomogeneous
              int bc_tg=((bc_s&2)>>1);   // (?1?) tg
              int bc_normal=(bc_s%2);    // (??1) normal

              double ndotu2=0.;  // non-homogeneous normal
              for(int kdim=0; kdim<ndim; kdim++)
                ndotu2 +=normal[kdim]*_data_eq[2].ub[ns_idx*NDOF_FEM+sur_toply[i]+kdim*NDOF_FEM];

              // Assemblying rhs // non-homogeneous ----------------------------
              if(mode == 1)  {

                int indx=(ns_idx+2)*NDOF_FEM;
                double old_vel= _data_eq[2].ub[ns_idx*NDOF_FEM+indx_sol];

#ifdef COUPLED_MESH // (press_mesh2 <- mesh2)
                old_vel=  old_vel+und_step*_dt*(p_old-p1_in)*((ivar==2)?1.:0);
// 		old_vel=-0.7;
                if(old_vel> 0.) old_vel=0.;
#endif

                FeM(indx_row) += bc_rhs*Ipenalty*( // Dirichlet -> non-homogeneous flag (1??)
                                   (1-bc_normal)*(1-bc_tg)*old_vel// _u_old[indx_sol]    // (100) single comp bc
                                   +bc_normal*ndotu2 *normal[ivar+_dir]              // (101)normal Dirichlet bc
                                   +bc_tg*(_data_eq[2].ub[ns_idx*NDOF_FEM+indx_sol]-ndotu2*normal[ivar+_dir])// (110) tg Dirichlet bc
                                 );

                if((1-bc_normal)*(1-bc_tg) ==1 && bc_rhs==1) printf(" %f %d  %f  %f %f \n", old_vel, und_step,p_old,p1_in,
                      _dt*(p_old-p1_in)
                                                                     );
              }

              // Assemblying Matrix ---------------------------------
              // u.n normale tg -(u.n)n
              KeM(indx_row,indx_row) += ((1-bc_normal)*(1-bc_tg)+bc_tg)*Ipenalty // (?00) single comp bc
                                        + Ipenalty*normal[ivar+_dir]*normal[ivar+_dir]*(
                                          +bc_normal   // (?01) Normal Dirichlet bc
                                          -bc_tg       // (?10) tg Dirichlet bc
                                        );

              for(int   jvar=ivar+_dir+1; jvar< ndim+_dir; jvar++)    {
                int      jvar1=jvar%DIMENSION;
#if NS_EQUATIONS==2
                FeM(indx_row) += -1.*Ipenalty*_data_eq[2].ub[(ns_idx+jvar1)*NDOF_FEM+sur_toply[i]]*
#else
                KeM(indx_row,sur_toply[i]+jvar1*el_ndof[2]) += Ipenalty*
#endif
                                 normal[jvar1]*normal[ivar+_dir]*(
                                   +bc_normal   // (?01) Normal Dirichlet bc
                                   -bc_tg       // (?10) tg Dirichlet bc
                                 );
              }// end  +u.n normale tg -(u.n)n

            }// i loop
          } // end if  Dirichlet  ***********************************

          //  Neumann boundary conditions  ***********************************
          else if(_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] !=0) {

            double alpha_t=0.;
#ifdef TBK_EQUATIONS
            double kappa_w=_kappa_g[0]+1.e-20;
            double yplus=rho*0.547722557505*_y_dist*sqrt(kappa_w)/_muf;
            double u_plus=(yplus< .1) ? yplus:2.5*log(9.*yplus);
            std::cout << " \n  yplus "<<  yplus << " kappa " << kappa_w << " y " << _y_dist;
            alpha_t=rho*0.547722557505*sqrt(kappa_w)/u_plus;
#endif

            for(int  qp=0; qp<  elb_ngauss; qp++) {  //gaussian integration loop (n_gauss)

              // quad/linear  [2]=quad [1]=linear------------------------------------
              det[2]     = _fe[2]->JacSur(qp,xxb_qnds);     // local coord _phi_g and jac
              JxW_g[2]=det[2]*_fe[2]->_weight1[ndim-2][qp]; // weight
              _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);    // global coord _phi_g

#ifdef AXISYM   // axisymmetric  (index ->0)
              interp_el_sol(_data_eq[2].ub,0,DIMENSION,_phi_g[2],elb_ndof[2],_ub_g[2]);
              JxW_g[2]  *=_ub_g[0];
#endif
              // ***********************************
              for(int i=0; i< elb_ndof[2]; i++) {  // Assemblying NS equation
                // set up row i
                const double phii_g=_phi_g[2][i];
                const int   indx_sol=sur_toply[i]+(ivar+_dir)*el_ndof[2];// volume dof index
                const int   indx_row=sur_toply[i]+ivar*el_ndof[2];// volume dof index
                // boundary flag
                int bc_s=(int)_bc_bd[indx_row];
                int bc_v=(int)_bc_vol[indx_row];
                double dtJxW_g=JxW_g[2]*((bc_v==0)?0:1);
                int bc_rhs   =((bc_s&4)>>2); // (1??) -> nonhomogeneous
                int bc_tg    =((bc_s&2)>>1); // (?1?) -> tg
                int bc_normal=(bc_s%2);      // (??1) -> normal

                // Assemblying rhs ----------------------------
                if(mode == 1)   {
//                   FeM(indx_row)  += bc_rhs*dtJxW_g*_phi_g[2][i]*normal[ivar+_dir]*(-pres/_refvalue[DIMENSION]);
//                   FeM(indx_row)  += bc_rhs*dtJxW_g*_phi_g[2][i]*((ivar==2)? 1:0)*(-pres/_refvalue[DIMENSION]);

                }

                // Assemblying Matrix ---------------------------------
                for(int j=0; j<elb_ndof[2]; j++) {  // quad -quad
                  const double phij_g= _phi_g[2][j];
                  KeM(indx_row,indx_row) +=dtJxW_g*bc_tg*alpha_t*phij_g*phii_g // (?1?) -> tg
                                           + dtJxW_g*alpha_t*normal[ivar+_dir]*normal[ivar+_dir]*phij_g*phii_g*(
                                             bc_normal // (??1) -> normal
                                             -bc_tg     // (?1?) -> tg
                                           );

                  for(int  jvar=ivar+_dir+1; jvar< ndim+_dir; jvar++)    {  // u.n normale tg -(u.n)n
                    int jvar1=jvar%DIMENSION;
#if NS_EQUATIONS==2
                    FeM(indx_row) +=   -_data_eq[2].ub[ns_idx*NDOF_FEM+sur_toply[j]+jvar1*el_ndof[2]]*
#else
                    KeM(indx_row,sur_toply[j]+jvar1*el_ndof[2]) +=
#endif
                                       dtJxW_g*alpha_t*normal[jvar1]*normal[ivar+_dir]*(
                                         +bc_normal   // (?01) Normal Dirichlet bc
                                         -bc_tg       // (?10) tg Dirichlet bc
                                       );
                  }

                }//  j-loop
              }// i-loop

            }// end gaussian  integration
          }  //  end if Neumann  ***********************************

        } // end if ivar   +++++++++++++++++++++++++++++++
      } //end if side
    } //end for iside

    // ==========================================================================
    // ====================== end boundary =====================================
    // ==========================================================================

    /// e) Add them to the global matrix ------------------------
    A[Level]->add_matrix(KeM,el_dof_indices);
    if(mode == 1) b[Level]->add_vector(FeM,el_dof_indices);

  } // end of element loop

  // clean ---------------------------------------------------
  el_dof_indices.clear();
  A[Level]->close();
  if(mode == 1) b[Level]->close();
//   A[Level]->print();

// #ifdef COUPLED_MESH  // print on file --------------------------
//   printf(" probe_count %d  \n",probe_count);
//   for(int   jvar=0; jvar< probe_count; jvar++)    {
// 
//     printf("  %d %d \n", jvar,_mgphys._node_glob_probe[jvar]);
// 
//   }
// 
// 
// 
// 
//   if(mid_t_count>0) {
//     FILE *cmyfile;  cmyfile = fopen("cmsh1.txt","w");
// // 		fprintf(cmyfile, "%f,%f\nTemperature  velocity ", mid_t/*/mid_t_count*/, mid_vel/*/mid_t_count*/);
//     fprintf(cmyfile, "%f,%f,%f\n Temperature,  Velocity, Pressure ",
//             mid_t/mid_t_count ,mid_vel/mid_t_count,mid_p/mid_t_count);
//     fclose(cmyfile);
//   }
// #endif             // end print on file --------------------------
#ifdef PRINT_INFO
  std::cout << " GenMatRhs(NS): Assembled  " << Level << " (Level) "  << std::endl;
#endif

  return;
}





// ========================================================================
void  MGSolNS::get_nod_probe(
) {// ===============================================================

  /// a) Set up
  // geometry -----------------------------------------------------------------------------------
  const int  ndim = DIMENSION;                                           //dimension
  const int  offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
  const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
  int        el_conn[NDOF_FEM], elb_conn[NDOF_FEMB];                     // element connectivity
  int        el_neigh[NDOF_FEM];                                         // element connectivity
  int        sur_toply[NDOF_FEMB];                                       // boundary topology
  double     xx_qnds[DIMENSION*NDOF_FEM], xxb_qnds[DIMENSION*NDOF_FEMB]; // element node coords
  double     normal[DIMENSION];                                          // normal to the boundary
  double     xm[DIMENSION];

  // Gauss integration ---------------------------------------------------------------------------
  const int  el_ngauss = _fe[2]->_NoGauss1[ndim-1];                //elem gauss points
  const int  elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];          //elem gauss points

  // Number of  element dof: constant[0]-linear[1]-quadratic[2] -----------------------------------
  int el_ndof[3];  el_ndof[0]=NDOF_K;              // number of element  dofs   [const(0)]
  for(int ideg=1; ideg<3; ideg++)  el_ndof[ideg]=((_nvars[ideg]>0)? _fe[ideg]->_NoShape[ndim-1]:0);
  int elb_ndof[3];      elb_ndof[0]=NDOF_BK;       // number of element boundary dofs [const(0)]
  elb_ndof[1]=NDOF_PB;  elb_ndof[2]=NDOF_FEMB;            //  [lin(1)+quad(2)]
  int el_mat_nrows =0;                             // Global dof -> Matrix
  for(int ideg=0; ideg<3; ideg++) el_mat_nrows +=_nvars[ideg]*_el_dof[ideg];
  int el_mat_ncols = el_mat_nrows;                        // square matrix
  std::vector<int> el_dof_indices(el_mat_ncols);          // dof index vector


  int probe_count=0;
  int Level=_NoLevels-1;


  for(int kproc=0; kproc<_mgmesh._n_subdom; kproc++) {


    int nel_lev=0;
    for(int isubd=0; isubd<kproc; isubd++) {
      nel_lev += _mgmesh._off_el[0][isubd*_NoLevels+Level+1]-_mgmesh._off_el[0][isubd*_NoLevels+Level];//const
    }
    int ndof_lev=0;
    for(int pr=0; pr <kproc; pr++) {
      int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
      ndof_lev +=delta;
    }
    //element type loop
    const int  nel_e = _mgmesh._off_el[0][Level+_NoLevels*kproc+1];
    const int  nel_b = _mgmesh._off_el[0][Level+_NoLevels*kproc];


    for(int iel=0; iel < (nel_e - nel_b); iel++) {
      /// b) Element  Loop over the volume (n_elem)

      // geometry element quantities ------------------------------------------------------------
      // Element connectivity  and coordinates (xx_qnds)
      _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds,kproc);  // get connectivity and coord
//      _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);// get neighboring element
      for(int  iside=0; iside<el_sides; iside++)    {
        el_neigh[iside] = _mgmesh._el_neighbor[0][(iel+_mgmesh._off_el[0][Level+_NoLevels*kproc])*el_sides+iside];
      }

      // local system fields
      get_el_dof_bc(Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

//     for(int ivar0=0; ivar0<_nvars[0]; ivar0++) {
//       for(int idof0=0; idof0<el_ndof[0]; idof0++) {
//         el_dof_indices[el_ndof[2]*_nvars[2]+idof0] =
//           _node_dof[Level][idof0+(iel+nel_lev)*el_ndof[0]+(ivar0+_nvars[2])*offset];
//       }
//     }
      // element nodes coordinates
      for(int idim=0; idim<DIMENSION; idim++) {// dimension loop
        xm[idim]=0.;                               // element middle point xm[idim]
        for(int d=0; d< NDOF_FEM; d++) {           // element point loop
//         _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];
          xm[idim] +=xx_qnds[idim*NDOF_FEM+d];     // element middle point
        }
        xm[idim] /= NDOF_FEM;                      // element middle point
      }

      // ====================== end volume (element) =======================================

      // ================================================================================
      //         boundary surface
      // ================================================================================

      
      
 

      
     
      
      for(int  iside=0; iside< el_sides; iside++) {
        if(el_neigh[iside] == -1) {
          // local matrix and rhs
          int flag_write=0;
          // setup boundary element -> connectivity+coordinates
          for(int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) {
            int lnode=_mgmesh._GeomEl._surf_top[lbnode+NDOF_FEMB*iside];// local nodes
            sur_toply[lbnode]=lnode;          // lbnode -> lnode
            elb_conn[lbnode]=el_conn[lnode];  // connctivity el_conn->elb_conn
            if(_bc_vol[sur_toply[lbnode]+3*NDOF_FEM] ==10) flag_write=1;
          }

          // (press <- mesh2)	===================================================
          if(flag_write==1) {
	    
        // normal
            for(int i=0; i< elb_ndof[2]; i++) {
              _mgphys._node_glob_probe[probe_count]=elb_conn[i];
              probe_count++;

            }
          }

        } //end if side
      } //end for iside

    } // end of element loop

  }
// print on file --------------------------
 for(int ivar=0; ivar< 10; ivar++)  _mgmesh._mesh2_data[ivar]=0.;
 
//   printf(" NS probe_count %d  \n",probe_count);
//   for(int   jvar=0; jvar< probe_count; jvar++)    {
//     printf(" %d  ",  _mgphys._node_glob_probe[jvar]);
//     if(jvar%20==0)  printf(" \n");
//   }
//   printf("  \n");
  _mgphys._flag_probe=probe_count;
#ifdef PRINT_INFO
  std::cout << "  MGSolNS::get_nod_probe: node count defined  " << Level << " (Level) with "  << _mgphys._flag_probe  << " (nodes) "  << std::endl;
#endif

  return;
}





// ========================================================================================
/// This function compute average surface  values on  over mesh[0] and
///  transfer data to mesh[1]
void MGSolNS::print_ext_data(
  double mesh2mesh_data[] // transfer vector from mesh[1]
) { // =====================================================================================
  int flag_probe=_mgphys._flag_probe;
  const int Level=_NoLevels-1;
  const int offset= _mgmesh._NoNodes[_NoLevels-1];

  // velocity (outlet from reactor) data
  double sum_u=0.;
  
//   for(int ivar=0; ivar< DIMENSION; ivar++) {
   
    for(int inode=0; inode< flag_probe; inode++) {
      double u_c=(*x_old[Level])(_node_dof[Level][_mgphys._node_glob_probe[inode]+0*offset])*_refvalue[0];
      double v_c=(*x_old[Level])(_node_dof[Level][_mgphys._node_glob_probe[inode]+1*offset])*_refvalue[1];
       double w_c=(*x_old[Level])(_node_dof[Level][_mgphys._node_glob_probe[inode]+2*offset])*_refvalue[2];
       if(fabs(u_c)>fabs(v_c)){sum_u +=   u_c*0.866025403784+  v_c*0.5+w_c*0. ;   }
       else{sum_u +=   v_c*0.866025403784+  u_c*0.5+w_c*0. ;}
//       sum_u[ivar]  += (*x_old[Level])(_node_dof[Level][_mgphys._node_glob_probe[inode]+ivar*offset])*_refvalue[ivar];
    }
     sum_u /=flag_probe; 
     _mgmesh._mesh2_data[0]=sum_u;// writing on mesh[0] normal dir
     mesh2mesh_data[0]=sum_u;     // writing on mesh[1] normal dir
      _mgmesh._mesh2_data[2]=0.;// writing on mesh[0]   tg 1 dir
     mesh2mesh_data[2]=0;     // writing on mesh[1]     tg 1 dir
      _mgmesh._mesh2_data[1]=0.;// writing on mesh[0]   tg 2 dir
     mesh2mesh_data[1]=0;     // writing on mesh[1]     tg 2 dir
//   }
 
  // pressure (outlet from reactor) data
  double press=0.;
  for(int inode=0; inode< flag_probe; inode++) {
    if(inode%9 <4)
      press  += (*x_old[Level])(_node_dof[Level][_mgphys._node_glob_probe[inode]+DIMENSION*offset])*_refvalue[DIMENSION];
  }
  // writing data in the system vector
  _mgmesh._mesh2_data[DIMENSION]=9.*press/(4.*flag_probe); // writing on mesh[0]
  mesh2mesh_data[DIMENSION]=9.*press/(4.*flag_probe);      // writing on mesh[1]

 
#ifdef PRINT_INFO
  std::cout << " probe NS computed  " << Level << " (Level) with "  << flag_probe  << " (nodes) "  << std::endl;
  
  std::cout << " 0n NS probe_count " << flag_probe << std::endl;
  std::cout << " vel = " << mesh2mesh_data[0] << "  " << 
        mesh2mesh_data[1] << "  "<< mesh2mesh_data[2] << "  "<<std::endl;
  std::cout << " pressure = " <<  mesh2mesh_data[DIMENSION]  << std::endl;
#endif

  return;
}







// =========================================================================================
/// This function controls the assembly and the solution of the NS_equation system:
void MGSolNS::MGTimeStep(
  const double time,  // time
  const int /*iter*/  // Number of max inter
) {
// =========================================================================================

/// A) Set up the time step
  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;
  if(_mgphys._flag_probe==0) get_nod_probe();

  /// B) Assemblying of the Matrix-Rhs
#if PRINT_TIME==1
  std::clock_t start_time=std::clock();
#endif
  GenMatRhs(time,_NoLevels-1,1);                                                // matrix and rhs
  for(int Level = 0 ; Level < _NoLevels-1; Level++)   GenMatRhs(time,Level,0);  // matrix
#if PRINT_TIME==1
  std::clock_t end_time=std::clock();
  std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif

  /// C) Solution of the linear MGsystem (MGSolNS::MGSolve).
  MGSolve(1.e-6,40);
#if PRINT_TIME==1
  end_time=std::clock();
  std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
            << "s "<< std::endl;
#endif

  /// D) Update of the old solution at the top Level  (MGSolNS::OldSol_update),
//  x_old[_NoLevels-1]->localize(*x_oold[_NoLevels-1]);

  x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);
  return;
}// =======================================================================================




#ifdef TBK_EQUATIONS
// =====================================================
void  MGSolNS::f_mu(double val[]) {
  // Turbulent viscosity
  if(_kappa_g[0]< 1.e-20) _kappa_g[0]= 1.e-20;  // kappa
  double tau_k=1.;  // turbulent time

#if (TBK_EQUATIONS==0)   //kappa  (Prandtl length)--------------------------
  tau_k=_y_dist/sqrt(_kappa_g[0]);
#endif  // end kappa -------------------------------------------------------

#if (TBK_EQUATIONS/2==1)  // kappa-epsilon (TBK_EQUATIONS=2 or 3) ---------- 
  tau_k= CMU*_kappa_g[0]/_kappa_g[1]; // tau=1/omega=CMU*kappa/epsilon
#endif   // kappa-epsilon --------------------------------------------------

#if (TBK_EQUATIONS/2==2)       // kappa-omega ------------------------------
  if(_kappa_g[1]> 1.)   tau_k=1/_kappa_g[1];  // tau=1/omega

#endif  // -----------------------------------------------------------------
  // turbulent viscosity
  double Re_t= _kappa_g[0]*tau_k/_IRe;
  // limits
  if(Re_t > MU_TOP)  Re_t =MU_TOP;
  if(Re_t < MU_LOW)  Re_t =MU_LOW;
  _mu_turb=Re_t;

  // double Re_t= exp(_kappa_g[0]-_kappa_g[1])*_rhof/_muf;
//   _mu_turb=Red_t;// val[4]=
//          double Re_t=_ub_g[_tb_idx]*_rhof/_muf;
//      if(Re_t < 1.) Re_t=1.; if (Re_t > MU_TOP) Re_t=MU_TOP;
//      const double alphastar_t=(0.024+ 0.166666666667* Re_t)/(1.+0.166666666667* Re_t);
//      Re_t *=alphastar_t;
#ifdef LOWRE
  // eff visc val[6](k) val[7](w) (all eq) -------------------
  _mu_turb= Re_t*(0.024+0.166666666667*Re_t)/(1.+0.166666666667*Re_t);
  val[4]

#endif
#ifdef SST
  // Phi2 - F2 ---------------
  double Phi2=500.*_muf/(_rhof*_y_dist*_y_dist*_omega_g);
  const double Phib=2.*sqrt(_kappa_g)/(0.09*_omega_g*_y_dist);
  if(Phi2 < Phib) Phi2=Phib;
  const  double F2=tanh(Phi2*Phi2);

  // eff visc coeff  fact_mu (all eq)
  double fact_mu=1;
  double fact_mub=0.31*_omega_g/(sqrt(0.5*_sP)* F2+1.e-20);
  if(fact_mub< 1.) fact_mu = fact_mub;
//      {
//        if(fact_mu > 10*fact_mub) {
// //        std::cout << fact_mu << "  bound corr " << fact_mub << "  " << std::endl;
//       fact_mu *=0.1;}
//       else
//     }
  _mu_turb=  Re_t*fact_mu(0.024+ 0.166666666667* Re_t)/(1.+0.166666666667* Re_t); // val[4]=

#endif

  return;
}

#endif

#if (NS_EQUATIONS%2==0)
// ===============================================================
// --------------   PRESSURE EQUATION [P_F] ------------------
// ===============================================================
// ==============================================================
// NS_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// NS_EQUATIONS==1 coupled    solver (u,v,w,p)
// NS_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================



// ==================================================================
/// This function constructs the 3d-2D MGSolNSP class I
MGSolNSP::MGSolNSP(MGEquationsMap& mg_equations_map_in,
                   std::string eqname_in,
                   std::string varname_in):
  MGSolDA(mg_equations_map_in,0,1,eqname_in,varname_in),
  /// A) reading parameters
  _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes
  _dt(_mgutils.get_par("dt")),       // parameter  dt
  _rhof(_mgphys.get_par("rho0")),    // parameter  density reference
  _uref(_mgphys.get_par("Uref")),    // parameter  vel reference
  _lref(_mgphys.get_par("Lref")) {
  //  ===============================================================
  /// B) setting class variables
  _var_names[0]="p";
  _refvalue[0]=_rhof*_uref*_uref;    // class variable names

  /// C ) setting solver type
  for(int  l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVER_NSP);

  /// D) setting nondimensional parameters
  return;
}//  ================================================================




// ===============================
// Initial and boundary conditions
// ===============================
/// This function generates the initial conditions for the energy equation:
// =================================================
void MGSolNSP::ic_read(double xp[],double u_value[]) {
// xp[]=(xp,yp) u_value[]=(u,v,p)
  u_value[0] =0.;
}

// ========================================================
/// This function  defines the boundary conditions for the energy equation:
void MGSolNSP::bc_read(double xp[],int bc_Neum[],int bc_flag[]) {
  // =================================================
  const double Lref = _mgphys.get_par("Lref");
  double ILref = 1./Lref;
#if DIMENSION==2
// xp[]=(xp,yp) bc_flag[T]=0-Dirichlet 1-Neumann
  //    boundary conditions box

  if(xp[0]< LXB*ILref+BDRY_TOLL) {
    bc_flag[0]=0;      // left
    bc_Neum[0]=1;
  }
  if(xp[0]> LXE*ILref-BDRY_TOLL) {
    bc_flag[0]=0;      // right
    bc_Neum[0]=1;
  }
  if(xp[1]< LYB*ILref+BDRY_TOLL) {
    bc_flag[0]=0;      // bottom
    bc_Neum[0]=1;
  }
  if(xp[1]> LYE*ILref-BDRY_TOLL) {
    bc_flag[0]=0;      // top
    bc_Neum[0]=0;
  }
#endif

#if DIMENSION==3
// =================================================
  // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
  //    boundary conditions box
  if(xp[0] < LXB*ILref+BDRY_TOLL) {
    bc_Neum[0]=1;
    bc_flag[0]=0;
  }
  if(xp[0] > LXE*ILref-BDRY_TOLL) {
    bc_Neum[0]=1;
    bc_flag[0]=0;
  }
  if(xp[1] < LYB*ILref+BDRY_TOLL) {
    bc_Neum[0]=1;
    bc_flag[0]=0;
  }
  if(xp[2] < LZB*ILref+BDRY_TOLL) {
    bc_Neum[0]=1;
    bc_flag[0]=0;
  }
  if(xp[2] > LZE*ILref-BDRY_TOLL) {
    bc_Neum[0]=1;
    bc_flag[0]=0;
  }
  if(xp[1] > LYE*ILref-BDRY_TOLL) {
    bc_Neum[0]=0;
    bc_flag[0]=0;
  }
#endif
  return;
} // end boundary conditions ==========================


//  ====================================================
/// This function assembles the matrix and the rhs:
//  ====================================================
void  MGSolNSP::GenMatRhs(
  const double /* time*/, // time  <-
  const int  Level,  // Level <-
  const  int mode    // mode  <- (1=rhs) (0=only matrix)
) {  // ===============================================

  ///A) Set up
  // geometry ---------------------------------------------------------------------------------------
  const int   ndim = DIMENSION;                      //dimension
  const int   offset = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
  double      xx_qnds[DIMENSION*NDOF_FEM];           // element node coords
  int         el_conn[NDOF_FEM];                    // element connectivity

  // gauss integration  ------------------------------------------------------------------------------
  double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];           // Jac, Jac*w Jacobean
  double dphijdx_g[3][DIMENSION];
  double dphiidx_g[3][DIMENSION]; // global derivatives at g point
  const int  el_ngauss = _fe[2]->_NoGauss1[ndim-1];                //elem gauss points

  // element dofs: costant[0]-linear[1]-quadratic[2]-------------------------------------------------
  int el_ndof[3];
  el_ndof[0]=1;                  // number of el dofs
  int el_mat_nrows =0;                            // number of mat rows (dofs)
  for(int ideg=1; ideg<3; ideg++) {
    el_ndof[ideg]=_fe[ideg]->_NoShape[ndim-1];
    el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
  };
  int el_mat_ncols = el_mat_nrows;                     //square matrix
  std::vector<int > el_dof_indices(el_mat_ncols);      // element dof vector

  //  fields -> Navier-Stokes  [NS_F]
  int ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NS_F]];

  // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) ---------------------------
  A[Level]->zero();
  if(mode ==1) b[Level]->zero();   // global matrix+rhs
  DenseMatrixM KeM;
  DenseVectorM FeM;                // local  matrix+rhs
  KeM.resize(el_mat_nrows,el_mat_ncols);
  FeM.resize(el_mat_nrows);        // resize  local  matrix+rhs


  const int  nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
  const int  nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
  for(int iel=0; iel < (nel_e - nel_b); iel++) {

    ///B) Element  Loop over the volume (n_elem)
    // set to zero matrix and rhs and center
    KeM.zero();
    FeM.zero();

    // geometry and element  fields ------------------------------------
    // Element Connectivity (el_conn)  and coordinates (xx_qnds)
    _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);

    // local system fields
    get_el_dof_bc(Level,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

    // element nodes coordinates
    for(int  d=0; d< NDOF_FEM; d++)
      for(int  idim=0; idim<DIMENSION; idim++)
        _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];

    // external fields
    for(int deg=0; deg<3; deg++) {
      for(int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
        _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],
                                             el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub);
      }
    }


    /// c) gaussian integration loop (n_gauss)
    for(int  qp=0; qp< el_ngauss; qp++) {  // +++++++++++++++++++++++++++++++++


      // shape functions at gaussian points -----------------------------------
      for(int  ideg=1; ideg<3; ideg++) {  // linear-quadratic  [1]=linear [2]=quad
        det[ideg]      = _fe[ideg]->Jac(qp,xx_qnds,InvJac[ideg]);     // Jacobian
        JxW_g[ideg] =det[ideg]*_fe[ideg]->_weight1[ndim-1][qp];       // weight
        _fe[ideg]->get_phi_gl_g(ndim,qp,_phi_g[ideg]);               // shape funct
        _fe[ideg]->get_dphi_gl_g(ndim,qp,InvJac[ideg],_dphi_g[ideg]); // global coord deriv
      }

      //external fields (quad only)
      interp_el_gdx(_data_eq[2].ub,ns_idx,DIMENSION,_dphi_g[2],el_ndof[2],_ub_dxg);

      /// d) Local (element) assemblying pressure equation
      for(int  i=0; i<el_ndof[1]; i++)     {
        // set up row i if(xp[1]<0.01) {
//         u_value[1] =4.;
//         if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) {  u_value[0] =.1; }
        const double phii_g=_phi_g[1][i];
        double Div_g=0.;
        for(int  idim=0; idim< ndim; idim++)  {
          dphiidx_g[1][idim]=_dphi_g[1][i+idim*el_ndof[1]];
          Div_g +=_ub_dxg[idim+idim*DIMENSION];      //divergence of velocity field
        }
        double dtxJxW_g=JxW_g[2]*_bc_vol[i];

        // Rhs Assemblying  ----------------------------
        if(mode == 1)  FeM(i) += -1.*dtxJxW_g*Div_g*phii_g/_dt;

        // Matrix Assemblying ---------------------------
        for(int  j=0; j<el_ndof[1]; j++) {

          double Lap=0.;
          for(int  idim=0; idim< ndim; idim++) {
            dphijdx_g[1][idim]=_dphi_g[1][j+idim*el_ndof[1]];
            Lap += dphijdx_g[1][idim]*dphiidx_g[1][idim]; // Laplacian
          }

          // pressure matrix assembling
          KeM(i,i) += 1-_bc_vol[i];
          KeM(i,j) += dtxJxW_g*_bc_vol[j]*Lap;
        } // ---------------------------------------------
      }
    } // end of the quadrature point qp-loop +++++++++++++++++++++++++

    /// e) Global assemblying pressure equation
    A[Level]->add_matrix(KeM,el_dof_indices);                  // global matrix
    if(mode == 1)   b[Level]->add_vector(FeM,el_dof_indices);  // global rhs
// 	_data_eq[0].mg_eqs[0]->get_el_oldsol(_nvars[2]+_nvars[1],_nvars[0],_el_dof[0],iel+ndof_lev,offset,0,_data_eq[0].ub);

  } // end of element loop
  // clean and close
  el_dof_indices.clear();
  A[Level]->close();
  if(mode == 1)  b[Level]->close();


// 	   	_data_eq[0].mg_eqs[0]->get_el_sol_piece(_nvars[2]+_nvars[1],_data_eq[0].indx_ub[0],_el_dof[0],iel+ndof_lev,offset,0,_data_eq[0].ub);
// int aaa[1];




// std::cout<< "\n Average Temp        ***** "<< mid_t/mid_t_count<<"   ***";
#ifdef PRINT_INFO
  std::cout<< " Matrix Assembled(P)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif

  return;
}

// =========================================================================================
/// This function controls the assembly and the solution of the P_equation system:
void MGSolP::MGTimeStep(
  const double time,  // time
  const int /*iter*/  // Number of max inter
) {
// =========================================================================================

/// A) Set up the time step
  std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;

  /// B) Assemblying of the Matrix-Rhs
#if PRINT_TIME==1
  std::clock_t start_time=std::clock();
#endif
  GenMatRhs(time,_NoLevels-1,1);                                                // matrix and rhs
  for(int Level = 0 ; Level < _NoLevels-1; Level++)   GenMatRhs(time,Level,0);  // matrix
#if PRINT_TIME==1
  std::clock_t end_time=std::clock();
  std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif

  /// C) Solution of the linear MGsystem (MGSolP::MGSolve).
  MGSolve(1.e-6,40);
#if PRINT_TIME==1
  end_time=std::clock();
  std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
            << "s "<< std::endl;
#endif

  /// D) Update of the old solution at the top Level  (MGSolP::OldSol_update),
  x[_NoLevels-1]->localize(*x_oold[_NoLevels-1]);  // dp
  x_old[_NoLevels-1]->add(1.,*x_oold[_NoLevels-1]);// p

  return;
}// =======================================================================================



#endif // ENDIF NS_EQUATIONS==0
#endif  //ENDIF NS_EQUATIONS

