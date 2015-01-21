// ===============================================================
// --------------   NAVIER-STOKES system [NS_F] ------------------
// ===============================================================
#include "Equations_conf.h"
#ifdef FSI_EQUATIONS
// ==============================================================
// FSI_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// FSI_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSI_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// ===============================================================

// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverFSI.h"       // Navier-Stokes class header file


// config file -------------------------------------------------
#include "MGFE_conf.h"        // FEM approximation
#include "MGGeomEl.h"        // FEM approximation
#include "Printinfo_conf.h"  // Print options

// local Femus class include -----------------------------------
// #include "MGMesh.h"          // Mesh class
#include "MeshExtended.h"
#include "MGFE.h"          // Mesh class
#include "MGSystem.h"        // System class
#include "MGEquationsSystem.h"  // Equation map class
#include "EquationSystemsExtendedM.h"  // Equation map class

// standard lib -----------------------------------------------
#include <string.h>          // string library

// local alg lib -----------------------------------------------
#include "dense_matrixM.h"   // algebra dense matrices
#include "sparse_matrixM.h"  // algebra sparse matrices
#include "dense_vectorM.h"   // algebra dense vectors
#include "numeric_vectorM.h" // algebra numerical vectors
#include "linear_solverM.h"  // algebra solvers
// ==============================================================
#ifdef HAVE_MED
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDLoader.hxx"
#include "InterfaceFunctionM.h"
#endif
// ==================================================================
/// This routine constructs the FSI class:
MGSolFSI::MGSolFSI(
    MGEquationsSystem& mg_equations_map_in,
    int             nvars_in[],
    std::string     eqname_in,
    std::string     varname_in
):  MGSolDA(mg_equations_map_in,nvars_in,eqname_in,varname_in),
    /// A) reading parameters
    _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes (top level)
    _dt(_mgutils.get_par("dt")),       // parameter  dt
    _uref(mg_equations_map_in.get_par("Uref")),    // parameter  u reference
    _lref(mg_equations_map_in.get_par("Lref")),    // parameter  l reference
    _rhof(mg_equations_map_in.get_par("rho0")),    // parameter density
    _muf(mg_equations_map_in.get_par("mu0")),
    _rhos(mg_equations_map_in.get_par("rhos")),    // parameter density
    _ni(mg_equations_map_in.get_par("nis")),
    _Emod(mg_equations_map_in.get_par("Es")),
    _hs(mg_equations_map_in.get_par("hs"))
{   // parameter viscosity
    // ================================================================

    /// B) setting class variables
    // class variable names
    _dir=0;
#if FSI_EQUATIONS==2       //   segregated ( P in NSP_EQUATIONS)
    if (!varname_in.compare("u")) _dir=0;  // u-equation
    if (!varname_in.compare("v")) _dir=1;  // v-equation
    if (!varname_in.compare("w")) _dir=2;  // w-equation
    _var_names[0]=varname_in;
    _refvalue[0]=_uref;
#else
    _var_names[0]="u";
    _refvalue[0]=_uref; // velocity 2D
    _var_names[1]="v";
    _refvalue[1]=_uref;
#if FSI_EQUATIONS==1       //  coupled  (+P)
    _var_names[DIMENSION]="p";
    _refvalue[DIMENSION]=_rhof*_uref*_uref;  // pressure
#endif
#if DIMENSION==3
    _var_names[2]="w";
    _refvalue[2]=_uref;// velocity 3D
#endif
#endif

    /// C ) setting solver type
    for (int l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVER_FSI);



    _dirg[0] = mg_equations_map_in.get_par("dirgx");    // x-gravity
    _dirg[1] = mg_equations_map_in.get_par("dirgy");    // y-gravity
    _dirg[2] = mg_equations_map_in.get_par("dirgz");    // z-gravity

//     std::cout << " Reynolds " << _IRe << " mus " <<  _mus << "\n";
    /// D) setting nondimensional parameters
    _IRe=_muf/(_rhof*_lref*_uref);          // Reynolds number
    _IFr=9.81*_lref/(_uref*_uref);          // Froud number
    _lambda=(_ni*_Emod)/(_rhos*(1+_ni)*(1-2*_ni));//.001;
    _mus=_Emod/(_rhos*2*(1+_ni)*_lref*_uref);

    std::cout << "lambda  " << _lambda << " mus " <<  _mus << "\n";
    return;
}//  =================================================================





// ====================================================================
/// This function assembles the matrix and the rhs:
//  ===================================================================
void  MGSolFSI::GenMatRhs(const double time, const int
                         Level,const  int mode) {
// FSI_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// FSI_EQUATIONS==1 coupled    solver (u,v,w,p)
// FSI_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
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
    double     normal[DIMENSION];
    double    mu_m;                                          // normal to the boundary

    // Gauss integration ---------------------------------------------------------------------------
    const int  el_ngauss = _fe[2]->_NoGauss1[ndim-1];                //elem gauss points
    const int  elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];          //elem gauss points
    double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];           // Jac, Jac*w Jacobean
    double dphijdx_g[3][DIMENSION];
    double dphiidx_g[3][DIMENSION]; // global derivatives at g point
  double dumping=1.030;
    
    
    
    // Number of  element dof: constant[0]-linear[1]-quadratic[2] -----------------------------------
    int el_ndof[3];
    el_ndof[0]=NDOF_K;
    int elb_ndof[3];
    elb_ndof[0]=NDOF_K; // number of el dofs
    int el_mat_nrows =0;                                            // number of mat rows (dofs)
    for (int ideg=1; ideg<3; ideg++) {                              //     ...
        el_ndof[ideg]=((_nvars[ideg]>0)?    _fe[ideg]->_NoShape[ndim-1]:0);                    //   computing
        elb_ndof[ideg]=((_nvars[ideg]>0)?_fe[ideg]->_NoShape[ndim-2]:0);                  //     ...
        el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
    };
    el_mat_nrows +=  el_ndof[0]*_nvars[0];


    int el_mat_ncols = el_mat_nrows;                                //square matrix
    std::vector<int> el_dof_indices(el_mat_ncols);                  // element dof vector

    // fields -> Navier-Stokes ----------------------------------------------------------------------


    int ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[FS_F]];      // NS equation [NS_F]]
        int sdx_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[SDSX_F]];      // NS equation [NS_F]]

    
     double vel_g[DIMENSION];
    double vel_gdx[DIMENSION*DIMENSION];   // velocity field
    double sx_gdx[1*DIMENSION];   // velocity field
    double sy_gdx[1*DIMENSION];   // velocity field
//   double sx_g[DIMENSION]; double sx_gdx[DIMENSION];   // velocity field
//   double sy_g[DIMENSION]; double sy_gdx[DIMENSION];   // velocity field
//   double sz_g[DIMENSION]; double sz_gdx[DIMENSION];   // velocity field
    //Cauchy components
    double sig_1_lin[DIMENSION]; //double sig_1_cons[DIMENSION];
    double sig_1_qua[DIMENSION];
    double sig_3_lin[DIMENSION]; //double  sig_3_cons[DIMENSION];
    double sigma_rhs[DIMENSION];
    for(int j=0; j < DIMENSION; j++) {
        sig_1_lin[j]=0;
        sig_1_qua[j]=0;
        //sig_3_cons[j]=0;  sig_1_cons[j]=0;
        sig_3_lin[j]=0;
        vel_gdx[j]=0;
        vel_gdx[3+j]=0;
        
//         vel_gdx[6+j]=0;
    }
    for(int i=0; i < DIMENSION; i++) for(int j=0; j < DIMENSION; j++)vel_gdx[j+i*DIMENSION]=0.;
    
    
//     double vel_g[DIMENSION];
//     double vel_gdx[DIMENSION*DIMENSION];   // velocity field
    double kappa_mg[2];
    double val_tbg[10];
    double h_eff[DIMENSION];                    // turbulence, h_eff
    double Pe_h[DIMENSION];
    double f_upwind[DIMENSION];             // local Peclet, upwind
    double src_value[DIMENSION];
 double bc_value[DIMENSION];
   int imesh=atoi(_mgutils.get_file("MESHNUMBER").c_str());
   
   
#ifdef  HAVE_MED
//     int interface_id=1;
//     EquationSystemsExtendedM *ext_es=dynamic_cast<EquationSystemsExtendedM *>(&_mgeqnmap);
//     MeshExtended  *ext_mesh=dynamic_cast<MeshExtended *>(&_mgmesh);
//     std::vector<int> mat=ext_mesh->_mat_id;
//     std::vector<int> bc_id=ext_mesh->_bc_id;
//     InterfaceFunctionM * source= ext_es->get_interface_fun(interface_id);
//     int * map_mg ;
//     int * map_med;
//     int n_map;
//     if(imesh==1) {
//         map_mg  = source->get_map_mg();
//         map_med = source->get_map_med();
//         n_map = source->get_n();
//     }

#endif

    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) --------------------------
    A[Level]->zero();
    if (mode ==1) b[Level]->zero();                // global matrix+rhs
    DenseMatrixM KeM;
    DenseVectorM FeM;                              // local  matrix+rhs
    KeM.resize(el_mat_nrows,el_mat_ncols);
    FeM.resize(el_mat_nrows); // resize  local  matrix+rhs


    int ndof_lev=0;
    for (int pr=0; pr <_mgmesh._iproc; pr++) {
        int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
        ndof_lev +=delta;
    }


    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
    for (int iel=0; iel < (nel_e - nel_b); iel++) {

        // set to zero matrix and rhs and center
        KeM.zero();
        FeM.zero();

        // geometry and element  fields ------------------------------------
        // Element Connectivity (el_conn)  and coordinates (xx_qnds)
        _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);
        _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);

        // set element-nodes variables  bc (bc_q_dofs)
        get_el_dof_bc(Level,iel+ndof_lev,_el_dof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);



        // element nodes coordinates
        for (int idim=0; idim<DIMENSION; idim++) {
            for (int d=0; d< NDOF_FEM; d++)  _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];
            // element grid distance
            h_eff[idim]=fabs(xx_qnds[idim*NDOF_FEM+2+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+0]);
            double h_1=fabs(xx_qnds[idim*NDOF_FEM+3+4*(DIMENSION-2)]-xx_qnds[idim*NDOF_FEM+1]);
            if (h_eff[idim] <  h_1) h_eff[idim]=h_1; // Max dx diagonal term
        }
        // external fields (from constant 0 to quadratic 2)
        for (int deg=0; deg<3; deg++) {
            for (int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
                _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-
                                                     _data_eq[deg].indx_ub[eq],el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],
                                                     _data_eq[deg].ub);
            }
        }
//        _data_eq[2].mg_eqs[0]->get_el_sol(DIMENSION,1,el_ndof[1],el_conn, offset,0,_data_eq[1].ub);
        ////-----------------------------------------------------------------------------------------------------------------------------------------------------
// #if FSI_EQUATIONS%2==0 // pressure as external field (projection or splitting)
//     _data_eq[1].mg_eqs[_data_eq[1].tab_eqs[P_F]]->get_el_oldsol(0,1,el_ndof[1],el_conn, offset,1,_data_eq[1].ub);
// #else
//     _data_eq[2].mg_eqs[_data_eq[2].tab_eqs[NS_F]]->get_el_sol(DIMENSION,1,el_ndof[1],el_conn, offset,0,_data_eq[1].ub);
// #endif
        ////-----------------------------------------------------------------------------------------------------------------------------------------------------

        // element fields ----------------------------------
        double phase=1.;
        double ls_interface=0.;   
        double rho =  1.;
        for(int i=0; i<NDOF_FEM; i++)  if(_bc_vol[i] > 3.5 && _bc_vol[i] < 9.5) {
                ls_interface=1.;
            }
        phase=(_bc_vol[NDOF_FEM-1]<1.5)?0:1;
#ifdef TBK_EQUATIONS

#ifdef DIST_FIX
        _y_dist=DIST_FIX;
#else
        _y_dist=_mgmesh._dist[ iel+nel_b];
#endif


#endif
        kappa_mg[0] =0;
        kappa_mg[1] =0;
        mu_m =0.;
        for(int i=0 ; i<ndim; i++) src_value[i]=0.;
#ifdef  HAVE_MED
//         if(imesh==1) {
//             if (mat[iel+nel_b]==interface_id) {
//                 int i_mg =0;
//                 int i_med=-1;
//                 for(int i_mg =0; i_mg<n_map; i_mg++) {
//                     if (map_mg[i_mg] ==  iel+nel_b)
//                     {
//                         i_med=i_mg;
//                         break;
//                     }
//                 }
//                 if(i_med >-1) {
//                     int elem_med= map_med[i_med];
//                     source->eval(elem_med, ndim, src_value);
//                 }
//             }
// 
//         }
#endif
// //==============================================================================================
// //==============================================================================================
// //===============================  __ _       _     _ ==========================================
// //=============================== / _| |_   _(_) __| |==========================================
// //===============================| |_| | | | | |/ _` |==========================================
// //===============================|  _| | |_| | | (_| |==========================================
// //===============================|_| |_|\__,_|_|\__,_|==========================================
// //==============================================================================================
// //==============================================================================================
// //==============================================================================================






        if (phase==0) {
//      if (ls_interface>0.4)
//              std::cout<<"I'm here";
            for (int qp=0; qp<  el_ngauss; qp++) {

                // shape functions at gaussian points --------------------------------------------------------------
                // quadratic continuous
                det[2]      = _fe[2]->Jac(qp,xx_qnds,InvJac[2]);     // Jacobian
                JxW_g[2] =det[2]*_fe[2]->_weight1[ndim-1][qp];       // weight
                _fe[2]->get_phi_gl_g(ndim,qp,_phi_g[2]);
                _fe[2]->get_dphi_gl_g(ndim,qp,InvJac[2],_dphi_g[2]); // global coord deriv
                // discontinuous and continuous linear
                for (int ideg=0; ideg<2; ideg++) if (_nvars[ideg]>0) _fe[ideg]->get_phi_gl_g(ndim,qp,_phi_g[ideg]);  // shape funct



                // linear and quadratic fields ---------------------------------------------------------------------
                // linear fields -> pressure
                // ------------------------------------------------------------------------------------------------------------
// #if (FSI_EQUATIONS%2==0)  // projection and segregated
//       interp_el_sol(_data_eq[1].ub,0,2,_phi_g[1],el_ndof[1],_ub_g[1]);
// #else                   // coupled
//       interp_el_sol(_data_eq[1].ub,0,1,_phi_g[1],el_ndof[1],_ub_g[1]);
//       _ub_g[1][1]=0.;
// #endif
                // --------------------------------------------------------------------------------------------

                // quadratic fields
                interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs], _phi_g[2],el_ndof[2],_ub_g[2]);                                // field
// displacement fields


// //   interp_el_sol(_data_eq[2].ub, sdx_idx,DIMENSION, _phi_g[2],el_ndof[2],_disp_g);    //actual displacement
//                 interp_el_sol(disp,0,DIMENSION, _phi_g[2],el_ndof[2],_disp_g);    //oold displacement


                // field
                interp_el_gdx(_data_eq[2].ub,ns_idx,DIMENSION,_dphi_g[2],el_ndof[2],vel_gdx); // derivatives
//         interp_el_gdx(_data_eq[2].ub,sdx_idx,1,_dphi_g[2],el_ndof[2],sx_gdx); // derivatives
//         interp_el_gdx(_data_eq[2].ub,sdx_idx+1,1,_dphi_g[2],el_ndof[2],sy_gdx); // derivatives
//       interp_el_gdx(_data_eq[2].ub,sdx_idx+2,1,_phi_g[2],el_ndof[2],sz_gdx); // derivatives
#ifdef AXISYM
                JxW_g[2]  *=_ub_g[2][0];
                JxW_g[1]  *=_ub_g[2][0];
#endif
                double IRe_eff=_IRe;


                // velocity[NS_F] -> (quad,_indx_eqs[NS_F+idim]) ---------------------------------------------------
                double mod2_vel=1.e-20;
                _sP=1.e-20;
double _ale_vel_g[3];
                for (int idim=0; idim< DIMENSION; idim++) {
                    // velocity field  -> gaussian points <- _ub_g[2][ns_idx+idim]
                    vel_g[idim]=_ub_g[2][ns_idx+idim];
                    mod2_vel +=vel_g[idim]*vel_g[idim];
                    _ale_vel_g[idim]=0.;//(_disp_g[idim]/_dt); //ale velocity
//        if (_ale_vel_g[idim]>1e-1)std::cout<<std::endl<<"\033[1;31m  iel \033[0m "<<iel<<" idim "<<idim<<" alevel  "<< _ale_vel_g[idim];
                    // turbulence production term  -> _sP
                    for (int jdim=0; jdim< DIMENSION; jdim++) {
                        vel_gdx[idim+jdim*DIMENSION]=0.;
                        const double sss=vel_gdx[idim+jdim*DIMENSION]+vel_gdx[jdim+idim*DIMENSION];
                        _sP += sss*sss;
                    }

                }

                mod2_vel =sqrt(mod2_vel);
                // upwind term -> f_upwind
                for (int idim=0; idim< DIMENSION; idim++) {
                    Pe_h[idim]=0.5*mod2_vel*h_eff[idim]/_IRe;             // local  Peclet
                    f_upwind[idim]=/*UP_WIND_FSI*/1.*rho*0.5*(1./tanh(Pe_h[idim])-1./Pe_h[idim])*h_eff[idim]/mod2_vel;
                }
                // -------------------- Temperature[T_F] -> (quad,_indx_eqs[T_F]) -------------------------------------
#ifdef T_EQUATIONS
#ifdef TEMP_CONST
                double Temp_g=_ub_g[_indx_ub[ _indx_eqs[1]]];
                rho *= _mgphys.adensity(Temp_g);
                mu *= _mgphys.aviscosity(Temp_g);
#endif
#endif
                // ----------------- Turbulent viscosity [K_F] -> (quad,_indx_eqs[K_F]) ---------------------------
                _mu_turb=0.;// eff visc turb at g pt
#ifdef TBK_EQUATIONS
                const int k_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[K_F]];  // kappa-equations
                _kappa_g[0]= _ub_g[2][k_idx];
                _kappa_g[1]= _ub_g[2][k_idx+1];
                f_mu(val_tbg);
                kappa_mg[0] +=_kappa_g[0];
                kappa_mg[1] +=_kappa_g[1];
#endif
                IRe_eff = _IRe*(1.+_mu_turb);//  visc_eff  at g pt
                mu_m +=IRe_eff;
                // ---------------------------------------------------------------------------------------------------


                /// d) Assemblying NS equation
                for (int i=0; i< el_ndof[2]; i++) { // +++++++++++
                    // set up row i
                    for (int
                            idim=0; idim<ndim; idim++) {
                        dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
                    }

                    const double phii_g=_phi_g[2][i];
                    // loop u,v,w
                    for (int  ivar=0; ivar< _nvars[2]; ivar++)    {

                        int    indx=i+ivar*el_ndof[2];//ivar=0;
//        double interface_node=1.; //check if is an interface node
//        if (_bc_vol[i+ivar*NDOF_FEM]==5) {interface_node=0.;}
                        double dtxJxW_g=JxW_g[2]*(_bc_vol[i+ivar*NDOF_FEM]%2)*_dt;


                        // Assemblying rhs ----------------------------
                        if (mode == 1)                {
                            FeM(indx)  +=  dtxJxW_g*rho*(
                                               vel_g[ivar+_dir]*phii_g/_dt +((ivar==0)? 1:0)*phii_g  *0 // time
//                              + rho*_IFr*_dirg[ivar+_dir]*phii_g  // x-gravity
// #if (FSI_EQUATIONS%2==0)      //       projection segregated
//                              +(_ub_g[1][0] + _ub_g[1][1])*dphiidx_g[2][ivar+_dir]*(1-_nvars[1])  //old pressure
// #endif
                                           );
                        }
                        // Assemblying Matrix quad ---------------------------------
                        // ---------------------------------------------------------
                        for (int j=0; j<el_ndof[2]; j++) {
                            const double phij_g= _phi_g[2][j];

                            // set up
                            double Lap_g=0.,Adv_g=0.,Div_g= 0;
#ifdef  AXISYM
                            Lap_g =2.*(1-ivar)*IRe_eff*
                                   phii_g*phij_g/(_ub_g[2][0]*_ub_g[2][0]);// axysimmetry
#endif
                            for (int kdim=0; kdim<ndim; kdim++) {
                                dphijdx_g[2][kdim] =_dphi_g[2][j+kdim*el_ndof[2]];
                                Adv_g +=(vel_g[kdim]-_ale_vel_g[kdim])*dphijdx_g[2][kdim]*phii_g;
                                Lap_g +=(IRe_eff+f_upwind[kdim]*vel_g[kdim]*vel_g[kdim])*(dphijdx_g[2][kdim]*dphiidx_g[2][kdim]);
//              Lap_g +=(1)*(dphijdx_g[2][kdim]*dphiidx_g[2][kdim]);
                                Div_g +=_ub_dxg[kdim+kdim*ndim];
                            }
                            int idimp1=(ivar+1+_dir)%ndim; // block +1 [2-6-7] ---------------
                            int idimp2=(ivar+2+_dir)%ndim; // block +2 [3-4-8] ------------
                            // diagonal blocks [1-5-9]
                            KeM(indx,j+ivar*el_ndof[2]) +=dtxJxW_g*rho*(
                                                              phij_g*phii_g/_dt  // time
                                                              + Adv_g            // advection
                                                              + Lap_g            // viscous Laplacian
//                                             +(dphijdx_g[2][ivar]*dphiidx_g[2][ivar]+dphijdx_g[2][idimp1]*dphiidx_g[2][idimp1]
//                                            +dphijdx_g[2][idimp2]*dphiidx_g[2][idimp2])

                                                              + (IRe_eff+0*f_upwind[ivar+_dir]*vel_g[ivar+_dir]*vel_g[ivar+_dir])// upwind
                                                              *(dphijdx_g[2][ivar+_dir]*dphiidx_g[2][ivar+_dir])   // viscous tensor
                                                          );


#if FSI_EQUATIONS==2    // splitting 
                            FeM(indx) += -_data_eq[2].ub[ns_idx*NDOF_FEM+j+idimp1*el_ndof[2]]*//_u_old[j+idimp1*el_ndof[2]]*
#else       // no splitting
                            KeM(indx,j+idimp1*el_ndof[2]) +=
#endif
                                         dtxJxW_g* rho*(IRe_eff+0*f_upwind[ivar+_dir]*vel_g[ivar+_dir]*vel_g[idimp1])*
                                         dphijdx_g[2][ivar+_dir]*dphiidx_g[2][idimp1];
//                          dtxJxW_g* rho*_dt*(
//                          _mus* dphijdx_g[2][ivar+_dir]*dphiidx_g[2][idimp1]+_lambda*dphijdx_g[2][idimp1]*dphiidx_g[2][ivar]);//sollid
#if DIMENSION==3
                            FeM(indx) += -0*_data_eq[2].ub[(sdx_idx)*NDOF_FEM+j+idimp1*el_ndof[2]]*
                                         dtxJxW_g* rho*(
                                             _mus* dphijdx_g[2][ivar+_dir]*dphiidx_g[2][idimp1]+_lambda*dphijdx_g[2][idimp1]*dphiidx_g[2][ivar]);

#if FSI_EQUATIONS==2    // splitting         
                            FeM(indx)+=  -_data_eq[2].ub[ns_idx*NDOF_FEM+j+idimp2*el_ndof[2]]*//          _u_old[j+idimp2*el_ndof[2]]*
#else        // no splitting
                            KeM(indx,j+idimp2*el_ndof[2])+=
#endif
                                         dtxJxW_g* rho*(IRe_eff+0*f_upwind[ivar+_dir]*vel_g[ivar+_dir]*vel_g[idimp2])*
                                         dphijdx_g[2][ivar+_dir]*dphiidx_g[2][idimp2];
//                          dtxJxW_g* rho*_dt*(
//                         _mus* dphijdx_g[2][ivar+_dir]*dphiidx_g[2][idimp2]+_lambda*dphijdx_g[2][idimp2]*dphiidx_g[2][ivar]);//solid equation

#endif                // ----------------------------------------------------------------------

                            FeM(indx)+= 0* -_data_eq[2].ub[(sdx_idx)*NDOF_FEM+j+idimp2*el_ndof[2]]*//          _u_old[j+idimp2*el_ndof[2]]*
                                        dtxJxW_g* rho*(
                                            _mus* dphijdx_g[2][ivar+_dir]*dphiidx_g[2][idimp2]+_lambda*dphijdx_g[2][idimp2]*dphiidx_g[2][ivar]);
                        } // end A element matrix quad -quad (end loop on j)---------------------------------------------

                        // ------------------------------------------------------------------
#if FSI_EQUATIONS==1      // B^T element matrix ( p*grad(v) )--------------------
                        for (int  ikl=0; ikl<2; ikl++) {  // ikl=0 discontinuous ikl=1 continuous pressure
                            for (int  ivarl=0; ivarl<_nvars[ikl]; ivarl++) {
                                for (int  j=0; j<el_ndof[ikl]; j++) {
                                    double psij_g=_phi_g[ikl][j];
                                    KeM(indx,j+ndim*el_ndof[2]) +=dtxJxW_g*(  // MPascal
                                                                      -psij_g*dphiidx_g[2][ivar]
#ifdef AXISYM
                                                                      -(1-ivar)*psij_g*phii_g/_ub_g[2][0]
#endif
                                                                  );
                                } // j
                            } // ivarl
                        }
// #endif  // end B^T element matrix ------------------------------------
                    } // end loop ivar
                } // end loop i

//------------------    QL    -------------------------------------------------------------------------------------------------------------------

#if FSI_EQUATIONS==1   // only coupled Assemblying Matrix linear-quadratic ---------------------------------
                for (int  ikl=1; ikl<2; ikl++) {  // ikl=0 discontinuous ikl=1 continuous pressure
                    for (int   ivarl=0; ivarl< _nvars[ikl]; ivarl++) {
                        for (int   i=0; i< el_ndof[ikl]; i++) { // +++++++++++
                            // set up row i
                            int  indx=i+el_ndof[2]*_nvars[2];//ivar=0;
                            const double psii_g=_phi_g[ikl][i];
//           if (mode == 1)   FeM(indx)  +=  JxW_g[2]*psii_g*_ub_g[1][0]*KOMP_NS; //    _u_old[DIMENSION]*KOMP_NS;
//             FeM(indx)  =.001;
//           KeM(indx,indx)=1;
                            // linear-quadratic Assemblying Matrix -----------------------------
                            for (int j=0; j<el_ndof[2]; j++)     {
                                const double phij_g= _phi_g[2][j];
//            if (ls_interface>0.4)
//              std::cout<<"I'm here"<< indx<<" "<< j +el_ndof[2];
                                for (int  jvar=0; jvar< _nvars[2]; jvar++) { // linear -quad
                                    // p-equation
                                    KeM(indx,j+jvar*el_ndof[2]) += JxW_g[2]*rho*(
                                                                       psii_g*_dphi_g[2][j+jvar*el_ndof[2]]  // div=0
#ifdef AXISYM
                                                                       +(1-jvar)*psii_g*phij_g/_ub_g[2][0]
#endif
                                                                   );
                                }// jvar
                            }  // j end linear-quad --------------------------------------------
#endif
                            // linear-linear Assemblying Matrix ------------------------------
//           for (int j=0; j<el_ndof[1]; j++)  {
//          const double psij_g=_phi_g[1][i];
//             KeM(i+el_ndof[2]*_nvars[2],j+ndim*el_ndof[2])  += JxW_g[2]*psii_g*psij_g*KOMP_NS;
//           } // end linear liner --------------------------------------------
                        }  // i
                    }// ivarl end linear +++++++++++
                }
#endif  // ------------------------------------------------------------------------------------------------

            } // end of the quadrature point qp-loop

            kappa_mg[0] /=el_ngauss;
            kappa_mg[1] /=el_ngauss;
            mu_m /=el_ngauss;

            // ====================== end volume (element) =======================================



            // ======================================================================
            // ========================= fluid boundary  =================================
            // ======================================================================
            for (int  iside=0; iside< el_sides; iside++) {
                if (el_neigh[iside] == -1) {

                    // local matrix and rhs
//         double IRe_eff=_IRe;

                    // setup boundary element -> connectivity+coordinates
                    for (int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) {
                        int lnode=_mgmesh._GeomEl._surf_top[lbnode+NDOF_FEMB*iside];// local nodes
                        sur_toply[lbnode]=lnode;          // lbnode -> lnode
//           elb_conn[lbnode]=el_conn[lnode];  // connctivity el_conn->elb_conn

                        for (int idim=0; idim<ndim; idim++) { // coordinates
                            double xyz=xx_qnds[idim*NDOF_FEM+lnode];
                            _data_eq[2].ub[idim*NDOF_FEM+lbnode]=xxb_qnds[idim*NDOF_FEMB+lbnode]=xyz;
                        }

                    }


                    // normal
                    _fe[2]->normal_g(xxb_qnds,normal);
                    // pressure flag
//         int flag_p=0;
//         for(int  k=0; k< elb_ndof[1]; k++) { flag_p += _bc_bd[sur_toply[k]+DIMENSION*NDOF_FEM];}

                    for (int ivar=0; ivar< _nvars[2]; ivar++)    {

                        // Dirichlet boundary conditions  ***********************
                        if (_bc_vol[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM]%2 ==0 ) {
//             int bc_s=(int)_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM];     // b cond
                            double det= _fe[2]->JacSur(elb_ngauss-1,xxb_qnds,InvJac[2]);              // jacobian
                            double Ipenalty=det/_dt;                        // Dirichlet bc flag

                            for (int i=0; i< elb_ndof[2]; i++) { // +++++++++++++++++++++++

                                const int  indx_sol=sur_toply[i]+(ivar+_dir)*el_ndof[2];//ivar=0;
                                const int  indx_row=sur_toply[i]+(ivar)*el_ndof[2];//ivar=0;

                                // flag setup (\int_[S bc_normal*(u \dot n)+bc_tg*(u-(u.n)n)+bc_rhs*val] dS)
                                int bc_s=(int)_bc_bd[indx_row]; //int bc_v=(int)_bc_vol[indx_row];
                                int bc_rhs=((bc_s&4)>>2);  // (1??) nonhomogeneous
                                int bc_tg=((bc_s&2)>>1);   // (?1?) tg
                                int bc_normal=(bc_s%2);    // (??1) normal

                                double ndotu2=0.;  // non-homogeneous normal
                                for (int kdim=0; kdim<ndim; kdim++)
                                    ndotu2 +=normal[kdim]*_data_eq[2].ub[ns_idx*NDOF_FEM+sur_toply[i]+kdim*NDOF_FEM];

                                // Assemblying rhs // non-homogeneous ----------------------------
                                if (mode == 1)  {
                                    FeM(indx_row) += bc_rhs*Ipenalty*( // Dirichlet -> non-homogeneous flag (1??)
//                 (1-bc_normal)*(1-bc_tg)*u1_in
                                                         (1-bc_normal)*(1-bc_tg)*_data_eq[2].ub[ns_idx*NDOF_FEM+indx_row]// _u_old[indx_sol]    // (100) single comp bc
                                                         +bc_normal*ndotu2 *normal[ivar+_dir]              // (101)normal Dirichlet bc
                                                         +bc_tg*(_data_eq[2].ub[ns_idx*NDOF_FEM+indx_sol]-ndotu2*normal[ivar+_dir])// (110) tg Dirichlet bc
                                                     );
                                    if ( FeM(indx_row)  !=  FeM(indx_row) ) std::cout<<" found NaN on FeM \n";//check if there are any NaN on FeM
//              std::cout<<" "<<FeM(indx_row);
                                }


                                // Assemblying Matrix ---------------------------------
                                // u.n normale tg -(u.n)n
                                KeM(indx_row,indx_row) += ((1-bc_normal)*(1-bc_tg)+bc_tg)*Ipenalty // (?00) single comp bc
                                                          + Ipenalty*normal[ivar+_dir]*normal[ivar+_dir]*(
                                                              +bc_normal   // (?01) Normal Dirichlet bc
                                                              -bc_tg       // (?10) tg Dirichlet bc
                                                          );

                                for (int   jvar=ivar+_dir+1; jvar< ndim+_dir; jvar++)    {
                                    int      jvar1=jvar%DIMENSION;
#if FSI_EQUATIONS==2
                                    FeM(indx_row) += -1.*Ipenalty*_data_eq[2].ub[(ns_idx+jvar1)*NDOF_FEM+sur_toply[i]]*
#else
                                    if(fabs(normal[jvar1]*normal[ivar+_dir]) >1.e-13 || bc_normal >= 1) {
//               double a= fabs(normal[jvar1]*normal[ivar+_dir]);
//               std::cout << a;
                                    }



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
                        else if (_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] !=0) {
                            double alpha_t=1.;
#ifdef TBK_EQUATIONS
                            double kappa_w=kappa_mg[0]+1.e-20;
                            double yplus=rho*0.547722557505*_y_dist*sqrt(kappa_w)/(_IRe);
                            double u_plus=(yplus< 11.66) ? yplus:2.5*log(9.*yplus);
//              u_plus= yplus;
                            alpha_t=rho*0.547722557505*sqrt(kappa_w)/u_plus;
//          std::cout << " \n  yplus "<<  yplus << " kappa " << kappa_w << " y " << _y_dist<< " alpha_t " << alpha_t;
#endif

                            for (int  qp=0; qp<  elb_ngauss; qp++) { //gaussian integration loop (n_gauss)

                                // quad/linear  [2]=quad [1]=linear------------------------------------
                                det[2]     = _fe[2]->JacSur(qp,xxb_qnds,InvJac[2]);   // local coord _phi_g and jac
                                JxW_g[2]=det[2]*_fe[2]->_weight1[(ndim-1)-1][qp]; // weight
                                _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);   // global coord _phi_g
//   interp_el_sol(_data_eq[2].ub,0,DIMENSION,_phi_g[2],elb_ndof[2],_ub_g[2]);
#ifdef AXISYM   // axisymmetric  (index ->0)
                                interp_el_sol(_data_eq[2].ub,0,DIMENSION,_phi_g[2],elb_ndof[2],_ub_g[2]);
                                JxW_g[2]  *=_ub_g[2][0];
#endif
                                // ***********************************
                                for (int i=0; i< elb_ndof[2]; i++) { // Assemblying NS equation
                                    // set up row i
                                    const double phii_g=_phi_g[2][i];
//                 const int   indx_sol=sur_toply[i]+(ivar+_dir)*el_ndof[2];// volume dof index
                                    const int   indx_row=sur_toply[i]+ivar*el_ndof[2];// volume dof index
                                    // boundary flag
                                    int bc_s=(int)_bc_bd[indx_row];
                                    int bc_v=(int)_bc_vol[indx_row]%2;

                                    double dtJxW_g=JxW_g[2]*((bc_v==0)?0:1)*_dt;
                                    int bc_rhs   =((bc_s&4)>>2); // (1??) -> nonhomogeneous
                                    int bc_tg    =((bc_s&2)>>1); // (?1?) -> tg
                                    int bc_normal=(bc_s%2);      // (??1) -> normal
                                    // Assemblying rhs ----------------------------
                                    if (mode == 1)   {
                                        double omega=0.1;
                                        FeM(indx_row)  += -bc_rhs*dtJxW_g*_phi_g[2][i]*normal[ivar+_dir]*(1500*
//                     sin(2*3.1415/omega*time)
                                                          1
                                                          /*((time<1.)? time:1.)*//_refvalue[DIMENSION]);
                                    }

                                    // Assemblying Matrix ---------------------------------
                                    for (int j=0; j<elb_ndof[2]; j++) { // quad -quad
                                        const double phij_g= _phi_g[2][j];
                                        KeM(indx_row,indx_row) +=dtJxW_g*bc_tg*alpha_t*phij_g*phii_g// (?1?) -> tg
                                                                 + dtJxW_g*alpha_t*normal[ivar+_dir]*normal[ivar+_dir]*phij_g*phii_g*(
                                                                     bc_normal // (??1) -> normal
                                                                     -bc_tg     // (?1?) -> tg
                                                                 );


                                        for (int  jvar=ivar+_dir+1; jvar< ndim+_dir; jvar++)    { // u.n normale tg -(u.n)n
                                            int jvar1=jvar%DIMENSION;
#if FSI_EQUATIONS==2
                                            FeM(indx_row) +=   -_data_eq[2].ub[ns_idx*NDOF_FEM+sur_toply[j]+jvar1*el_ndof[2]]*
#else
                                            KeM(indx_row,sur_toply[j]+jvar1*el_ndof[2]) +=
#endif
                                                               0*dtJxW_g*alpha_t*normal[jvar1]*normal[ivar+_dir]*phij_g*phii_g*(
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

        }//end phase==0
        
        
        

// ==========================================================================
// ==========================================================================
//=======================             _ _     _ =============================
//=======================   ___  ___ | (_) __| |=============================
//=======================  / __|/ _ \| | |/ _` |=============================
//=======================  \__ \ (_) | | | (_| |=============================
//=======================  |___/\___/|_|_|\__,_|=============================
//===========================================================================
// ==========================================================================
// ==========================================================================



        /// c) gaussian integration loop (n_gauss)
        if (phase==1) {
//              if (ls_interface>0.4)
//              std::cout<<"I'm here";
            for (int qp=0; qp<  el_ngauss; qp++) {
                // shape functions at gaussian points --------------------------------------------------------------
                // quadratic continuous
                det[2]      = _fe[2]->Jac(qp,xx_qnds,InvJac[2]);     // Jacobian
                JxW_g[2] =det[2]*_fe[2]->_weight1[ndim-1][qp];       // weight
                _fe[2]->get_phi_gl_g(ndim,qp,_phi_g[2]);
                _fe[2]->get_dphi_gl_g(ndim,qp,InvJac[2],_dphi_g[2]); // global coord deriv
                // discontinuous and continuous linear
                for (int ideg=0; ideg<2; ideg++) if (_nvars[ideg]>0) _fe[ideg]->get_phi_gl_g(ndim,qp,_phi_g[ideg]);  // shape funct



                // linear and quadratic fields ---------------------------------------------------------------------
                // linear fields -> pressure
                // ------------------------------------------------------------------------------------------------------------
// #if (FSI_EQUATIONS%2==0)  // projection and segregated
//       interp_el_sol(_data_eq[1].ub,0,2,_phi_g[1],el_ndof[1],_ub_g[1]);
// #else                   // coupled
//       interp_el_sol(_data_eq[1].ub,0,1,_phi_g[1],el_ndof[1],_ub_g[1]);
//       _ub_g[1][1]=0.;
// #endif
                // -----------------------------------------------------------------------------------------------------------

                // quadratic fields
                interp_el_sol(_data_eq[2].ub,0,_data_eq[2].indx_ub[_data_eq[2].n_eqs],_phi_g[2],el_ndof[2],_ub_g[2]);  // field
                interp_el_gdx(_data_eq[2].ub,ns_idx,DIMENSION,_dphi_g[2],el_ndof[2],vel_gdx); // derivatives


                //calculation deformation gradient
#if DIMENSION==2
                double J= (1+vel_gdx[0+0*DIMENSION]*_dt)*(1+vel_gdx[1+1*DIMENSION]*_dt)-( vel_gdx[1+0*DIMENSION]*_dt*vel_gdx[0+1*DIMENSION]*_dt);
// 1 + vel_gdx[0+0*DIMENSION] + vel_gdx[1+1*DIMENSION];//-( vel_gdx[1+0*DIMENSION]*vel_gdx[0+1*DIMENSION]); //third invariant
// J= 1;
// std::cout<< 1 + vel_gdx[0+0*DIMENSION] + vel_gdx[1+1*DIMENSION]<< "   "<< (1+vel_gdx[0+0*DIMENSION])*(1+vel_gdx[1+1*DIMENSION])-( vel_gdx[1+0*DIMENSION]*vel_gdx[0+1*DIMENSION])<<std::endl;
                double I= pow((vel_gdx[0+0*DIMENSION]*_dt+1),2)+ pow((vel_gdx[0+1*DIMENSION]*_dt),2)+
                          pow((vel_gdx[1+1*DIMENSION]*_dt+1),2)+ pow((vel_gdx[1+0*DIMENSION]*_dt),2); //first invariant of the trasformation

#endif
#if DIMENSION==3
                double J= (
                              (1+vel_gdx[0+0*DIMENSION]*_dt)*(1+vel_gdx[1+1*DIMENSION]*_dt)*(1+vel_gdx[2+2*DIMENSION]*_dt)
                              +  vel_gdx[0+1*DIMENSION]*_dt* vel_gdx[1+2*DIMENSION]*_dt* vel_gdx[2+0*DIMENSION]*_dt
                              +  vel_gdx[0+2*DIMENSION]*_dt* vel_gdx[1+0*DIMENSION]*_dt* vel_gdx[2+1*DIMENSION]*_dt

                          )
                          -(
                              vel_gdx[2+0*DIMENSION]*_dt*(1+vel_gdx[1+1*DIMENSION]*_dt)*vel_gdx[0+2*DIMENSION]*_dt
                              +  vel_gdx[2+1*DIMENSION]*vel_gdx[1+2*DIMENSION]* (1+vel_gdx[0+0*DIMENSION]*_dt)
                              + (1+vel_gdx[2+2*DIMENSION]*_dt) *vel_gdx[1+0*DIMENSION]*_dt* vel_gdx[0+1*DIMENSION]*_dt

                          );
// 1 + vel_gdx[0+0*DIMENSION] + vel_gdx[1+1*DIMENSION];//-( vel_gdx[1+0*DIMENSION]*vel_gdx[0+1*DIMENSION]); //third invariant
// J= 1;
// std::cout<< 1 + vel_gdx[0+0*DIMENSION] + vel_gdx[1+1*DIMENSION]<< "   "<< (1+vel_gdx[0+0*DIMENSION])*(1+vel_gdx[1+1*DIMENSION])-( vel_gdx[1+0*DIMENSION]*vel_gdx[0+1*DIMENSION])<<std::endl;
                double I= pow((vel_gdx[0+0*DIMENSION]*_dt+1),2)+ pow((vel_gdx[0+1*DIMENSION]*_dt),2)+pow((vel_gdx[0+2*DIMENSION]*_dt),2)+
                          pow((vel_gdx[1+1*DIMENSION]*_dt+1),2)+ pow((vel_gdx[1+0*DIMENSION]*_dt),2)+ pow((vel_gdx[1+2*DIMENSION]*_dt),2)+ //first invariant of the trasformation
                          pow((vel_gdx[2+2*DIMENSION]*_dt+1),2)+ pow((vel_gdx[2+0*DIMENSION]*_dt),2)+ pow((vel_gdx[2+1*DIMENSION]*_dt),2) ;
#endif

#ifdef DS_EQUATIONS
                // derivatives of the material displacement(with respect to the material coordinate)
                interp_el_gdx(_data_eq[2].ub,sdx_idx,DIMENSION,_dphi_g[2],el_ndof[2],_ub_dxg);
#endif
#ifdef AXISYM
                JxW_g[2]  *=_ub_g[2][0];
                JxW_g[1]  *=_ub_g[2][0];
#endif
                double IRe_eff=_IRe;


                // velocity[NS_F] -> (quad,_indx_eqs[NS_F+idim]) ---------------------------------------------------
                double mod2_vel=1.e-20;
                _sP=1.e-20;

                for (int idim=0; idim< DIMENSION; idim++) {
                    // velocity field  -> gaussian points <- _ub_g[2][ns_idx+idim]
                    vel_g[idim]=_ub_g[2][ns_idx+idim];
                    mod2_vel +=vel_g[idim]*vel_g[idim];
                    // turbulence production term  -> _sP
                    for (int jdim=0; jdim< DIMENSION; jdim++) {
//        vel_gdx[idim+jdim*DIMENSION]=0.;
//           const double sss=vel_gdx[idim+jdim*DIMENSION]+vel_gdx[jdim+idim*DIMENSION];
//           _sP += sss*sss;
//
//      printf("%d =dim %g \n ",_iproc,vel_gdx[idim+jdim*DIMENSION]);
                    }
                }
//       printf("%d %g \n ",_iproc,_sP);
                mod2_vel =sqrt(mod2_vel);
                // upwind term -> f_upwind
                for (int idim=0; idim< DIMENSION; idim++) {
                    Pe_h[idim]=0.5*mod2_vel*h_eff[idim]/_IRe;             // local  Peclet
                    f_upwind[idim]=UP_WIND_FSI*rho*0.5*(1./tanh(Pe_h[idim])-1./Pe_h[idim])*h_eff[idim]/mod2_vel;
                }
                // -------------------- Temperature[T_F] -> (quad,_indx_eqs[T_F]) -------------------------------------
#ifdef T_EQUATIONS
#ifdef TEMP_CONST
                double Temp_g=_ub_g[_indx_ub[ _indx_eqs[1]]];
                rho *= _mgphys.adensity(Temp_g);
                mu *= _mgphys.aviscosity(Temp_g);
#endif
#endif
                // ----------------- Turbulent viscosity [K_F] -> (quad,_indx_eqs[K_F]) ---------------------------
                _mu_turb=0.;// eff visc turb at g pt
#ifdef TBK_EQUATIONS
                const int k_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[K_F]];  // kappa-equations
                _kappa_g[0]= _ub_g[2][k_idx];
                _kappa_g[1]= _ub_g[2][k_idx+1];
                f_mu(val_tbg);
                kappa_mg[0] +=_kappa_g[0];
                kappa_mg[1] +=_kappa_g[1];
#endif
                IRe_eff = _IRe*(1.+_mu_turb);//  visc_eff  at g pt
                mu_m +=IRe_eff;
                // ---------------------------------------------------------------------------------------------------
//  for (int i=0; i< el_ndof[2]; i++) {

//         for (int i=0; i< el_ndof[2]; i++) {
//         for (int  ivar=0; ivar< _nvars[2]; ivar++)    {
//    for (int j=0; j<el_ndof[2]; j++) {
//          FeM(j)+=2;
//         KeM(j,j)+=1;}}};


                /// d) Assemblying NS equation
                for (int i=0; i< el_ndof[2]; i++) { // +++++++++++
                    // set up row i
                    for (int
                            idim=0; idim<ndim; idim++) {
                        dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
                    }

                    const double phii_g=_phi_g[2][i];
                    // loop u,v,w
                    for (int  ivar=0; ivar< _nvars[2]; ivar++)    {

                        int    indx=i+ivar*el_ndof[2];//ivar=0;
// std::cout<< "\n "<<(_bc_vol[i+ivar*NDOF_FEM]%2)<< " "<<(_bc_vol[i+ivar*NDOF_FEM]);
                        double dtxJxW_g=JxW_g[2]*(_bc_vol[i+ivar*NDOF_FEM]%2)*_dt;//*(( ls_interface> 0.5  )? 0:1);
// _mus=0;
#ifndef COMPRESSIBLE
                        _lambda=0;
                        J=1.;
#endif

                        // Assemblying rhs ----------------------------
                        if (mode == 1)                {
                            FeM(indx)  +=  dtxJxW_g*rho*(
                                               vel_g[ivar+_dir]*phii_g/_dt
//                                +((ivar==0)? -1.:0.)*10*phii_g   // time
//                              + rho*_IFr*_dirg[ivar+_dir]*phii_g  // x-gravity
// #if (FSI_EQUATIONS%2==0)      //       projection segregated
//                              +(_ub_g[1][0] + _ub_g[1][1])*dphiidx_g[2][ivar+_dir]*(1-_nvars[1])  //old pressure
// #endif
                                           );
                        }
                        // Assemblying Matrix quad ---------------------------------
                        // ---------------------------------------------------------
                        for (int j=0; j<el_ndof[2]; j++) {
                            const double phij_g= _phi_g[2][j];

                            // set up
                            double Lap_g=0.,Adv_g=0.,Div_g= 0;
#ifdef  AXISYM
                            Lap_g =2.*(1-ivar)*IRe_eff*
                                   phii_g*phij_g/(_ub_g[2][0]*_ub_g[2][0]);// axysimmetry
#endif
                            for (int kdim=0; kdim<ndim; kdim++) {
                                dphijdx_g[2][kdim] =_dphi_g[2][j+kdim*el_ndof[2]];
                                Adv_g +=vel_g[kdim]*dphijdx_g[2][kdim]*phii_g;
//               Lap_g +=(IRe_eff+f_upwind[kdim]*vel_g[kdim]*vel_g[kdim])*(dphijdx_g[2][kdim]*dphiidx_g[2][kdim]
//               );
                                Lap_g +=(_mus)*(dphijdx_g[2][kdim]*dphiidx_g[2][kdim]
                                               );
                                Div_g +=_ub_dxg[kdim+kdim*ndim];
                            }
                            int idimp1=(ivar+1+_dir)%ndim; // block +1 [2-6-7] ---------------
                            int idimp2=(ivar+2+_dir)%ndim; // block +2 [3-4-8] ------------
// // Building the Cauchy Tensor




#ifdef LINEAR_ELASTICITY
/// Linear tensor
                            sig_1_lin[ivar]=_mus*(2*dphijdx_g[2][ivar]*dphiidx_g[2][ivar]+
                                                  dphijdx_g[2][idimp1]*dphiidx_g[2][idimp1]
#if DIMENSION==3
                                                  + dphijdx_g[2][idimp2]*dphiidx_g[2][idimp2]
#endif
                                                 );
                            sig_1_lin[idimp1]= _mus*(dphijdx_g[2][ivar+_dir]*dphiidx_g[2][idimp1]);
#if DIMENSION==3

                            sig_1_lin[idimp2]= _mus*(dphijdx_g[2][ivar]*dphiidx_g[2][idimp2])
#endif
                                               ;
                            sig_3_lin[ivar]=_lambda*dphijdx_g[2][ivar]*dphiidx_g[2][ivar];
#ifdef AXISYM
                            sig_3_lin[ivar]+=_lambda*(1-ivar)*(dphijdx_g[2][ivar]*phii_g+dphiidx_g[2][ivar]*phij_g+phij_g*phii_g/_ub_g[2][0])/_ub_g[2][0];
#endif
                            sig_3_lin[idimp1]=_lambda*dphijdx_g[2][idimp1]*dphiidx_g[2][ivar];
#ifdef AXISYM
                            sig_3_lin[idimp1]+=_lambda*(1-ivar)*(phij_g*dphiidx_g[2][ivar])/_ub_g[2][0];
                            sig_3_lin[idimp1]+=_lambda*(ivar)*(phii_g*dphijdx_g[2][ivar])/_ub_g[2][0];
#endif
#if DIMENSION==3
                            sig_3_lin[idimp2]= _lambda*dphijdx_g[2][idimp2]*dphiidx_g[2][ivar];

#endif

#endif //LINEAR_ELASTICITY





#ifdef NEO_HOOKIAN
/// Linear tensor
                            double kappa=_lambda+2/3*_mus;
                            if (J<0.0000001) {
                                J=0.0000001;/*std::cout<<"J to little"<<std::endl;*/
                            }
                            sig_1_lin[ivar]=_mus*(1/pow(J,5/3))*(2*dphijdx_g[2][ivar]*dphiidx_g[2][ivar]+
                                                                 dphijdx_g[2][idimp1]*dphiidx_g[2][idimp1]
#if DIMENSION==3
                                                                 + dphijdx_g[2][idimp2]*dphiidx_g[2][idimp2]
#endif
                                                                );
                            sig_1_lin[idimp1]= _mus*(1/pow(J,5/3))*(dphijdx_g[2][ivar+_dir]*dphiidx_g[2][idimp1]);
#if DIMENSION==3

                            sig_1_lin[idimp2]= _mus*(1/pow(J,5/3))*(dphijdx_g[2][ivar]*dphiidx_g[2][idimp2])
#endif
                                               ;
        #ifdef COMPRESSIBLE                                    
                            sig_3_lin[ivar]=0.5*kappa*dphijdx_g[2][ivar]*dphiidx_g[2][ivar];
#ifdef AXISYM
                            sig_3_lin[ivar]+=0.5*kappa*(1-ivar)*(dphijdx_g[2][ivar]*phii_g+dphiidx_g[2][ivar]*phij_g+phij_g*phii_g/_ub_g[2][0])/_ub_g[2][0];
#endif
                            sig_3_lin[idimp1]=0.5*kappa*dphijdx_g[2][idimp1]*dphiidx_g[2][ivar];
#ifdef AXISYM
                            sig_3_lin[idimp1]+=0.5*kappa*(1-ivar)*(phij_g*dphiidx_g[2][ivar])/_ub_g[2][0];
                            sig_3_lin[idimp1]+=0.5*kappa*(ivar)*(phii_g*dphijdx_g[2][ivar])/_ub_g[2][0];
#endif
#if DIMENSION==3
                            sig_3_lin[idimp2]=0.5*kappa* dphijdx_g[2][idimp2]*dphiidx_g[2][ivar];

#endif
#endif
//            sigma_rhs[ivar]+=dphiidx_g[2][ivar]*(2*_mus*J*I+3*pow(J,5/3)*kappa-6*_mus*J)/(6*pow(J,8/3));
                            
                            
                            
                            sigma_rhs[ivar]=0.;

                            sig_1_qua[ivar]=_mus*(1/pow(J,5/3))*(dphijdx_g[2][ivar]*_ub_dxg[ivar+DIMENSION*ivar]*dphiidx_g[2][ivar]+
                                            dphijdx_g[2][ivar]*_ub_dxg[ivar+DIMENSION*idimp1]*dphiidx_g[2][idimp1]);

                            sig_1_qua[idimp1]=_mus*(1/pow(J,5/3))*(dphijdx_g[2][idimp1]*_ub_dxg[idimp1+DIMENSION*idimp1]*dphiidx_g[2][ivar]+
                                              dphijdx_g[2][ivar]*_ub_dxg[idimp1+DIMENSION*idimp1]*dphiidx_g[2][idimp1]);


                            sigma_rhs[idimp1]=0.;
#if DIMENSION==3

                            sig_1_qua[idimp2]=0.;

                            sigma_rhs[idimp2]=0.;
#endif
                            
//                                                      sig_1_qua[idimp2]=_mus*(
//                                                   dphiidx_g[2][idimp2]*dphijdx_g[2][idimp2]*_ub_dxg[ivar+DIMENSION*idimp2]+
//                                                   dphiidx_g[2][idimp2]*dphijdx_g[2][idimp1]*_ub_dxg[ivar+DIMENSION*idimp1]+
//                                                   dphiidx_g[2][idimp2]*dphijdx_g[2][ivar]*_ub_dxg[ivar+DIMENSION*ivar]
//                                               );
// 
//                             sigma_rhs[idimp2]=_mus*(
//                                                   dphiidx_g[2][idimp2]*dphijdx_g[2][idimp2]*vel_gdx[ivar+DIMENSION*idimp2]+
//                                                   dphiidx_g[2][idimp2]*dphijdx_g[2][idimp1]*vel_gdx[ivar+DIMENSION*idimp1]+
//                                                   dphiidx_g[2][idimp2]*dphijdx_g[2][ivar]*vel_gdx[ivar+DIMENSION*ivar]
//                                               );
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
#endif //NEO_HOOKIAN   










#ifdef VENANT_KIRCHH
/// Linear tensor
                            sig_1_lin[ivar]=_mus*(2*dphijdx_g[2][ivar]*dphiidx_g[2][ivar]);
                            sig_3_lin[ivar]=dphijdx_g[2][ivar]*dphiidx_g[2][ivar];
#ifdef AXISYM
                            sig_3_lin[ivar]+=(1-ivar)*(dphijdx_g[2][ivar]*phii_g+dphiidx_g[2][ivar]*phij_g+phij_g*phii_g/_ub_g[2][0])/_ub_g[2][0];
#endif
                            sig_3_lin[idimp1]=dphijdx_g[2][idimp1]*dphiidx_g[2][ivar];
#ifdef AXISYM
                            sig_3_lin[idimp1]+=(1-ivar)*(phij_g*dphiidx_g[2][ivar])/_ub_g[2][0];
                            sig_3_lin[idimp1]+=(ivar)*(phii_g*dphijdx_g[2][ivar])/_ub_g[2][0];
#endif
#if DIMENSION==3
                            sig_3_lin[idimp2]= dphijdx_g[2][idimp2]*dphiidx_g[2][ivar];

#endif

#endif //VENANT_KIRCHH

                            // diagonal blocks [1-5-9]
//             printf("%d =iproc  %g \n ",_iproc,vel_gdx[idimp1+3*ivar]);
                            KeM(indx,j+ivar*el_ndof[2]) +=dtxJxW_g*rho*(1+dumping)*(
                                                              phij_g*phii_g/_dt  // time
                                            + Adv_g            // advection
//                                             + Lap_g            // viscous Laplacian
                                                              +_dt*(sig_1_lin[ivar]
                                                                    +sig_1_qua[ivar])
                                                              +_dt*sig_3_lin[ivar]

                                                          );
                            FeM(indx) +=_data_eq[2].ub[(sdx_idx)*NDOF_FEM+j+ivar*el_ndof[2]]*dtxJxW_g*rho*(
                                            -sig_1_lin[ivar]
                                            -sig_1_qua[ivar]
                                            -sig_3_lin[ivar]
                                            +sigma_rhs[ivar]
                                        );

#if FSI_EQUATIONS==2    // splitting 
                            FeM(indx) += -_data_eq[2].ub[ns_idx*NDOF_FEM+j+idimp1*el_ndof[2]]*//_u_old[j+idimp1*el_ndof[2]]*
#else       // no splitting
                            KeM(indx,j+idimp1*el_ndof[2]) +=
#endif
                                         dtxJxW_g* rho*_dt*(
                                             sig_1_lin[idimp1]
                                             +sig_1_qua[idimp1]
                                             +sig_3_lin[idimp1]
                                         );
                            FeM(indx) += -_data_eq[2].ub[(sdx_idx)*NDOF_FEM+j+idimp1*el_ndof[2]]*
                                         dtxJxW_g* rho*(
                                             sig_1_lin[idimp1]
                                             +sig_1_qua[idimp1]
                                             +sig_3_lin[idimp1]
                                             +sigma_rhs[idimp1]
                                         );
#if DIMENSION==3

#if FSI_EQUATIONS==2    // splitting         
                            FeM(indx)+=  -_data_eq[2].ub[ns_idx*NDOF_FEM+j+idimp2*el_ndof[2]]*//          _u_old[j+idimp2*el_ndof[2]]*
#else        // no splitting
                            KeM(indx,j+idimp2*el_ndof[2])+=
#endif
                                         dtxJxW_g* rho*_dt*(
                                             sig_1_lin[idimp2]
                                             +sig_1_qua[idimp2]
                                             +sig_3_lin[idimp2]
                                         );

                            // ----------------------------------------------------------------------

                            FeM(indx)+=  -_data_eq[2].ub[(sdx_idx)*NDOF_FEM+j+idimp2*el_ndof[2]]*//          _u_old[j+idimp2*el_ndof[2]]*
                                         dtxJxW_g* rho*(
                                             sig_1_lin[idimp2]+sig_1_qua[idimp2] +sig_3_lin[idimp2]  +sigma_rhs[idimp2]
                                         );
#endif
                        } // end A element matrix quad -quad (end loop on j)---------------------------------------------

#ifndef COMPRESSIBLE
                        // ------------------------------------------------------------------
#if FSI_EQUATIONS==1      // B^T element matrix ( p*div(v) )--------------------
                        for (int  ikl=0; ikl<2; ikl++) {  // ikl=0 discontinuous ikl=1 continuous pressure
                            for (int  ivarl=0; ivarl<_nvars[ikl]; ivarl++) {
                                for (int  j=0; j<el_ndof[ikl]; j++) {
                                    double psij_g=_phi_g[ikl][j];
                                    KeM(indx,j+ndim*el_ndof[2]) +=dtxJxW_g*(  // MPascal
                                                                      -psij_g*dphiidx_g[2][ivar]
#ifdef AXISYM
                                                                      -(1-ivar)*psij_g*phii_g/_ub_g[2][0]
#endif
                                                                  );
                                } // j
                            } // ivarl
                        }
#endif  // end B^T element matrix ------------------------------------

#endif
                    } // end loop ivar
                } // end loop i

//------------------    QL    -------------------------------------------------------------------------------------------------------------------
#define KOMP_NS (0.07)
#if FSI_EQUATIONS==1   // only coupled Assemblying Matrix linear ---------------------------------
                for (int  ikl=0; ikl<2; ikl++) {  // ikl=0 discontinuous ikl=1 continuous pressure
                    for (int   ivarl=0; ivarl< _nvars[ikl]; ivarl++) {
                        for (int   i=0; i< el_ndof[ikl]; i++) { // +++++++++++
                            // set up row i
                            int  indx=i+el_ndof[2]*_nvars[2];//ivar=0;
                            const double psii_g=_phi_g[ikl][i];
                            if (mode == 1)   FeM(indx)  +=  JxW_g[2]*psii_g*_ub_g[1][0]*KOMP_NS; //    _u_old[DIMENSION]*KOMP_NS;
//               if (mode == 1)    FeM(indx) += 1*((ls_interface>0.5)? 0:1)*0.;
//            FeM(indx)  =0.;
//          KeM(indx,indx)+=1*((ls_interface> 0.5)? 1:1);
//          if (ls_interface>0.5) interface_element_cout++;
//               KeM(indx,indx) +=((ls_interface>0.5)? 0:1)*1;
#ifdef COMPRESSIBLE
                            if (_bc_vol[i+_nvars[2]*NDOF_FEM]<3.5) KeM(indx,indx)  =10.;
#endif
//             ((_bc_vol[i+_nvars[2]*NDOF_FEM]>4.5)? 0.00001:1);

                            // linear-quadratic Assemblying Matrix -----------------------------
#ifndef COMPRESSIBLE
                            for (int j=0; j<el_ndof[2]; j++)     {
                                const double phij_g= _phi_g[2][j];
                                for (int  jvar=0; jvar< _nvars[2]; jvar++) { // linear -quad
                                    // p-equation
                                    KeM(indx,j+jvar*el_ndof[2]) += JxW_g[2]*rho*(
                                                                       psii_g*_dphi_g[2][j+jvar*el_ndof[2]]  // div=0
#ifdef AXISYM
                                                                       +(1-jvar)*psii_g*phij_g/_ub_g[2][0]
#endif
                                                                   );
                                }// jvar
                            }  // j end linear-quad --------------------------------------------

//               linear-linear Assemblying Matrix ------------------------------
                            for (int j=0; j<el_ndof[1]; j++)  {
                                const double psij_g=_phi_g[1][i];
                                KeM(i+el_ndof[2]*_nvars[2],j+ndim*el_ndof[2])  += JxW_g[2]*psii_g*psij_g*KOMP_NS;
                            } // end linear liner --------------------------------------------
#endif
// //
//
                        }  // i
                    }// ivarl end linear +++++++++++
                }
#endif  // ------------------------------------------------------------------------------------------------
//
            } // end of the quadrature point qp-loop

            kappa_mg[0] /=el_ngauss;
            kappa_mg[1] /=el_ngauss;
            mu_m /=el_ngauss;

            // ====================== end volume (element) =======================================



            // ======================================================================
            // ========================= solid boundary  =================================
            // ======================================================================
            for (int  iside=0; iside< el_sides; iside++) {
                if (el_neigh[iside] == -1) {

                    // local matrix and rhs
//         double IRe_eff=_IRe;

                    // setup boundary element -> connectivity+coordinates
                    for (int  lbnode=0; lbnode<NDOF_FEMB; lbnode++) {
                        int lnode=_mgmesh._GeomEl._surf_top[lbnode+NDOF_FEMB*iside];// local nodes
                        sur_toply[lbnode]=lnode;          // lbnode -> lnode
//           elb_conn[lbnode]=el_conn[lnode];  // connctivity el_conn->elb_conn

                        for (int idim=0; idim<ndim; idim++) { // coordinates
                            double xyz=xx_qnds[idim*NDOF_FEM+lnode];
                            _data_eq[2].ub[idim*NDOF_FEM+lbnode]=xxb_qnds[idim*NDOF_FEMB+lbnode]=xyz;
                        }
                    }

                    // normal
                    _fe[2]->normal_g(xxb_qnds,normal);
                    // pressure flag
//         int flag_p=0;
//         for(int  k=0; k< elb_ndof[1]; k++) { flag_p += _bc_bd[sur_toply[k]+DIMENSION*NDOF_FEM];}
//            std::cout<<" "<<normal[0]<<" "<<normal[1]<<" "<<normal[2]<< " "<< xxb_qnds[4+0*NDOF_FEMB]<< " "<< xxb_qnds[4+1*NDOF_FEMB]<< " "<< xxb_qnds[4+2*NDOF_FEMB]<<"\n";
                    for (int ivar=0; ivar< _nvars[2]; ivar++)    {

                        // Dirichlet boundary conditions  ***********************
                        if (_bc_vol[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM]%2 ==0 ) {

//             int bc_s=(int)_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM];     // b cond
                            double det= _fe[2]->JacSur(elb_ngauss-1,xxb_qnds,InvJac[2]);              // jacobian
                            double Ipenalty=det/_dt;                        // Dirichlet bc flag

                            for (int i=0; i< elb_ndof[2]; i++) { // +++++++++++++++++++++++

                                const int  indx_sol=sur_toply[i]+(ivar+_dir)*el_ndof[2];//ivar=0;
                                const int  indx_row=sur_toply[i]+(ivar)*el_ndof[2];//ivar=0;

                                // flag setup (\int_[S bc_normal*(u \dot n)+bc_tg*(u-(u.n)n)+bc_rhs*val] dS)
                                int bc_s=(int)_bc_bd[indx_row]; //int bc_v=(int)_bc_vol[indx_row];
                                int bc_rhs=((bc_s&4)>>2);  // (1??) nonhomogeneous
                                int bc_tg=((bc_s&2)>>1);   // (?1?) tg
                                int bc_normal=(bc_s%2);    // (??1) normal

                                double ndotu2=0.;  // non-homogeneous normal
                                for (int kdim=0; kdim<ndim; kdim++)
                                    ndotu2 +=normal[kdim]*_data_eq[2].ub[ns_idx*NDOF_FEM+sur_toply[i]+kdim*NDOF_FEM];

                                // Assemblying rhs // non-homogeneous ----------------------------
                                if (mode == 1)  {
                                    FeM(indx_row) += bc_rhs*Ipenalty*( // Dirichlet -> non-homogeneous flag (1??)
                                                         (1-bc_normal)*(1-bc_tg)*_data_eq[2].ub[ns_idx*NDOF_FEM+indx_sol]// _u_old[indx_sol]    // (100) single comp bc
                                                         +bc_normal*ndotu2 *normal[ivar+_dir]              // (101)normal Dirichlet bc
                                                         +bc_tg*(_data_eq[2].ub[ns_idx*NDOF_FEM+indx_sol]-ndotu2*normal[ivar+_dir])// (110) tg Dirichlet bc
                                                     );
                                    if ( FeM(indx_row)  !=  FeM(indx_row) ) std::cout<<" found NaN on FeM \n";//check if there are any NaN on FeM
                                }


                                // Assemblying Matrix ---------------------------------
                                // u.n normale tg -(u.n)n
                                KeM(indx_row,indx_row) += ((1-bc_normal)*(1-bc_tg)+bc_tg)*Ipenalty // (?00) single comp bc
                                                          + Ipenalty*normal[ivar+_dir]*normal[ivar+_dir]*(
                                                              +bc_normal   // (?01) Normal Dirichlet bc
                                                              -bc_tg       // (?10) tg Dirichlet bc
                                                          );

                                for (int   jvar=ivar+_dir+1; jvar< ndim+_dir; jvar++)    {
                                    int      jvar1=jvar%DIMENSION;
#if FSI_EQUATIONS==2
                                    FeM(indx_row) += -1.*Ipenalty*_data_eq[2].ub[(ns_idx+jvar1)*NDOF_FEM+sur_toply[i]]*
#else
                                    if(fabs(normal[jvar1]*normal[ivar+_dir]) >1.e-13 || bc_normal >= 1) {
//               double a= fabs(normal[jvar1]*normal[ivar+_dir]);
//               std::cout << a;
                                    }



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
                        else if (_bc_bd[sur_toply[NDOF_FEMB-1]+ivar*NDOF_FEM] !=0) {
                            double alpha_t=1.;
#ifdef TBK_EQUATIONS
                            double kappa_w=kappa_mg[0]+1.e-20;
                            double yplus=rho*0.547722557505*_y_dist*sqrt(kappa_w)/(_IRe);
                            double u_plus=(yplus< 11.66) ? yplus:2.5*log(9.*yplus);
//              u_plus= yplus;
                            alpha_t=rho*0.547722557505*sqrt(kappa_w)/u_plus;
//          std::cout << " \n  yplus "<<  yplus << " kappa " << kappa_w << " y " << _y_dist<< " alpha_t " << alpha_t;
#endif

                            for (int  qp=0; qp<  elb_ngauss; qp++) { //gaussian integration loop (n_gauss)

                                // quad/linear  [2]=quad [1]=linear------------------------------------
                                det[2]     = _fe[2]->JacSur(qp,xxb_qnds,InvJac[2]);   // local coord _phi_g and jac
                                JxW_g[2]=det[2]*_fe[2]->_weight1[ndim-2][qp]; // weight
                                _fe[2]->get_phi_gl_g(ndim-1,qp,_phi_g[2]);   // global coord _phi_g
                                interp_el_sol(_data_eq[2].ub,0,DIMENSION,_phi_g[2],elb_ndof[2],_ub_g[2]);
#ifdef AXISYM   // axisymmetric  (index ->0)
                                interp_el_sol(_data_eq[2].ub,0,DIMENSION,_phi_g[2],elb_ndof[2],_ub_g[2]);
                                JxW_g[2]  *=_ub_g[2][0];
#endif
                                // ***********************************

                                for (int i=0; i< elb_ndof[2]; i++) { // Assemblying NS equation
                                    // set up row i
                                    const double phii_g=_phi_g[2][i];
//                 const int   indx_sol=sur_toply[i]+(ivar+_dir)*el_ndof[2];// volume dof index
                                    const int   indx_row=sur_toply[i]+ivar*el_ndof[2];// volume dof index
                                    // boundary flag
                                    int bc_s=(int)_bc_bd[indx_row];
                                    int bc_v=(int)_bc_vol[indx_row]%2;
                                    double dtJxW_g=JxW_g[2]*((bc_v==0)?0:1)*_dt;
                                    int bc_rhs   =((bc_s&4)>>2); // (1??) -> nonhomogeneous
                                    int bc_tg    =((bc_s&2)>>1); // (?1?) -> tg
                                    int bc_normal=(bc_s%2);      // (??1) -> normal

                                    // Assemblying rhs ----------------------------

                                    if (mode == 1)   {
//                double corrector=1.;
//                if(xxb_qnds[i+2*NDOF_FEMB]>0.5 && ivar==2) {corrector=-1;}
                                        FeM(indx_row)  += -bc_rhs*dtJxW_g*_phi_g[2][i]*normal[ivar+_dir]*(1/_refvalue[DIMENSION])*(1-ls_interface);//*corrector;
                                    }
                                    // Assemblying Matrix ---------------------------------
                                    for (int j=0; j<elb_ndof[2]; j++) { // quad -quad
                                        const double phij_g= _phi_g[2][j];
                                        KeM(indx_row,indx_row) +=dtJxW_g*bc_tg*alpha_t*phij_g*phii_g// (?1?) -> tg
                                                                 + dtJxW_g*alpha_t*normal[ivar+_dir]*normal[ivar+_dir]*phij_g*phii_g*(
                                                                     bc_normal // (??1) -> normal
                                                                     -bc_tg     // (?1?) -> tg
                                                                 );

                                        for (int  jvar=ivar+_dir+1; jvar< ndim+_dir; jvar++)    { // u.n normale tg -(u.n)n
                                            int jvar1=jvar%DIMENSION;
#if FSI_EQUATIONS==2
                                            FeM(indx_row) +=   -_data_eq[2].ub[ns_idx*NDOF_FEM+sur_toply[j]+jvar1*el_ndof[2]]*
#else
                                            KeM(indx_row,sur_toply[j]+jvar1*el_ndof[2]) +=
#endif
                                                               0*dtJxW_g*alpha_t*normal[jvar1]*normal[ivar+_dir]*phij_g*phii_g*(
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
            // ====================== end solid boundary ================================
            // ==========================================================================

        }//end phase==1

        // ==========================================================================
        // ====================== end solid ========================================
        // ==========================================================================
        
        
        
        
        
        /// e) Add them to the global matrix ------------------------
        A[Level]->add_matrix(KeM,el_dof_indices);
        if (mode == 1) b[Level]->add_vector(FeM,el_dof_indices);


//   if(iel==0) {
//    cout << KeM << endl;
//   }

    } // end of element loop
    // clean ---------------------------------------------------


    el_dof_indices.clear();
    A[Level]->close();
    if (mode == 1) b[Level]->close();
//   A[Level]->print();

    // ----------------------------------------------------------
#ifdef PRINT_INFO
    std::cout << " GenMatRhs(NS): Assembled  " << Level << " (Level) "  << std::endl;
#endif

    return;
}


// =========================================================================================
/// This function controls the assembly and the solution of the NS_equation system:
void MGSolFSI::MGTimeStep(
    const double time,  // time
    const int /*iter*/  // Number of max inter
) {

    std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;

    /// A) Assemblying of the MG tangent matrices and rhs (MGSolFSI::GenMatRhs)

#if PRINT_TIME==1
    std::clock_t start_time=std::clock();
#endif

    GenMatRhs(time,_NoLevels-1,1);
    // Probe the position that will comunicate
    A[_NoLevels-1]->close();

    for(int Level = 0 ; Level < _NoLevels-1; Level++) {
        GenMatRhs(time,Level,0); // matrix
        A[Level]->close();
    }

#if PRINT_TIME==1
    std::clock_t end_time=std::clock();
    std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif

    /// B) Solution of the linear MGsystem (MGSolFSI::MGSolve).
    MGSolve(1.e-6,40,1,12,40,12);

#if PRINT_TIME==1
    end_time=std::clock();
    std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif

//     const double delta_t = _mgutils.get_par("dt");
//     const double a1 = 1./(0.5*delta_t);
//     const double a2 = (1. - 0.5)/0.5;


    /// C) Update of the old solution at the top Level
    //update of the old solution (old velocity)
    (*x_oold[_NoLevels-1])=(*x_old[_NoLevels-1]);
    //update of the solution (velocity)
    x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);

//   //update of the old acceleration
//   (*xdd_oold[_NoLevels-1])=(*xdd_old[_NoLevels-1]);
//
//   //update of the acceleration
//   (*xdd_old[_NoLevels-1]).scale(-a2);
//   (*xdd_old[_NoLevels-1]).add(a1,*x_old[_NoLevels-1]);
//   (*xdd_old[_NoLevels-1]).add(-a1,*x_oold[_NoLevels-1]);

    return;
}




#ifdef TBK_EQUATIONS
// =====================================================
void  MGSolFSI::f_mu(double val[]) {
    // Turbulent viscosity
    if (_kappa_g[0]< 1.e-20) _kappa_g[0]= 1.e-20; // kappa
    if (_kappa_g[1]< 1.e-20) _kappa_g[1]= 1.e-20; // epsilon/omega
    double tau_k=1.;  // turbulent time

#if (TBK_EQUATIONS==0)   //kappa  (Prandtl length)--------------------------
    tau_k=_y_dist/sqrt(_kappa_g[0]);
#endif  // end kappa -------------------------------------------------------

#if (TBK_EQUATIONS/2==1)  // kappa-epsilon (TBK_EQUATIONS=2 or 3) ---------- 
    tau_k= CMU*_kappa_g[0]/_kappa_g[1]; // tau=1/omega=CMU*kappa/epsilon
#endif   // kappa-epsilon --------------------------------------------------

#if (TBK_EQUATIONS/2==2)       // kappa-omega ------------------------------
    tau_k=1/_kappa_g[1]; // tau=1/omega

#endif  // -----------------------------------------------------------------
    if (tau_k> MAX_TAU) tau_k= MAX_TAU;
// turbulent viscosity
    double Re_t= _kappa_g[0]*tau_k/_IRe;
    // limits
    if (Re_t > MU_TOP)  Re_t =MU_TOP;
    if (Re_t < MU_LOW)  Re_t =MU_LOW;

    _mu_turb=Re_t;

//     // Boundary corrections
//    double R_t=Re_t/CMU;                                             // turbulent Reynolds number
//    double R_eps =_y_dist*sqrt(_kappa_g[0]/sqrt(R_t))/_IRe;
//    double f_mu =(1.-exp(-1.*R_eps/(14.)))*(1.-exp(-1.*R_eps/14.));
//    _mu_turb *=f_mu*(1.+5./pow(R_t,0.75)*exp(-1.*R_t*R_t/40000.));    // correction

// #ifdef LOWRE     /// B) Low-Reynolds  correction model
//     const double alpha_star=1.*(0.024+ 0.166666666667*Re_t)/(1.+0.166666666667*Re_t);
//     _mu_turb *=alpha_star;

// #endif


#ifdef SST

    // F1 and F2 coeff calculation
    double F1,F2;
    double F1_first  = 500.*_IRe*tau_k/(_y_dist*_y_dist);
    double F2_first= 2.*sqrt(_kappa_g[0])*tau_k/(BETASTAR*_y_dist);
    if (F1_first > F2_first) F2_first=F1_first;
    F2=tanh(F2_first*F2_first);

    // nu_t calculation
//       double alpha_star        = 1.;//(0.024+ Re_t/6.)/(1.+ Re_t/6.);
    double alpha_starb = sqrt(_sP)*F2*tau_k/0.31;
    if (alpha_starb < 1.) {
        alpha_starb=1.; /*printf("alpha star is %f and and %f and %f\n",alpha_starb, 1./tau_k, F2);*/
    }
    _mu_turb /= alpha_starb;

// _mu_turb = Re_t;
    if (_mu_turb > MU_TOP)  _mu_turb =MU_TOP;
    if (_mu_turb < MU_LOW)  _mu_turb =MU_LOW;
//       printf("mu_turb is %f \n",_mu_turb);

#endif

    return;
}

#endif

// #if (FSI_EQUATIONS%2==0)
// // ===============================================================
// // --------------   PRESSURE EQUATION [P_F] ------------------
// // ===============================================================
// // ==============================================================
// // FSI_EQUATIONS==0 projection solver (u,v,w) ( P in NSP_EQUATIONS)
// // FSI_EQUATIONS==1 coupled    solver (u,v,w,p)
// // FSI_EQUATIONS==2 segregated solver (u,v,w) ( P in NSP_EQUATIONS)
// // ===============================================================
//
//
//
// // ==================================================================
// /// This function constructs the 3d-2D MGSolFSIP class I
// MGSolFSIP::MGSolFSIP(MGEquationsMap& mg_equations_map_in,
//                    std::string eqname_in,
//                    std::string varname_in):
//     MGSolDA(mg_equations_map_in,0,1,eqname_in,varname_in),
//     /// A) reading parameters
//     _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes
//     _dt(_mgutils.get_par("dt")),       // parameter  dt
//     _rhof(_mgphys.get_par("rho0")),    // parameter  density reference
//     _uref(_mgphys.get_par("Uref")),    // parameter  vel reference
//     _lref(_mgphys.get_par("Lref")) {
//   //  ===============================================================
//   /// B) setting class variables
//   _var_names[0]="p";  _refvalue[0]=_rhof*_uref*_uref;    // class variable names
//
//   /// C ) setting solver type
//   for (int  l=0; l<_NoLevels; l++)  _solver[l]->set_solver_type(SOLVER_NSP);
//
//   /// D) setting nondimensional parameters
//   return;
// }//  ================================================================
//
//
//
//
// // ===============================
// // Initial and boundary conditions
// // ===============================
// /// This function generates the initial conditions for the energy equation:
// // =================================================
// void MGSolFSIP::ic_read(double xp[],double u_value[]) {
// // xp[]=(xp,yp) u_value[]=(u,v,p)
//   u_value[0] =0.;
// }
//
// // ========================================================
// /// This function  defines the boundary conditions for the energy equation:
// void MGSolFSIP::bc_read(double xp[],int bc_Neum[],int bc_flag[]) {
//   // =================================================
//   const double Lref = _mgphys.get_par("Lref");
//   double ILref = 1./Lref;
// #if DIMENSION==2
// // xp[]=(xp,yp) bc_flag[T]=0-Dirichlet 1-Neumann
//   //    boundary conditions box
//
//   if (xp[0]< LXB*ILref+BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }  // left
//   if (xp[0]> LXE*ILref-BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }  // right
//   if (xp[1]< LYB*ILref+BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=1;  }  // bottom
//   if (xp[1]> LYE*ILref-BDRY_TOLL) {  bc_flag[0]=0;   bc_Neum[0]=0;  }  // top
// #endif
//
// #if DIMENSION==3
// // =================================================
//   // xp[]=(xp,yp) bc_flag[u,v,p]=0-Dirichlet 1-Neumann
//   //    boundary conditions box
//   if (xp[0] < LXB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
//   if (xp[0] > LXE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
//   if (xp[1] < LYB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
//   if (xp[2] < LZB*ILref+BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
//   if (xp[2] > LZE*ILref-BDRY_TOLL) {    bc_Neum[0]=1;    bc_flag[0]=0;  }
//   if (xp[1] > LYE*ILref-BDRY_TOLL) {    bc_Neum[0]=0;    bc_flag[0]=0;  }
// #endif
//   return;
// } // end boundary conditions ==========================
//
//
// //  ====================================================
// /// This function assembles the matrix and the rhs:
// //  ====================================================
// void  MGSolFSIP::GenMatRhs(
//   const double /* time*/, // time  <-
//   const int  Level,  // Level <-
//   const  int mode    // mode  <- (1=rhs) (0=only matrix)
// ) {  // ===============================================
//
//   ///A) Set up
//   // geometry ---------------------------------------------------------------------------------------
//   const int   ndim = DIMENSION;                      //dimension
//   const int   offset = _mgmesh._NoNodes[_NoLevels-1];// mesh nodes
//   double      xx_qnds[DIMENSION*NDOF_FEM];           // element node coords
//   int         el_conn[NDOF_FEM];                    // element connectivity
//
//   // gauss integration  ------------------------------------------------------------------------------
//   double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];           // Jac, Jac*w Jacobean
//   double dphijdx_g[3][DIMENSION];  double dphiidx_g[3][DIMENSION]; // global derivatives at g point
//   const int  el_ngauss = _fe[2]->_NoGauss1[ndim-1];                //elem gauss points
//
//   // element dofs: costant[0]-linear[1]-quadratic[2]-------------------------------------------------
//   int el_ndof[3];  el_ndof[0]=1;                  // number of el dofs
//   int el_mat_nrows =0;                            // number of mat rows (dofs)
//   for (int ideg=1; ideg<3; ideg++) {
//     el_ndof[ideg]=_fe[ideg]->_NoShape[ndim-1];
//     el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
//   };
//   int el_mat_ncols = el_mat_nrows;                     //square matrix
//   std::vector<int > el_dof_indices(el_mat_ncols);      // element dof vector
//
//   //  fields -> Navier-Stokes  [NS_F]
//   int ns_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[NS_F]];
//
//   // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) ---------------------------
//   A[Level]->zero();                       if (mode ==1) b[Level]->zero();  // global matrix+rhs
//   DenseMatrixM KeM;                       DenseVectorM FeM;                // local  matrix+rhs
//   KeM.resize(el_mat_nrows,el_mat_ncols);  FeM.resize(el_mat_nrows);        // resize  local  matrix+rhs
//
//
//   const int  nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
//   const int  nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
//   for (int iel=0; iel < (nel_e - nel_b); iel++) {
//
//     ///B) Element  Loop over the volume (n_elem)
//     // set to zero matrix and rhs and center
//     KeM.zero();         FeM.zero();
//
//     // geometry and element  fields ------------------------------------
//     // Element Connectivity (el_conn)  and coordinates (xx_qnds)
//     _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);
//
//     // local system fields
//     get_el_dof_bc(Level,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);
//
//     // element nodes coordinates
//     for (int  d=0; d< NDOF_FEM; d++)
//       for (int  idim=0; idim<DIMENSION; idim++)
//         _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];
//
//     // external fields
//     for (int deg=0; deg<3; deg++) {
//       for (int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
//         _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],
//                                              el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub);
//       }
//     }
//
//
//     /// c) gaussian integration loop (n_gauss)
//     for (int  qp=0; qp< el_ngauss; qp++) { // +++++++++++++++++++++++++++++++++
//
//
//       // shape functions at gaussian points -----------------------------------
//       for (int  ideg=1; ideg<3; ideg++) { // linear-quadratic  [1]=linear [2]=quad
//         det[ideg]      = _fe[ideg]->Jac(qp,xx_qnds,InvJac[ideg]);     // Jacobian
//         JxW_g[ideg] =det[ideg]*_fe[ideg]->_weight1[ndim-1][qp];       // weight
//         _fe[ideg]->get_phi_gl_g(ndim,qp,_phi_g[ideg]);               // shape funct
//         _fe[ideg]->get_dphi_gl_g(ndim,qp,InvJac[ideg],_dphi_g[ideg]); // global coord deriv
//       }
//
//       //external fields (quad only)
//       interp_el_gdx(_data_eq[2].ub,ns_idx,DIMENSION,_dphi_g[2],el_ndof[2],_ub_dxg);
//
//       /// d) Local (element) assemblying pressure equation
//       for (int  i=0; i<el_ndof[1]; i++)     {
//         // set up row i if(xp[1]<0.01) {
// //         u_value[1] =4.;
// //         if(xp[0] >LXB*ILref+BDRY_TOLL && xp[0] <LXE*ILref-BDRY_TOLL) {  u_value[0] =.1; }
//         const double phii_g=_phi_g[1][i];
//         double Div_g=0.;
//         for (int  idim=0; idim< ndim; idim++)  {
//           dphiidx_g[1][idim]=_dphi_g[1][i+idim*el_ndof[1]];
//           Div_g +=_ub_dxg[idim+idim*DIMENSION];      //divergence of velocity field
//         }
//         double dtxJxW_g=JxW_g[2]*_bc_vol[i];
//
//         // Rhs Assemblying  ----------------------------
//         if (mode == 1)  FeM(i) += -1.*dtxJxW_g*Div_g*phii_g/_dt;
//
//         // Matrix Assemblying ---------------------------
//         for (int  j=0; j<el_ndof[1]; j++) {
//
//           double Lap=0.;
//           for (int  idim=0; idim< ndim; idim++) {
//             dphijdx_g[1][idim]=_dphi_g[1][j+idim*el_ndof[1]];
//             Lap += dphijdx_g[1][idim]*dphiidx_g[1][idim]; // Laplacian
//           }
//
//           // pressure matrix assembling
//           KeM(i,i) += 1-_bc_vol[i];
//           KeM(i,j) += dtxJxW_g*_bc_vol[j]*Lap;
//         } // ---------------------------------------------
//       }
//     } // end of the quadrature point qp-loop +++++++++++++++++++++++++
//
//     /// e) Global assemblying pressure equation
//     A[Level]->add_matrix(KeM,el_dof_indices);                  // global matrix
//     if (mode == 1)   b[Level]->add_vector(FeM,el_dof_indices); // global rhs
//
//   } // end of element loop
//   // clean and close
//   el_dof_indices.clear();
//   A[Level]->close();    if (mode == 1)  b[Level]->close();
//
// #ifdef PRINT_INFO
//   std::cout<< " Matrix Assembled(P)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
// #endif
//
//   return;
// }
//
// // =========================================================================================
// /// This function controls the assembly and the solution of the P_equation system:
// void MGSolP::MGTimeStep(
//   const double time,  // time
//   const int /*iter*/  // Number of max inter
// ) {
// // =========================================================================================
//
// /// A) Set up the time step
//   std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;
//
//   /// B) Assemblying of the Matrix-Rhs
// #if PRINT_TIME==1
//   std::clock_t start_time=std::clock();
// #endif
//   GenMatRhs(time,_NoLevels-1,1);                                                // matrix and rhs
//   for (int Level = 0 ; Level < _NoLevels-1; Level++)   GenMatRhs(time,Level,0); // matrix
// #if PRINT_TIME==1
//   std::clock_t end_time=std::clock();
//   std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
// #endif
//
//   /// C) Solution of the linear MGsystem (MGSolP::MGSolve).
//   MGSolve(1.e-6,40);
// #if PRINT_TIME==1
//   end_time=std::clock();
//   std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
//             << "s "<< std::endl;
// #endif
//
//   /// D) Update of the old solution at the top Level  (MGSolP::OldSol_update),
//   x[_NoLevels-1]->localize(*x_oold[_NoLevels-1]);  // dp
//   x_old[_NoLevels-1]->add(1.,*x_oold[_NoLevels-1]);// p
//
//   return;
// }// =======================================================================================
//
//
//
// #endif // ENDIF FSI_EQUATIONS==0
#endif  //ENDIF FSI_EQUATIONS
// #endif  // NS_equation is personal

