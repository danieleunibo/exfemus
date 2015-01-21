#include "Equations_conf.h"
// ============================================
#ifdef DS_EQUATIONS // 3D-2D Energy equation
// ============================================
// class files --------------------------------------------------
#include "MGSclass_conf.h"        // Navier-Stokes class conf file
#include "MGSolverDS.h"

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


// ======================================================
/// This function constructs the 3d-2D MGSolDSX_Y_Z class
MGSolDS::MGSolDS(MGEquationsSystem& mg_equations_map_in,
                 const int vars_in[],
                 std::string eqname_in,
                 std::string varname_in):
    MGSolDA(mg_equations_map_in,vars_in,eqname_in,varname_in),
    // mesh params ------------
    _offset(_mgmesh._NoNodes[_NoLevels-1]), // mesh nodes
    // phys params ------------
    _dt(_mgutils.get_par("dt")),       // parameter  dt
    _rhof(mg_equations_map_in.get_par("rho0")),    // parameter  density reference
    _uref(mg_equations_map_in.get_par("Uref")),    // parameter  vel reference
    _lref(mg_equations_map_in.get_par("Lref")),    // parameter  length reference
    _Tref(mg_equations_map_in.get_par("Tref")),    // parameter  temperature reference
    _rhos(mg_equations_map_in.get_par("rhos"))     // parameter solid density
{   //  =================================================
    // class variable names
    _dir=0;

    if(!varname_in.compare("dx")) _dir=0;
    if(!varname_in.compare("dy")) _dir=1;
    if(!varname_in.compare("dz")) _dir=2;
#if DIMENSION ==2
    if (_dir==2) std::cout<<"Too many Dimension!!\n";
#endif
    _var_names[0]=varname_in;
    _refvalue[0]=_lref;

    for (int l=0; l<NDOF_FEM; l++) {
        _bc_vol[l]=-1;
        _bc_bd[l]=-1;
    }

    // class solver type (SOLVERT  in MGSTconf.h)
//   for (int l=0;l<_NoLevels;l++) _solver[l]->set_solver_type(CGM);
    for (int l=0; l<_NoLevels; l++) _solver[l]->set_solver_type(GMRESM);

    return;
}

//  ====================================================
/// This function assembles the matrix and the rhs:
//  ====================================================
void  MGSolDS::MoveMesh(
    const int Level  // Level <-
) {  // ===============================================

    //update dx_old ========================================
    /// E) mesh update
    const int n_nodes=_mgmesh._NoNodes[Level];
    const int n_elem=_mgmesh._NoElements[0][Level];
    const int offsetp=_dir*n_nodes;
    const int  offset = _mgmesh._NoNodes[Level];                     // mesh nodes
    int        el_conn[NDOF_FEM];//, elb_conn[NDOF_FEMB];                     // element connectivity
    double     xx_qnds[DIMENSION*NDOF_FEM];
    int el_mat_nrows =0;
    int el_ndof[3];
    el_ndof[0]=0;
    const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
    int        sur_toply[NDOF_FEMB]; // boundary topology

    for (int ideg=1; ideg<3; ideg++) {                              //     ...
        el_ndof[ideg]=((_nvars[ideg]>0)?    _fe[ideg]->_NoShape[DIMENSION-1]:0);                    //   computing
        el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
    };
//  _mgmesh._node_map[0][0];
    int el_mat_ncols = el_mat_nrows;                     //square matrixbc_Neum[0]=0;
    std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector


    int Level_c=0;//_NoLevels-1;
    DenseVectorM LoD;//Local Displacement
    LoD.resize(el_ndof[2]);
    for (int proc=0; proc<_mgmesh._n_subdom; proc++) {
        int ndof_lev=0;
        for(int pr=0; pr <proc; pr++) {
            int delta =_mgmesh._off_el[0][pr*_NoLevels+Level_c+1]-_mgmesh._off_el[0][pr*_NoLevels+Level_c];
            ndof_lev +=delta;
        }


        //find the coarse nodes
        const int coarse_nel_e =_mgmesh._off_el[0][Level_c+_NoLevels*proc+1]; // start element
        const int coarse_nel_b =_mgmesh._off_el[0][Level_c+_NoLevels*proc];   // stop element
        int        coarse_node[NDOF_FEM];//, elb_conn[NDOF_FEMB];                     // element connectivity
        std::cout;
        //coarse elemnt loop
        for(int iel=0; iel < (coarse_nel_e - coarse_nel_b); iel++) {
            LoD.zero();
            //getting connectivity and coordinates
            _mgmesh.get_el_nod_conn(0,Level_c,iel,el_conn,xx_qnds,proc);
            //getting boundary contions
            get_el_dof_bc(Level_c,iel+ndof_lev,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

            int count=0;
            int written[NDOF_FEM]; //needed to have fine displacement on the solid
            for(int inode=0; inode<el_ndof[2]; inode++)     {
                coarse_node[inode]=  _mgmesh._node_map[Level_c][el_dof_indices[inode]];
                int glob_indx=_mgmesh._node_map[Level_c][el_dof_indices[inode]];
                if (_bc_vol[inode]>1.5) {
                    LoD(inode) =(*x_old[Level])(glob_indx)-(*x_oold[Level])(glob_indx);
                    written[count]=inode;
                    count++;
                }//end solid
                else { //liquid
                    if (inode < NDOF_P ) LoD(inode)= (*x_old[Level])(glob_indx)-(*x_oold[Level])(glob_indx);
                }

            }//end loop node
#if DIMENSION==2
            for(int iside=0; iside< el_sides; iside++)  {
                int check=0;
                for(int idof=0; idof<NDOF_FEMB; idof++)  sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
                check=0;
                for (int i=0; i<count; i++) if (written[i]==sur_toply[NDOF_FEMB-1]) check =1;
                if (check==0) {
                    LoD(sur_toply[NDOF_FEMB-1]) = 0.5*(LoD(sur_toply[NDOF_FEMB-3]))
                                                  + 0.5*(LoD(sur_toply[NDOF_FEMB-2]));
                }

            }//end iside
#endif
#if DIMENSION==3
            for(int iside=0; iside< el_sides; iside++)  {
//                 int check=0;
                for(int idof=0; idof<NDOF_FEMB; idof++)  sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[5]) check =1;
//                 if (check==0) {
                   if (_bc_vol[sur_toply[4]]<1.5)  LoD(sur_toply[4]) = 0.5*(LoD(sur_toply[0]))
                                        + 0.5*(LoD(sur_toply[1]));
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[6]) check =1;
//                 if (check==0) {
                      if (_bc_vol[sur_toply[5]]<1.5) LoD(sur_toply[5]) = 0.5*(LoD(sur_toply[1]))
                                        + 0.5*(LoD(sur_toply[2]));
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[7]) check =1;
//                 if (check==0) {
                     if (_bc_vol[sur_toply[6]]<1.5)  LoD(sur_toply[6]) = 0.5*(LoD(sur_toply[2]))
                                        + 0.5*(LoD(sur_toply[3]));
		      if (_bc_vol[sur_toply[7]]<1.5)  LoD(sur_toply[7]) = 0.5*(LoD(sur_toply[3]))
                                        + 0.5*(LoD(sur_toply[0]));
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[8]) check =1;
//                 if (check==0) {
                    if (_bc_vol[sur_toply[8]]<1.5)   LoD(sur_toply[8]) = 0.25*(LoD(sur_toply[0]))
                                        +0.25*(LoD(sur_toply[1]))
                                        + 0.25*(LoD(sur_toply[2]))
                                        + 0.25*(LoD(sur_toply[3]));
//                 }

            }//end iside
#endif
            int check=0;
            for (int i=0; i<count; i++) if (written[i]==NDOF_FEM-1) check =1; //last node
            if (check==0) {
                LoD(NDOF_FEM-1)=0;
                for (int jnode=0; jnode<NDOF_P; jnode++)  LoD(NDOF_FEM-1)+=LoD(jnode)/NDOF_P;
            }

            for (int inode=0; inode < el_ndof[2] ; inode++) disp[Level]->set(_mgmesh._node_map[Level_c][el_dof_indices[inode]],LoD(inode)); //setting the global displacement
        }//end loop iel
        el_dof_indices.clear();
    }//end iproc
// 	    end coarse elemnt loop
//              end find the coarse nodes
// //finer level
    for (int ilev=Level_c; ilev<_NoLevels; ilev++) {
        for (int proc=0; proc<_mgmesh._n_subdom; proc++) {

            int ndof_lev=0;
            for(int pr=0; pr <proc; pr++) {
                int delta =_mgmesh._off_el[0][pr*_NoLevels+ilev+1]-_mgmesh._off_el[0][pr*_NoLevels+ilev];
                ndof_lev +=delta;
            }
            const int nel_e =_mgmesh._off_el[0][ilev+_NoLevels*proc+1]; // start element
            const int nel_b =_mgmesh._off_el[0][ilev+_NoLevels*proc];   // stop element
            for(int iel=0; iel < (nel_e - nel_b); iel++) {
                LoD.zero();
                //getting connectivity and coordinates
                _mgmesh.get_el_nod_conn(0,ilev,iel,el_conn,xx_qnds,proc);
                //getting boundary contions
                get_el_dof_bc(ilev,iel+ndof_lev,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);

                int count=0;
                int written[NDOF_FEM]; //needed to have fine displacement on the solid
                int countg=0;
                for(int inode=0; inode<el_ndof[2]; inode++)     {
                    int dof_idx=_mgmesh._node_map[ilev][el_dof_indices[inode]];
                    if (_bc_vol[inode]>1.5) {
                        written[count]=inode;
                        count++;
//                      if (inode > NDOF_P )
                        {
                            LoD(inode) =(*x_old[Level])(dof_idx)-(*x_oold[Level])(dof_idx);

                        }

                    }//end solid
                    else { //liquid
                        if (inode < NDOF_P ) LoD(inode)=(*disp[Level])(dof_idx);
                    }

                }//end loop node
     #if DIMENSION==2
                for(int iside=0; iside< el_sides; iside++)  {
                    int check=0;
                    for(int idof=0; idof<NDOF_FEMB; idof++)  sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
                    check=0;
                    for (int i=0; i<count; i++) if (written[i]==sur_toply[NDOF_FEMB-1]) check =1;
                    if (check==0) {
                        LoD(sur_toply[NDOF_FEMB-1]) = 0.5*(LoD(sur_toply[NDOF_FEMB-3]))
                                                      + 0.5*(LoD(sur_toply[NDOF_FEMB-2]));

                    }

                }
#endif
                
                #if DIMENSION==3
            for(int iside=0; iside< el_sides; iside++)  {
//                 int check=0;
                for(int idof=0; idof<NDOF_FEMB; idof++)  sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[5]) check =1;
//                 if (check==0) {
                   if (_bc_vol[sur_toply[4]]<1.5)  LoD(sur_toply[4]) = 0.5*(LoD(sur_toply[0]))
                                        + 0.5*(LoD(sur_toply[1]));
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[6]) check =1;
//                 if (check==0) {
                      if (_bc_vol[sur_toply[5]]<1.5) LoD(sur_toply[5]) = 0.5*(LoD(sur_toply[1]))
                                        + 0.5*(LoD(sur_toply[2]));
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[7]) check =1;
//                 if (check==0) {
                     if (_bc_vol[sur_toply[6]]<1.5)  LoD(sur_toply[6]) = 0.5*(LoD(sur_toply[2]))
                                        + 0.5*(LoD(sur_toply[3]));
		      if (_bc_vol[sur_toply[7]]<1.5)  LoD(sur_toply[7]) = 0.5*(LoD(sur_toply[3]))
                                        + 0.5*(LoD(sur_toply[0]));
//                 }
//                 check=0;
//                 for (int i=0; i<count; i++) if (written[i]==sur_toply[8]) check =1;
//                 if (check==0) {
                    if (_bc_vol[sur_toply[8]]<1.5)   LoD(sur_toply[8]) = 0.25*(LoD(sur_toply[0]))
                                        +0.25*(LoD(sur_toply[1]))
                                        + 0.25*(LoD(sur_toply[2]))
                                        + 0.25*(LoD(sur_toply[3]));
//                 }

            }//end iside
#endif

                int check=0;
                for (int i=0; i<count; i++) if (written[i]==NDOF_FEM-1) check =1; //last node
                if (check==0) {
                    for (int jnode=0; jnode<NDOF_P; jnode++)  LoD(NDOF_FEM-1)+=LoD(jnode)/NDOF_P;

                }
                int check2=0;
                for (int inode=0; inode < el_ndof[2] ; inode++) {
                    disp[Level]->set(_mgmesh._node_map[ilev][el_dof_indices[inode]],LoD(inode));
                }
            }//end loop iel

            el_dof_indices.clear();
        }//end iproc
    }//end ilev
    disp[Level]->close();
    return;
}
//  ====================================================
/// This function assembles the matrix and the rhs:
//  ====================================================
void  MGSolDS::GenMatRhs(
    const double /* time*/, // time  <-
    const int Level,  // Level <-
    const  int mode    // mode  <- (1=rhs) (0=only matrix)
) {  // ===============================================

    /// Set up
    // geometry -----
    const int  ndim = DIMENSION;                                           //dimension
    double     xx_qnds[DIMENSION*NDOF_FEM], xxb_qnds[DIMENSION*NDOF_FEMB]; // element node coords
    int        el_conn[NDOF_FEM];//, elb_conn[NDOF_FEMB];                     // element connectivity
    int        el_neigh[NDOF_FEM];                                         // element connectivity
    const int  offset = _mgmesh._NoNodes[_NoLevels-1];                     // mesh nodes
    const int  el_sides= _mgmesh._GeomEl._n_sides[0];                      // element nodes
    int        sur_toply[NDOF_FEMB]; // boundary topology
//   double     vel_g[DIMENSION];  for(int idim=0;idim<DIMENSION;idim++) vel_g[idim] =0.;
//   double  val_tbg[3];
    int fl_int;
    int phase;
//   int interface;
//   double rho;
    double u_old;
    // gauss integration  ------------------
    const int el_ngauss = _fe[2]->_NoGauss1[ndim-1];         //elem gauss points
    const int elb_ngauss = _fe[2]->_NoGauss1[DIMENSION-2];   //elem gauss points
    double det[3],JxW_g[3],InvJac[3][DIMENSION*DIMENSION];    // Jac, Jac*w Jacobean
    double dphijdx_g[3][DIMENSION];
    double dphiidx_g[3][DIMENSION];                           // global derivatives at g point

    //number of constant[0]-linear[1]-quadratic[2] element dof
    int el_ndof[3];
    el_ndof[0]=0;
    int elb_ndof[3];
    elb_ndof[0]=0;
    int el_mat_nrows =0;

    for (int ideg=1; ideg<3; ideg++) {                              //     ...
        el_ndof[ideg]=((_nvars[ideg]>0)?    _fe[ideg]->_NoShape[ndim-1]:0);                    //   computing
        elb_ndof[ideg]=((_nvars[ideg]>0)?_fe[ideg]->_NoShape[ndim-2]:0);                  //     ...
        el_mat_nrows +=_nvars[ideg]*el_ndof[ideg];
    };

    int el_mat_ncols = el_mat_nrows;                     //square matrixbc_Neum[0]=0;
    std::vector<int> el_dof_indices(el_mat_ncols);      // element dof vector

    const int fsi_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[FS_F]];       // Fluid-Structure Interaction equation
    const int disp_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[SDSX_F]];    // Displacement equation


    // element matrix and rhs  (mode 0= matrix only or mode 1=matrix +rhs) -----------
    A[Level]->zero();
    if(mode ==1) b[Level]->zero();
    DenseMatrixM KeM;
    DenseVectorM FeM;
    KeM.resize(el_mat_nrows,el_mat_ncols);
    FeM.resize(el_mat_nrows);
    int ndof_lev=0;
    for (int pr=0; pr <_mgmesh._iproc; pr++) {
        int delta =_mgmesh._off_el[0][pr*_NoLevels+Level+1]-_mgmesh._off_el[0][pr*_NoLevels+Level];
        ndof_lev +=delta;
    }
    /// b) Element  Loop over the volume (n_elem)
    const int nel_e = _mgmesh._off_el[0][Level+_NoLevels*_iproc+1]; // start element
    const int nel_b = _mgmesh._off_el[0][Level+_NoLevels*_iproc];   // stop element
    for(int iel=0; iel < (nel_e - nel_b); iel++) {

        // set to zero matrix and rhs and center
        KeM.zero();
        FeM.zero();

        // geometry and element  fields ------------------------------------
        // Element Connectivity (el_conn)  and coordinates (xx_qnds)
        _mgmesh.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds);
        _mgmesh.get_el_neighbor(el_sides,0,Level,iel,el_neigh);

        // set element-nodes variables  bc (bc_q_dofs)
//     get_el_dof_bc(Level,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);
        get_el_dof_bc(Level,iel+ndof_lev,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);
//         get_el_dof_bc(Level,iel,el_ndof,el_conn,offset,el_dof_indices,_bc_vol,_bc_bd);
        // external fields
        for(int d=0; d< NDOF_FEM; d++) {   // element nodes xxg (DIM)
            for(int idim=0; idim<DIMENSION; idim++) {
                _data_eq[2].ub[idim*NDOF_FEM+d]=xx_qnds[idim*NDOF_FEM+d];
            }
        }

        for(int deg=2; deg<3; deg++) {
            for(int eq=0; eq<_data_eq[deg].n_eqs; eq++) {
                _data_eq[deg].mg_eqs[eq]->get_el_sol(0,_data_eq[deg].indx_ub[eq+1]-_data_eq[deg].indx_ub[eq],el_ndof[deg],el_conn,offset,_data_eq[deg].indx_ub[eq],_data_eq[deg].ub);
            }
        }

        //  external cell properties -------------------------------------
        phase=(_bc_vol[NDOF_FEM-1]<1.5)?0:1;
//     for(int i=1;i<el_sides+1;i++)  if(_bc_vol[NDOF_FEM-i] > 4.5 ) interface=1.;
//     rho = 1.;

// liquid phase
        if(phase == 0) {
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


            // ======================================================================
            // Volume =============================================================
            // ======================================================================

            /// c) gaussian integration loop (n_gauss)
            // --------------------------------------------
            for(int qp=0; qp< el_ngauss; qp++) {

                // shape functions at gaussian points -----------------------------------
                det[2]      = _fe[2]->Jac(qp,xx_qnds,InvJac[2]);     // Jacobian
                JxW_g[2] =det[2]*_fe[2]->_weight1[ndim-1][qp];       // weight
                _fe[2]->get_phi_gl_g(ndim,qp,_phi_g[2]);               // shape funct
                _fe[2]->get_dphi_gl_g(ndim,qp,InvJac[2],_dphi_g[2]); // global coord deriv


                /// d) Local (element) assemblying energy equation
                // *********************** *******************************
                for(int i=0; i<el_ndof[2]; i++)     {

                    // set up row i
                    const double phii_g=_phi_g[2][i];

                    for(int idim=0; idim< ndim; idim++)  dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];

// 	fl_int=1;  if(interface==1) fl_int=0;
                    fl_int= ((_bc_vol[i]>3.5)? 0:1);
//        if(interface==1) std::cout<< _bc_vol[i]<<" \n";
                    double dtxJxW_g=JxW_g[2]*(_bc_vol[i]%2)*fl_int*_dt;
                    u_old=_data_eq[2].ub[(disp_idx+_dir)*NDOF_FEM+i];
                    FeM(i)=dtxJxW_g*(u_old/_dt*phii_g)*0;
                    // Matrix Assemblying ---------------------------
                    for(int j=0; j<el_ndof[2]; j++) {



#ifdef ALE_MOVING

// // // // 	  ALE MOVING MESH
//          double phij_g= _phi_g[2][j];
                        double Lap=0.;
                        // set up test function
                        const double phij_g=_phi_g[2][j];
                        for(int idim=0; idim< ndim; idim++) {
                            dphijdx_g[2][idim]=_dphi_g[2][j+idim*el_ndof[2]];
                            Lap +=dphijdx_g[2][idim]*dphiidx_g[2][idim]; // Laplacian
                        }
                        // energy-equation
                        KeM(i,j) +=dtxJxW_g*(phij_g*phii_g/_dt*0+
                                             100*Lap)
                                   ;






// // // 	 END ALE MOVING MESH
#endif
#ifdef ELASTIC_MOVING
                        double  _mus=0;
                        double _lambda=0;
                        double rho=1.;
                        int sdx_idx=_data_eq[2].indx_ub[_data_eq[2].tab_eqs[SDSX_F]];      // NS equation [NS_F]]
                        double sig_1_lin[DIMENSION], sig_3_lin[DIMENSION]; //double sig_1_cons[DIMENSION];



                        // set up row i
                        for (int idim=0; idim<ndim; idim++) {
                            dphiidx_g[2][idim]=_dphi_g[2][i+idim*el_ndof[2]];
                        }

                        const double phii_g=_phi_g[2][i];
                        // loop u,v,w
                        for (int  ivar=0; ivar< _nvars[2]; ivar++)    {

                            int    indx=i+ivar*el_ndof[2];//ivar=0;
// std::cout<< "\n "<<(_bc_vol[i+ivar*NDOF_FEM]%2)<< " "<<(_bc_vol[i+ivar*NDOF_FEM]);
                            double dtxJxW_g=JxW_g[2]*(_bc_vol[i+ivar*NDOF_FEM]%2)*_dt;//*(( ls_interface> 0.5  )? 0:1);

//             if (mode == 1)                {
//               FeM(indx)  +=  dtxJxW_g*rho*(
//                                vel_g[ivar+_dir]*phii_g/_dt // time
//                                +phii_g*0
//                              );
//             }

                            // Assemblying Matrix quad ---------------------------------
                            // ---------------------------------------------------------
                            for (int j=0; j<el_ndof[2]; j++) {
                                const double phij_g= _phi_g[2][j];

                                int idimp1=(ivar+1+_dir)%ndim; // block +1 [2-6-7] ---------------
                                int idimp2=(ivar+2+_dir)%ndim; // block +2 [3-4-8] ------------
// // Building the Cauchy Tensor



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

                                // diagonal blocks [1-5-9]
//             printf("%d =iproc  %g \n ",_iproc,vel_gdx[idimp1+3*ivar]);
                                KeM(indx,j+ivar*el_ndof[2]) +=dtxJxW_g*rho*(
//                                               phij_g*phii_g/_dt  // time
                                                                  +_dt*sig_1_lin[ivar]
                                                                  +_dt*_lambda*sig_3_lin[ivar]

                                                              );
                                FeM(indx) +=_data_eq[2].ub[(sdx_idx)*NDOF_FEM+j+ivar*el_ndof[2]]*dtxJxW_g*rho*(
                                                -sig_1_lin[ivar]
                                                -_lambda*sig_3_lin[ivar]
                                            );



                                FeM(indx) += -_data_eq[2].ub[(sdx_idx)*NDOF_FEM+j+idimp1*el_ndof[2]]*
                                             dtxJxW_g* rho*(
                                                 sig_1_lin[idimp1]
                                                 +_lambda*sig_3_lin[idimp1]
                                             );
#if DIMENSION==3




                                FeM(indx)+=  -_data_eq[2].ub[(sdx_idx)*NDOF_FEM+j+idimp2*el_ndof[2]]*//          _u_old[j+idimp2*el_ndof[2]]*
                                             dtxJxW_g* rho*(
                                                 sig_1_lin[idimp2]+sig_1_qua[idimp2] +_lambda*sig_3_lin[idimp2]  +sigma_rhs[idimp2]
                                             );
#endif
                            } // end A element matrix quad -quad (end loop on j)---------------------------------------------


                        } // end loop ivar




#endif  // ELASTIC_MOVING





                    }
                } // ----------------------------------------
            } // end of the quadrature point qp-loop ***********************


            // ======================================================================
            // ====================== fluid boundary ======================================
            // ======================================================================

            for(int iside=0; iside< el_sides; iside++)  {
                if(el_neigh[iside] == -1) {

                    for(int idof=0; idof<NDOF_FEMB; idof++) {
                        sur_toply[idof]=_mgmesh._GeomEl._surf_top[idof+NDOF_FEMB*iside];// local nodes
                        int idofb=sur_toply[idof];
//           elb_conn[idof]=el_conn[idofb];                 // connctivity vector
                        for(int idim=0; idim<DIMENSION; idim++) {
                            xxb_qnds[idim*NDOF_FEMB+idof]=xx_qnds[idim*NDOF_FEM+idofb]; // coordinates
                            _data_eq[2].ub[idim*NDOF_FEM+idof]=xxb_qnds[idim*NDOF_FEMB+idof];
                        }
                    }

                    // Dirichlet boundary conditions  ***********************************
                    if((_bc_vol[sur_toply[NDOF_FEMB-1]]%2) ==0) {

                        //[NDOF_FEMB-1] is the midpoint of a quadaratic FEM element (HEX9 (2D) or HEX27 (3D))
                        int bc_s=(int)_bc_bd[sur_toply[NDOF_FEMB-1]];     // b cond
                        double det= _fe[2]->JacSur(elb_ngauss-1,xxb_qnds,InvJac[2]);// jacobian
                        double Ipenalty=det/_dt;                               // Dirichlet bc flag
                        // local boundary loop   ---------------------------------------
                        for(int lb_node=0; lb_node< elb_ndof[2]; lb_node++) {

                            int lv_node= sur_toply[lb_node]; // local vol index
                            // flag setup (\int bc_var*T+bc_var*val)
//            int  bc_val = (int)((bc_s&2)>>1);  // (1?)non homogeneous
                            int  bc_var = (int)(bc_s%2);       // (?1) variable
                            // Assemblying  Matrix & rhs
                            fl_int= ((_bc_vol[lv_node]>3.5)? 0:1);
// 	    if (interface==1 ) std::cout<< "\n "<< _bc_vol[lv_node];
                            if(mode == 1) FeM(lv_node) += fl_int*bc_var*Ipenalty*( _data_eq[2].ub[fsi_idx*NDOF_FEM+lv_node]*_dt +  _data_eq[2].ub[disp_idx*NDOF_FEM+lv_node] );
                            KeM(lv_node,lv_node) += fl_int*Ipenalty;  //  Dirichlet bc
                        }// lb_node -end  local boundary loop -------------------------
                    } // end if Dirichlet  boundary conditions
                    // **********************************************************************

                } //end if side
            } // ======================  end for fluid boundary =======================================
        }


// ==========================================================================
// ==========================================================================
// ======================             _ _     _ =============================
// ======================   ___  ___ | (_) __| |=============================
// ======================  / __|/ _ \| | |/ _` |=============================
// ======================  \__ \ (_) | | | (_| |=============================
// ======================  |___/\___/|_|_|\__,_|=============================
// ==========================================================================
// ==========================================================================
// ==========================================================================

        else {

            /// d) Local (element) assemblying energy equation
            // ************************************************
            for(int i=0; i<el_ndof[2]; i++)     {
                // set up row i
                // Rhs Assemblying  ----------------------------
                FeM(i) = _data_eq[2].ub[(fsi_idx+_dir)*NDOF_FEM+i]*_dt + _data_eq[2].ub[(disp_idx+_dir)*NDOF_FEM+i];

                // Matrix Assemblying ---------------------------
                KeM(i,i)=1.;

            }  // end of the i-loop
        }


        /// e) Global assemblying energy equation
        A[Level]->add_matrix(KeM,el_dof_indices);                  // global matrix
        if(mode == 1)   b[Level]->add_vector(FeM,el_dof_indices);  // global rhs
//   cout << endl;
//   cout << KeM << endl;
//   cout << endl;
//   cout << endl;



    } // end of element loop
    // clean
    el_dof_indices.clear();
#ifdef PRINT_INFO
    std::cout<< " Matrix Assembled(T)  for  Level "<< Level << " dofs " << A[Level]->n() <<"\n";
#endif



    return;
}



// =========================================================================================
/// This function controls the assembly and the solution of the T_equation system:
void MGSolDS::MGTimeStep(const double time, const int /*iter*/) { // ------------------------------------

// #ifdef AXISYM
// if(_dir==1) return;
// #endif

    std::cout  << std::endl << "  " << _eqname.c_str() << " solution "  << std::endl;
    /// B) Assemblying of the rhs (top level in surface and volume with MGSolNS::GenRhs,GenRhsB),
#if PRINT_TIME==1
    std::clock_t start_time=std::clock();
#endif
    GenMatRhs(time,_NoLevels-1,1);
    A[_NoLevels-1]->close();
    /// C) Assemblying of the  MGmatrices (MGSolNS::GenMatrix),
    for(int Level = 0 ; Level < _NoLevels-1; Level++) {
        GenMatRhs(time,Level,0); // matrix
        A[Level]->close();
    }
#if PRINT_TIME==1
    std::clock_t end_time=std::clock();
    std::cout << "  Assembly time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC << " s "<< std::endl;
#endif

    /// E) Solution of the linear MGsystem (MGSolNS::MGSolve).
    MGSolve(1.e-6,40);
#if PRINT_TIME==1
    end_time=std::clock();
    std::cout << " Assembly+solution time -----> ="<< double(end_time- start_time) / CLOCKS_PER_SEC
              << "s "<< std::endl;
#endif
    /// A) Update of the old solution at the top Level  (MGSolNS::OldSol_update),
    //update x_oold
    x_oold[_NoLevels-1]->zero();
    x_old[_NoLevels-1]->localize(*x_oold[_NoLevels-1]);
//   //update x_old
    x[_NoLevels-1]->localize(*x_old[_NoLevels-1]);

    const int flag_moving_mesh = _mgutils.get_par("moving_mesh");
    const int n_nodes=_mgmesh._NoNodes[_NoLevels-1];
    const int n_elem=_mgmesh._NoElements[0][_NoLevels-1];
    const int offsetp=_dir*n_nodes;
    if(flag_moving_mesh) {
#ifdef FINE_MOVEMENT
        /// E) mesh update
        const int n_nodes=_mgmesh._NoNodes[_NoLevels-1];
        const int offsetp=_dir*n_nodes;
        for (int inode=0; inode<n_nodes; inode++) {
            double disp=(*x_old[_NoLevels-1])(inode)-(*x_oold[_NoLevels-1])(inode);
            _mgmesh._xyz[inode+offsetp] += disp;
            _mgmesh._dxdydz[inode+offsetp] = disp;
        }
#endif
#ifdef COARSE_MOVEMENT
        MoveMesh(_NoLevels-1);
// disp[_NoLevels-1]->zero();
        for (int inode=0; inode<n_nodes; inode++) {
            _mgmesh._xyz[inode+offsetp] +=  (*disp[_NoLevels-1])(inode);
            _mgmesh._dxdydz[inode+offsetp]= (*disp[_NoLevels-1])(inode);
        }
//   MoveMesh(_NoLevels-1);
#endif
    }
    // ==============================================================
    return;
}














#endif


