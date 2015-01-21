#include "ReactData.h"

#include <cmath>

#include "MGMesh.h"
#include "MGSystem.h"
#include "MGUtils.h"

const int ReactData::mat_zone[8][8]={
  {INNER, INNER, INNER, INNER, INNER, INTER, OUTER, DUMMY}, // A1_1 -A1_7
  {INNER, INNER, INNER, INNER, INTER, INTER, OUTER, OUTER}, // A2_1 -A2_8
  {INNER, INNER, INNER, INTER, CTRLR, INTER, OUTER, DUMMY}, // A3_1 -A3_7
  {INNER, INNER, INNER, INNER, INTER, INTER, OUTER, DUMMY}, // A4_1 -A4_7
  {INNER, INTER, INTER, INTER, OUTER, OUTER, DUMMY, DUMMY}, // A5_1 -A5_6
  {INNER, INTER, CTRLR, INTER, OUTER, OUTER, DUMMY, DUMMY}, // A6_1 -A6_6
  {INTER, OUTER, OUTER, OUTER, DUMMY, DUMMY, DUMMY, DUMMY}, // A7_1 -A7_4
  {OUTER, OUTER, OUTER, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY}  // A8_1 -A8_3
};



// ---------------------------------
// Power factor
// ---------------------------------

// //              _                 _  __                           
//     _ __   ___ | |_   _   _ _ __ (_)/ _| ___  _ __ _ __ ___   ___ 
// // | '_ \ / _ \| __| | | | | '_ \| | |_ / _ \| '__| '_ ` _ \ / _ \
// // | |_) | (_) | |_  | |_| | | | | |  _| (_) | |  | | | | | |  __/
// // | .__/ \___/ \__|  \__,_|_| |_|_|_|  \___/|_|  |_| |_| |_|\___|
// // |_|          
// 
// const double ReactData::mat_pf[8][8]={
//   { 1., 1.,1., 1., 1., 1., 1., 0.   }, // A1_1 -A1_7
//   { 1., 1., 1., 1., 1., 1., 1., 1.}, // A2_1 -A2_8
//   { 1., 1., 1., 1., CTL_R, 1., 1., 0.   }, // A3_1 -A3_7
//   { 1., 1., 1., 1., 1., 1., 1., 0.   }, // A4_1 -A4_7
//   { 1., 1., 1., 1., 1., 1., 0.,    0.   }, // A5_1 -A5_6
//   { 1., 1., CTL_R, 1., 1., 1., 0.,    0.   }, // A6_1 -A6_6
//   { 1., 1., 1., 1., 0.,    0.,    0.,    0.   }, // A7_1 -A7_4
//   { 1., 1.,1., 0.,    0.,    0.,    0.,    0.   }  // A8_1 -A8_3
// };
// 
// 
// 
// const double ReactData::axpf[10][3]={
// 
// {8.48697E-01,	8.48697E-01,	8.48697E-01},
// {9.63034E-01,	9.63034E-01,	9.63034E-01},
// {1.06399E-0,	1.06399E-0,	1.06399E-0},
// {1.12959E-0,	1.12959E-0,	1.12959E-0},
// {1.16319E-0,	1.16319E-0,	1.16319E-0},
// {1.15623E-0,	1.15623E-0,	1.15623E-0},
// {1.10313E-0,	1.10313E-0,	1.10313E-0},
// {1.00872E-0,	1.00872E-0,	1.00872E-0},
// {8.55509E-01,	8.55509E-01,	8.55509E-01},
// {7.07905E-01,	7.07905E-01,	7.07905E-01}
// };
// 
// //              _                 _  __                           
// //  _ __   ___ | |_   _   _ _ __ (_)/ _| ___  _ __ _ __ ___   ___ 
// // | '_ \ / _ \| __| | | | | '_ \| | |_ / _ \| '__| '_ ` _ \ / _ \
// // | |_) | (_) | |_  | |_| | | | | |  _| (_) | |  | | | | | |  __/
// // | .__/ \___/ \__|  \__,_|_| |_|_|_|  \___/|_|  |_| |_| |_|\___|
// // |_|          
// 





const double ReactData::mat_pf[8][8]={
  { 0.941, 0.962, 0.989, 1.017, 1.045, 1.178, 1.097, 0.   }, // A1_1 -A1_7
  { 0.949, 0.958, 0.978, 0.996, 1.114, 1.123, 1.193, 0.857}, // A2_1 -A2_8
  { 0.963, 0.975, 0.992, 1.087, CTL_R, 1.011, 0.925, 0.   }, // A3_1 -A3_7
  { 0.967, 0.973, 0.989, 1.008, 1.108, 1.028, 0.931, 0.   }, // A4_1 -A4_7
  { 0.961, 1.060, 1.078, 1.137, 1.198, 0.943, 0.,    0.   }, // A5_1 -A5_6
  { 0.954, 1.029, CTL_R, 0.990, 1.068, 0.850, 0.,    0.   }, // A6_1 -A6_6
  { 1.019, 1.066, 0.966, 0.842, 0.,    0.,    0.,    0.   }, // A7_1 -A7_4
  { 0.878, 0.853, 0.757, 0.,    0.,    0.,    0.,    0.   }  // A8_1 -A8_3
};



const double ReactData::axpf[10][3]={
{8.60089E-01,	8.48697E-01,	8.33685E-01},
{9.32998E-01,	9.63034E-01,	9.52704E-01},
{1.03749E-0,	1.06399E-0,	1.06801E-0},
{1.10010E-0,	1.12959E-0,	1.14484E-0},
{1.14410E-0,	1.16319E-0,	1.17922E-0},
{1.13892E-0,	1.15623E-0,	1.16983E-0},
{1.09049E-0,	1.10313E-0,	1.10469E-0},
{1.01844E-0,	1.00872E-0,	1.00793E-0},
{9.05869E-01,	8.55509E-01,	8.63532E-01},
{7.71503E-01,	7.07905E-01,	6.75547E-01}
};


// ----------------------------------
// axial pressure loss
// ----------------------------------

#ifdef GRIDS_LOSS_OFF

    const double ReactData::axlpf[10+4][3]={
    // Inlet
    {0.5,   0.5,   0.5  },
    {0.,    0.,    0.   },
    // Core
    {0.,    0.,    0.   },
    {0.,    0.,    0.   },
    {0.,    0.,    0.   },
    {0.,    0.,    0.   },
    {0.,    0.,    0.   },
    {0.,    0.,    0.   },
    {0.,    0.,    0.   },
    {0.,    0.,    0.   },
    {0.,    0.,    0.   },
    {0.,    0.,    0.   },
    // Mixing zone
    {1.0,   1.0,   1.0  },
    {0.,    0.,    0.   },
    };

#endif


#ifdef GRIDS_LOSS_ON

    const double ReactData::axlpf[10+4][3]={  // 
    // // Inlet
    {0.5,   0.5,   0.5  }, // inlet abrupt area change
    {0.,    0.,    0.   },
    // // Core
    {0.,    0.,    0.   },
    {0.52,  0.52,  0.52 }, // grid spacer
    {0.,    0.,    0.   },
    {0.52,  0.52,  0.52 }, // grid spacer
    {0.,    0.,    0.   },
    {0.52,  0.52,  0.52 }, // grid spacer
    {0.,    0.,    0.   },
    {0.52,  0.52,  0.52 }, // grid spacer
    {0.,    0.,    0.   },
    {0.52,  0.52,  0.52 }, // grid spacer
    // // Mixing zone
    {1.0,   1.0,   1.0  }, // outlet abrupt area change
    {0.,    0.,    0.   },
    };
		
#endif

    
const double ReactData::axllf[10+4]={  // mesh elements height in the core
// Inlet                               // (to calculate accidentals pressure loss
0.095,                                 // like distributed ones!)
0.095,                                 // BEWARE: YOU HAVE TO SET ALSO THE FUNCTION LOSSF!
// Core
0.09,
0.09, 
0.09,  
0.09,  
0.09,  
0.09,  
0.09,  
0.09, 
0.09,  
0.09,  
// Mixing zone
0.095, 
0.095,  
};
      

// ======================================================
//  This function returns the planar fuel zone
int ReactData::zone(
  double xm, 
  double ym
) {// ==================================================

    int i_ass = (int)((ym+HR)/(2*HR));
    int j_ass = (int)((xm+HR*((i_ass%2)))/(2*HR));
    int zone_f= mat_zone[i_ass][j_ass];
    if (zone_f == DUMMY) std::cerr<< "warn dummy zone \n" << std::endl;
    return zone_f;
}// ==================================================

// ======================================================
//    This function returns the planar fuel power peak
double ReactData::asby_power(
  double xm, 
  double ym
) {// ==================================================
    double peak_f;
    int i_ass = (int)((ym+HR)/(2*HR));
    int j_ass = (int)((xm+HR*(i_ass%2))/(2*HR));
    peak_f =mat_pf[i_ass][j_ass];
    return peak_f;
}// ==================================================

// ======================================================
 //    This function returns the axial  fuel power peak
// --------------------------------
double ReactData::axial_pf(
  double x,
  double y, 
  double z
) {// ==================================================

    if (z >HIN && z< HOUT) {
        int izone= zone(x,y);
        int k=(int)(NZC*(z-HIN)/(HOUT-HIN));
        if (izone > 2) return 1.;
        else return axpf[k][izone];//*0.8316;// /1.2025;//1.2025 e' il valore della bessel j0 qst divisione va fatta solo per il calore uniforme;
    }
    else return 0.;
}// ===================================================

// ======================================================
//    This function returns the presssure  loss
//    over a plane
double ReactData::lossf(
  double x,
  double y, 
  double z
) {// ===================================================

    // int izone= zone(x,y);
    int k = 0;
//     if (z < HIN) k =(int)((2*z)/HIN);
    if (z < HIN) k =(z>=(C_IN) && z<=(C_IN+0.095) ? 0:1);
    else if (z < HOUT+0.1) k = 2+(int)(NZC*(z-HIN)/(HOUT-HIN));
    else k = 12+(int)(2*(z-HOUT)/(2-HOUT));

    return axlpf[k][0];
}// // ===================================================

// ======================================================
//    This function returns the  Axial  pressure  loss
// ----------------------------------
double ReactData::axial_lf(
  double z
) {// ===================================================
    int k = 0;
    if (z < HIN) k =(int)((2*z)/HIN);
    else if (z < HOUT) k = 2+(int)(NZC*(z-HIN)/(HOUT-HIN));
    else k = 12+(int)(2*(z-HOUT)/(2-HOUT));

    return axllf[k];

} // ====================================================


void  ReactData::datagen(
    const MGMesh & mesh_in,  // mesh class ->
    const MGSystem & phys_in,  // system class ->
    const MGUtils & _mgutils
) {  // ===============================================

    /// Set up
    // geometry ---------------------------------------------------------------------------------------
    const int  ndim = DIMENSION;            //dimension
    int        el_conn[NDOF_FEM];           // element connectivity
    double     xx_qnds[DIMENSION*NDOF_FEM]; // element node coords
    double  xm[3];                          // cell midpoint

    for (int Level=0; Level <mesh_in._NoLevels; Level++) {
        /// b) Element  Loop over the volume (n_elem)

        for (int iproc=0; iproc <mesh_in._n_subdom ; iproc++) {
            const int nel_e = mesh_in._off_el[0][Level+mesh_in._NoLevels*iproc+1]; // start element
            const int nel_b = mesh_in._off_el[0][Level+mesh_in._NoLevels*iproc];   // stop element
            for (int iel=0; iel < (nel_e - nel_b); iel++) {

                // Element Connectivity (el_conn)  and coordinates (xx_qnds)
                mesh_in.get_el_nod_conn(0,Level,iel,el_conn,xx_qnds,iproc);

                for (int idim=0; idim<DIMENSION; idim++) {
                    double sum=0.;
                    for (int d=0; d< NDOF_FEM; d++) sum += xx_qnds[idim*NDOF_FEM+d];// element nodes xxg (DIM)
                    xm[idim]=sum/NDOF_FEM;
                }
                
	     double peak_f = 0.;
	     double loss_f = 0.;
	     
	     double PIN = HIN; double POUT = HOUT;
	     
#ifdef IN_OUT_LOSS
		PIN = C_IN; POUT = C_OUT;
#endif
	     
                if (xm[2]>=HIN &&  xm[2]<=HOUT) {
		  
		    std::string msh_name=_mgutils.get_file("MESHNAME");
		    
		    if (msh_name=="box"){
		  
			double pig = 3.14159265359;
			double mod = (xm[0]*xm[0])+(xm[1]*xm[1]);
// 			peak_f = (1+((-0.4530790*pow(mod,1.5))+(0.3217150*mod)+(0.5218949*sqrt(mod))+(-0.3)))*(cos((pig*(xm[2]-2.7))/1.5))
// 				*(0.9/0.797123458); // simple power factor distribution for the box mesh
				
    		        peak_f = 1.;    // to setting uniform power factor distribution (for code validation)
				
			phys_in._sys_data[0][iel+nel_b]=peak_f;
		    
		    }
		    
		    if ((msh_name=="elsy_core")||(msh_name=="elsy")){
		      
			peak_f = asby_power(xm[0],xm[1])*axial_pf(xm[0],xm[1],xm[2]); // actual ELSY power factor distribution
			phys_in._sys_data[0][iel+nel_b]=peak_f;
		      
		    }
		    
                }
                else {
                    phys_in._sys_data[0][iel+nel_b]=0.;
                }
		
		if (xm[2]>=PIN &&  xm[2]<=POUT) {
		  
	            loss_f = lossf(xm[0],xm[1],xm[2])/axial_lf(xm[2]);
                    phys_in._sys_data[1][iel+nel_b]=loss_f;
                }
                else {
                    phys_in._sys_data[1][iel+nel_b]=0.;
                }

            } // end of element loop
        }  // proc loop
#ifdef PRINT_INFO
    std::cout<< " Fuel distribution  for  Level "<< Level <<"\n";
#endif  
      
    } // level loop


    return;
}