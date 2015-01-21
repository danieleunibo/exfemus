 #ifndef _equationsconf_                                      
 #define _equationsconf_                                      
                                                              
 enum  FIELDS{
    NS_F  =0,     // [0] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
    NSX_F =0,     // [0] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
    NSY_F =1,     // [1] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)
    NSZ_F =2,     // [2] -> Navier-Stokes (quadratic (2),NS_EQUATIONS)    
    P_F   =3,     // [3] -> Pressure (linear (1),NS_EQUATIONS==0 or 2)
    T_F   =4,     // [4] -> Temperature   (quadratic (2),T_EQUATIONS)
    K_F   =5,     // [5] -> Turbulence K  (quadratic (2),TB K_EQUATIONS)
    EW_F  =6,     // [6] -> Turbulence W  (quadratic (2),TB W_EQUATIONS)
    KTT_F =7,     // [7] -> Turbulence K  (quadratic (2),TB K_EQUATIONS)
    EWTT_F =8,    // [8] -> Turbulence W  (quadratic (2),TB W_EQUATIONS)
    SDSX_F =9,    // [9] -> Displacement (quadratic (2), DS_EQUATIONS)
    SDSY_F =10,   // [10]-> Displacement (quadratic (2), DS_EQUATIONS)
    SDSZ_F =11,   // [11]-> Displacement (quadratic (2), DS_EQUATIONS)
    DA_F   =12,   // [12]-> DA solver (quadratic (2), DA_EQUATIONS)
    DA_P   =13,   // [13]-> DA solver (piecewise, DA_EQUATIONS)
    TA_F   =14,   // [14]-> Temp adjoint
    FS_F   =15,   // [15]-> Temp adjoint g or FSI equations (Dimension)
    CO_F   =18    // [16]-> Color function for FSI equations
  };                                                             
 // ==================================================        
 //                   Fluid mechanics                        
 // ==================================================        
                                                              
 //  Navier-Stokes -----------------------------------       
 //  NS projection -> 0 coupled -> 1  splitted -> 2           
 //  -------------------------------------------------        
  
#define  NS_EQUATIONS (1)   
#define  NSA_EQUATIONS (1)   
                                                              
 //  Turbulence --------------------------------------       
 //  Length   ->0     Spalart-Allmaras -> 1                   
 //  KE  coupled -> 3  splitted -> 2                          
 //  KW  coupled -> 5  splitted -> 4                          
 // --------------------------------------------------        
#define  TBK_EQUATIONS  (0)   
#define  TBKA_EQUATIONS  (0)  
                                                              
                                                              
 // ==================================================        
 //               ENERGY                                     
 // ==================================================        
                                                              
 //  Temperature equation ---------------------------        
 //  -------------------------------------------------        
//     #define T_EQUATIONS  (1)
                                                              
                                                              
 //  Energy Turbulence   ----------------------------        
 //  KE  coupled -> 3  splitted -> 2    
 //  KW  coupled -> 5  splitted -> 4
 //  -------------------------------------------------        
//    #define TTBK_EQUATIONS    (5)                            
                                                              
                                                              
 // ==================================================        
 //          STRUCTURAL MECHANICS                            
 // ==================================================        
 //  SM projection -> 0 coupled -> 1  splitted -> 2          
 //  -------------------------------------------------        
 //  #define SM_EQUATIONS                                 
                                                              
                                                              
 // ==================================================        
 //               NEUTRONICS                                
 // ==================================================        
 //  NEUTRON DIFFUSION                                       
 //  -------------------------------------------------        
                                                              
 //  #define NEU_EQUATIONS                                 
                                                              
                                                              
 // ==================================================        
 //                ELECTROSTATICS                               
 // ==================================================        
 //  POTENTIAL EQUATION                                       
 //  -------------------------------------------------        
                                                              
 //  #define V_EQUATIONS                                 
                                                              
                                                              
 // ==================================================        
 //             FLUID STRUCTURE  FSI                         
 // ==================================================        
                                                              
 //  FSI equations                                           
 //  FSI projection -> 0 coupled -> 1  splitted -> 2          
 //  -------------------------------------------------        
 //  #define FSI_EQUATIONS                                 
                                                              
 // ==================================================        
 //                   MHD                         
 // ==================================================        
                                                              
 // FSI equations                                            
 // FSI projection -> 0 coupled -> 1  splitted -> 2          
 //  -------------------------------------------------        
 // #define MHD_EQUATIONS                                 
                                                              
 // ==================================================        
 //                 TWO_PHASE                            
 // ==================================================        
 //  --------------------------------------                   
 // #define TWOPH_EQUATIONS                                 
                         
                         
//  #define DA1_EQUATIONS      // solver DA3D, test system                       
 // ================================================          
 // =========== only if you know how to do =========          
 // ================================================          
//  #define DA_EQUATIONS         // solver DA3D, test system     
                                                              
 #ifdef FSI_EQUATIONS  
 /// Displacement                                             
 #define  DS_EQUATIONS                                         
 #if (FSI_EQUATIONS%2==0)       // projection method          
  #define FSIP_EQUATIONS (1)     // need separated P equations
 #endif                                                       
 #endif                                                       
                                                              
 #ifdef SM_EQUATION                                           
 #if (SM_EQUATIONS%2==0)       // projection method           
 #define SMP_EQUATIONS (1)     // need separated P equations    
 #endif                                                       
 #endif                                                       
 // ================================================          
                                                              
                                                              
 #endif  // end file _equationsconf_                          
