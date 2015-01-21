 #ifndef _equationsconf_                                      
 #define _equationsconf_                                      
                                                              
                                                              
 // ==================================================        
 //                   Fluid mechanics                        
 // ==================================================        
                                                              
 //  Navier-Stokes -----------------------------------       
 //  NS projection -> 0 coupled -> 1  splitted -> 2           
 //  -------------------------------------------------        
     #define NS_EQUATIONS  (1)
     #define NS_IS_PERSONAL 
                                                              
                                                              
 //  Turbulence --------------------------------------       
 //  Length   ->0     Spalart-Allmaras -> 1                   
 //  KE  coupled -> 3  splitted -> 2                          
 //  KW  coupled -> 5  splitted -> 4                          
 // --------------------------------------------------        
 //  #define  TBK_EQUATIONS                                   
                                                              
                                                              
 // ==================================================        
 //               ENERGY                                     
 // ==================================================        
                                                              
 //  Temperature equation ---------------------------        
 //  -------------------------------------------------        
     #define T_EQUATIONS  (1)
     #define T_IS_PERSONAL                 

 //  Energy Turbulence   ----------------------------        
 //  KE  coupled -> 3  splitted -> 2                          
 //  -------------------------------------------------        
 //   #define TTBK_EQUATIONS                               
                                                              
                                                              
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
                                                              
 // ================================================          
 // =========== only if you know how to do =========          
 // ================================================          
 #define DA_EQUATIONS          // needed by all the systems   
                                                              
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
