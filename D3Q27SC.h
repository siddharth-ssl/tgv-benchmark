#ifndef  _LBM_RD3Q27_H_
#define  _LBM_RD3Q27_H_
enum dvRD3Q37{
DV_ZERO_ZERO_ZERO,
DV_P1_ZERO_ZERO,
DV_M1_ZERO_ZERO,
DV_ZERO_P1_ZERO,
DV_ZERO_M1_ZERO,
DV_ZERO_ZERO_P1,
DV_ZERO_ZERO_M1,
DV_P1_P1_ZERO,                                                                                         
DV_M1_M1_ZERO,                                                                                         
DV_M1_P1_ZERO,                                                                                         
DV_P1_M1_ZERO,                                                                                         
DV_P1_ZERO_P1,                                                                                         
DV_M1_ZERO_M1,                                                                                         
DV_M1_ZERO_P1,                                                                                         
DV_P1_ZERO_M1,                                                                                         
DV_ZERO_P1_P1,                                                                                         
DV_ZERO_M1_M1,                                                                                         
DV_ZERO_M1_P1,                                                                                         
DV_ZERO_P1_M1,                                                                                       
DV_P1_P1_P1,                                                                                        
DV_M1_M1_M1,                                                                                        
DV_P1_M1_P1,                                                                                        
DV_M1_P1_M1,                                                                                        
DV_M1_P1_P1,                                                                                        
DV_P1_M1_M1,                                                                                        
DV_M1_M1_P1,                                                                                        
DV_P1_P1_M1                                                                                    
};                                                                                                     
                                                                                                       
int dvfreeSlipY[27]={};  

int dvOpp[27]= {
DV_ZERO_ZERO_ZERO,
DV_M1_ZERO_ZERO,
DV_P1_ZERO_ZERO,
DV_ZERO_M1_ZERO,
DV_ZERO_P1_ZERO,
DV_ZERO_ZERO_M1,
DV_ZERO_ZERO_P1,                                                                                       
DV_M1_M1_ZERO,  
DV_P1_P1_ZERO,                                                                                         
DV_P1_M1_ZERO,                                                                                         
DV_M1_P1_ZERO,                                                                                         
DV_M1_ZERO_M1,                                                                                         
DV_P1_ZERO_P1,                                                                                         
DV_P1_ZERO_M1,                                                                                         
DV_M1_ZERO_P1,                                                                                         
DV_ZERO_M1_M1,                                                                                         
DV_ZERO_P1_P1,                                                                                         
DV_ZERO_P1_M1,                                                                                         
DV_ZERO_M1_P1,                                                                                      
DV_M1_M1_M1,                                                                                         
DV_P1_P1_P1,                                                                                        
DV_M1_P1_M1,                                                                                        
DV_P1_M1_P1,                                                                                        
DV_P1_M1_M1,                                                                                        
DV_M1_P1_P1,                                                                                        
DV_P1_P1_M1,                                                                                    
DV_M1_M1_P1    
};
                                                                                                       
double T0Inv=3.0; 
double theta0 = 1./3.;                                                                                     
double a2 = T0Inv;                                                                                     
double a4 = a2*a2 ;                                                                                    
double b2 = T0Inv;                                                                                     
double d2 = T0Inv;                                                                                     
                                                                                                       
double wSC  =  (b2-1.0)/(a4*b2);                                                                       
double wFCC =  1.0/(2.0*b2*b2*b2);
double wBCC =  1.0/(8.0*d2*d2*d2);

double sum = 6.0*wSC + 12.0*wFCC + 8.0*wBCC;
double wZ   = 1.0 - sum; 
  
double cx[27]={0.0, 1.0, -1.0, 0.0,  0.0,  0.0,  0.0, 1.0, -1.0, -1.0,  1.0, 1.0, -1.0, -1.0,  1.0, 0.0,  0.0,  0.0,  0.0, 1.0, -1.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0}; 
double cy[27]={0.0, 0.0,  0.0, 1.0, -1.0,  0.0,  0.0, 1.0, -1.0,  1.0, -1.0, 0.0,  0.0,  0.0,  0.0, 1.0, -1.0, -1.0,  1.0, 1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0};
double cz[27]={0.0, 0.0,  0.0, 0.0,  0.0,  1.0, -1.0, 0.0,  0.0,  0.0,  0.0, 1.0, -1.0,  1.0, -1.0, 1.0, -1.0,  1.0, -1.0, 1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0};
double wt[27]={wZ , wSC,  wSC, wSC,  wSC,  wSC,  wSC,wFCC, wFCC, wFCC, wFCC, wFCC, wFCC, wFCC, wFCC, wFCC, wFCC,wFCC,wFCC,wBCC, wBCC, wBCC, wBCC, wBCC, wBCC, wBCC, wBCC}; 

#include<D3Q27SC.cpp>
#endif
