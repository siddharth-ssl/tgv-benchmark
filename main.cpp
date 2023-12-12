/*********************************************************************************************
 *   Copyright (c) <2015>, <Santosh Ansumali@JNCASR>                                         *
 *   All rights reserved.                                                                    *
 *   Redistribution and use in source and binary forms, with or without modification, are    *
 *   permitted provided that the following conditions are met:                               *
 *                                                                                           *
 *    1. Redistributions of source code must retain the above copyright notice, this list of *
 *       conditions and the following disclaimer.                                            *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list *
 *       of conditions and the following disclaimer in the documentation and/or other        *
 *       materials provided with the distribution.                                           *
 *    3. Neither the name of the <JNCASR> nor the names of its contributors may be used to   *
 *       endorse or promote products derived from this software without specific prior       *
 *       written permission.                                                                 *
 *                                                                                           *
 *       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND     *
 *       ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED       *
 *       WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  *
 *       IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,    *
 *       INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,      *
 *       BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,       *
 *       DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF     *
 *       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE     *
 *       OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED   *
 *       OF THE POSSIBILITY OF SUCH DAMAGE.                                                  *
 *                                                                                           *
 *       Suggestions:          ansumali@jncasr.ac.in                                         *
 *       Bugs:                 ansumali@jncasr.ac.in                                         *
 *                                                                                           *
 *********************************************************************************************/
#include<math.h>
#include<stdio.h>
#include<iostream>
#include <fstream>
#include<iomanip>
#include "grid3D.h"
#include "gridBC.cpp"
#include "D3Q27SC.h"
#include "advect.cpp"
#include "collide.cpp"
#include "init.cpp"
#include "print.cpp"
#include "tgv_analytical.cpp"

// n1, n2, n3, Ma, Re, iterations, print_freq
typedef double myReal;
int main(int argc, char** argv)
{
 //grid point on node and cell in each direction                         
  int n1 = std::atoi(argv[1]);                                                            
  int n2 = std::atoi(argv[2]);                                                            
  int n3 = std::atoi(argv[3]);
                                                                         
  myReal cs2 = 1.0/3.0;  
  myReal Ma = std::atof(argv[4]);
  myReal U0 = Ma * std::sqrt(cs2);                              
  myReal Re = std::atof(argv[5]);    

  myReal len = 2.0*M_PI;                                                 
  myReal charLen = 1.0;                                                   
  myReal kinVisc = U0*charLen/Re ;                                        
  myReal tau = kinVisc/cs2 ;                                              
  myReal dx1 = 1.0;                                                    
  myReal dt = dx1 ;                                                       
  myReal Kn = sqrt(3.0)*U0/(Re);                                       
  myReal ratioTau = 1.0;                            
  tau = tau/dt;                                                           
  myReal beta = 1.0/(2.0*tau + 1.0);                                     
  int hNode =1;
  gridBCC3D<27,myReal> lb_grid(n1, n2, n3, hNode,hNode,hNode);

  myReal convectionTime = n1/(2.*M_PI*U0);

  printDetails(lb_grid, Re, convectionTime, beta, Ma, U0);  
 
  InitializeTGV(lb_grid, U0);

  int iterations = (int)(std::atof(argv[6])*convectionTime);//(5.0);  
  int print_freq = (int)(std::atof(argv[7])*convectionTime);//(int)(0.1*convectionTime);  

  int count=0;
  for(int t=0;t<iterations;t++)
  {  
    if(t%print_freq == 0) {
      l2norm(lb_grid, U0, kinVisc, t);
      printVtk(lb_grid, count);
      count++;
    }

    collide(lb_grid,beta);
    periodicAll(lb_grid);
    advection(lb_grid);

  }
  return 0;
}
