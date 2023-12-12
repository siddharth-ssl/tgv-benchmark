template <int numfield, typename dataType> 
void printDetails(gridBCC3D< numfield,dataType>  &Grid1,dataType Re,dataType convectionTime,dataType beta, dataType Ma, dataType U0)
{
 std::cout<<"Grid Size:       "<<Grid1.nE1<<"x"<<Grid1.nE2<<"x"<<Grid1.nE3<<std::endl;   
 std::cout<<"Ma:              "<<Ma<<std::endl;
 std::cout<<"U0:              "<<U0<<std::endl; 
 std::cout<<"Re:              "<<Re<<std::endl;
 std::cout<<"Convection Time: "<<convectionTime<<std::endl;
 std::cout<<"beta:            "<<beta<<std::endl;
} 

template <int numfield, typename dataType>
void printVtk(gridBCC3D<numfield,dataType> &gridLB, int step)
{
	std::ofstream file ;
	char fileName[250] ;
	int check = std::system("mkdir -p results");
	sprintf(fileName,"results/velocity_%d.vtk",step) ;
	file.open(fileName) ;

	file<<"# vtk DataFile Version 3.0"<<std::endl<<"Velocity"<<std::endl<<"ASCII"<<std::endl<<"DATASET STRUCTURED_GRID"<<std::endl ;
	file<<"DIMENSIONS "<<(gridLB.m1)<<" "<<(gridLB.m2)<<"  "<<(gridLB.m3)<<std::endl ;
	file<<"POINTS "<< gridLB.actualGridSize <<" double"<<std::endl ;

	for (int k=gridLB.nB3; k<= gridLB.nE3;  k++)
	{
	  for (int j=gridLB.nB2; j<= gridLB.nE2;  j++)
	  {
		  for (int i=gridLB.nB1; i<= gridLB.nE1;  i++)
		  {
			  file<<i<<" "<<j<<" "<<k<<std::endl ;
		  }
	  }
  }

	file<<"POINT_DATA "<< gridLB.actualGridSize <<std::endl ;
	file<<"VECTORS"<<" "<<"velocity"<<" "<<"double"<<std::endl ;

  dataType fTemp[numfield], rho, ux, uy, uz, theta;
	for (int k=gridLB.nB3; k<= gridLB.nE3;  k++)
	{
	  for (int j=gridLB.nB2; j<= gridLB.nE2;  j++)
	  {
		  for (int i=gridLB.nB1; i<= gridLB.nE1;  i++)
		  {  
        copyFromGrid(gridLB, fTemp, i, j, k);
        getMoments(fTemp, rho, ux, uy, uz, theta); 
			  file << ux << " " << uy << " " << uz << std::endl ;
		  }
	  }
  }

	file.close();
}


template <int numfield, typename dataType> 
void printEnstrophy(gridBCC3D< numfield,dataType>  &lbGridR27, dataType beta, dataType dt, int t, dataType convectionTime, dataType U0)
{
 dataType rho, ux, uy, uz, theta, pxx, pyy , pzz, pxy, pxz, pyz, sXX, sYY, sZZ, sXY, sYZ, sXZ;
 dataType fEq[27], fTemp[numfield];
 dataType energy = 0.0, rhosum = 0.0, internalEnergy(0.0), kineticEnergy(0.0), thetasum(0.0), enstrophySum(0.0), enstrophy(0.0), enstrophyMax(0.0);
  
 for (int i3=lbGridR27.nB3 ; i3<=lbGridR27.nE3 ;  i3++)
 {  
  for (int i2=lbGridR27.nB2 ; i2<=lbGridR27.nE2;  i2++)
  { 
   for (int i1=lbGridR27.nB1 ; i1<= lbGridR27.nE1;  i1++)
   {
    copyFromGrid(lbGridR27, fTemp, i1, i2, i3);

    getEnstrophyMoments(fTemp, rho, ux, uy, uz, theta, pxx, pyy , pzz, pxy, pxz, pyz, sXX, sYY, sZZ, sXY, sYZ, sXZ, beta, dt);
    rhosum += rho;
    thetasum += theta;
    internalEnergy += 1.5*rho*theta; 
    kineticEnergy += rho*(ux*ux + uy*uy +uz*uz)*0.5;
    enstrophy = sXX*sXX + sYY*sYY + sZZ*sZZ + 2.0*sXY*sXY + 2.0*sXZ*sXZ + 2.0*sYZ*sYZ; 
    enstrophySum += enstrophy; 
   }
  }
 }
 kineticEnergy  /= lbGridR27.nE1*lbGridR27.nE2*lbGridR27.nE3*U0*U0;
 internalEnergy /= lbGridR27.nE1*lbGridR27.nE2*lbGridR27.nE3*U0*U0;
 enstrophySum   /= lbGridR27.nE1*lbGridR27.nE2*lbGridR27.nE3*U0*U0;
 rhosum         /= lbGridR27.nE1*lbGridR27.nE2*lbGridR27.nE3;
 thetasum       /= lbGridR27.nE1*lbGridR27.nE2*lbGridR27.nE3; 
  std::cout<< (dataType)t/convectionTime << std::setw(16)<<rhosum <<std::setw(20)<< kineticEnergy<<std::setw(18)<<enstrophySum  << std::endl ;
}

