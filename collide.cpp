///////////////// 
//   COLLIDE   //
/////////////////
template <int numfield, typename dataType> 
void copyFromGrid(gridBCC3D< numfield,dataType>  &Grid, dataType* fTemp, int i, int j, int k)
{
  for (int dv=0;dv<numfield;dv++)
    fTemp[dv] = Grid.Node(i,j,k,dv);
}

template <int numfield, typename dataType> 
void copyToGrid(gridBCC3D< numfield,dataType>  &Grid, dataType* fTemp, int i, int j, int k)
{
  for (int dv=0;dv<numfield;dv++)
    Grid.Node(i,j,k,dv) = fTemp[dv];
}

template <int numfield, typename dataType> 
void collide(gridBCC3D< numfield,dataType>  &Grid, double beta)
{
 double fEq[27], fTemp[27], feq[27];
 double rho, ux, uy, uz, theta;
 double twoBeta = 2.0*beta ;
 double factor = 1.0-twoBeta;

 for (int i3=Grid.nB3 ; i3<=Grid.nE3 ;  i3++)
 {  
  for (int i2=Grid.nB2 ; i2<=Grid.nE2;  i2++)
  { 
   for (int i1=Grid.nB1 ; i1<= Grid.nE1;  i1++)
   {
    // getMomentsNode(Grid, i1,i2,i3,rho, ux, uy, uz, theta); 
    copyFromGrid(Grid, fTemp, i1, i2, i3);
    getMoments(fTemp, rho, ux, uy, uz, theta); 
    getF(feq,rho,ux,uy,uz,theta);

    for(int dv=0;dv<numfield;dv++)
     fTemp[dv] = fTemp[dv] + twoBeta*(feq[dv]-fTemp[dv]);  
    
    copyToGrid(Grid, fTemp, i1, i2, i3);     
   }
  } 
 }
}