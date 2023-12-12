/////////////////////////
// Initial Conditions  //
/////////////////////////
template <int numfield, typename dataType>
void InitializeTGV(gridBCC3D<numfield,dataType> &gridLB, dataType U0)
{
  dataType ke = 0;
  dataType x, y, rho , u1, u2, u3, theta, feq[numfield];
  dataType kL= 2.*M_PI/(gridLB.m1);

  for (int k=gridLB.nB3; k<=gridLB.nE3; k++) {
    for(int j=gridLB.nB2; j<=gridLB.nE2; j++) {
      for(int i=gridLB.nB1; i<=gridLB.nE1; i++) {
			  x = (i);
			  y = (j); 

	      rho= 1.0 - (0.75*U0*U0*(cos(2*kL*x) + cos(2*kL*y)));
        u1 = U0*cos(kL*x)*sin(kL*y);
	      u2 = -U0*sin(kL*x)*cos(kL*y);
        u3 = 0.;
        theta = 1./3.;

        getF(feq, rho, u1, u2, u3, theta);
      
        copyToGrid(gridLB, feq, i, j, k);  
      }
    }
  }
}

template <int numfield, typename dataType> 
void InitializeKida(gridBCC3D< numfield,dataType>  &lbGridR27, dataType U0)
{
 dataType feq[27];
 dataType theta(1.0/3.0), x, y, z;
 int n3 = lbGridR27.nE3;
 int n2 = lbGridR27.nE2;
 int n1 = lbGridR27.nE1; 
 dataType rho(1.0), u1, u2, u3;

 for (int k=lbGridR27.nB3; k<= lbGridR27.nE3;  k++) { 
    for (int j=lbGridR27.nB2; j<= lbGridR27.nE2;  j++) {
      for (int i=lbGridR27.nB1; i<= lbGridR27.nE1;  i++) { 
        x = (2.0*M_PI*(i-0.5)/(dataType)n1);
        y = (2.0*M_PI*(j-0.5)/(dataType)n2);
        z = (2.0*M_PI*(k-0.5)/(dataType)n3);

        u1 = U0*sin(x)*(cos(3.0*y)*cos(z) - cos(y)*cos(3.0*z)) ;
        u2 = U0*sin(y)*(cos(3.0*z)*cos(x) - cos(z)*cos(3.0*x)) ;
        u3 = U0*sin(z)*(cos(3.0*x)*cos(y) - cos(x)*cos(3.0*y)) ;

        getfEqIsothermal(feq, rho, u1,u2,u3, theta);
    
        copyToGrid(lbGridR27, feq, i, j, k);  
      }
    }
  }
}
