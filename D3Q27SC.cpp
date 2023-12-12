
template <int numfield, typename dataType> 
inline void getGradHeatOptSplit(gridBCC3D< numfield,dataType>  &lbGridR27, int i1, int i2, int i3, double &rho, double &u1, double &u2, double &u3, double &theta, double &q1, double &q2,double &q3)
{
 double csq,rhoInv;
 rho=0.0; u1=0.0; u2=0.0; u3=0.0;theta=0.0;q1=0.0;q2=0.0;q3=0.0;
 for(int dv=0;dv<27;dv++)
 {
  rho += lbGridR27.Node(i1,i2,i3,dv);
  u1 += lbGridR27.Node(i1,i2,i3,dv)*cx[dv];
  u2 += lbGridR27.Node(i1,i2,i3,dv)*cy[dv];
  u3 += lbGridR27.Node(i1,i2,i3,dv)*cz[dv]; 
  csq = cx[dv]*cx[dv] + cy[dv]*cy[dv] + cz[dv]*cz[dv];
  theta += lbGridR27.Node(i1,i2,i3,dv)*csq;
 }
 rhoInv = 1.0/rho;
 u1 *= rhoInv; u2 *= rhoInv; u3 *= rhoInv;
 theta -= rho*(u1*u1 + u2*u2 + u3*u3);
 theta /= 3.0*rho;    
 for(int dv=0;dv<27;dv++)
 {
  csq = cx[dv]*cx[dv] + cy[dv]*cy[dv] + cz[dv]*cz[dv];     
  q1 += lbGridR27.Node(i1,i2,i3,dv)*csq*cx[dv];   
  q2 += lbGridR27.Node(i1,i2,i3,dv)*csq*cy[dv]; 
  q3 += lbGridR27.Node(i1,i2,i3,dv)*csq*cz[dv];      
 }  
 double factor=(1.0+2.0*theta);
 q1 = q1*rhoInv-u1*factor;
 q2 = q2*rhoInv-u2*factor;
 q3 = q3*rhoInv-u3*factor; 
} 

/////////////
// Get Feq //
/////////////
inline void getfEqIsothermal(double *Feq,double rho, double u1, double u2, double u3, double theta)
{
 double v1, v2, v3, tmp, dot, T0Inv;
 T0Inv = 3.0;
 v1 = u1*T0Inv;
 v2 = u2*T0Inv;
 v3 = u3*T0Inv;
 double usq = u1*v1 + u2*v2+v3*u3;
 tmp =1.0-usq*0.5;
 Feq[0] = wt[0]*rho*tmp;
 for (int dv=1;dv<27;dv++)
 {
  dot = (v1*cx[dv] + v2*cy[dv] + v3*cz[dv]);
  Feq[dv]= wt[dv]*rho*(tmp +dot*(1.0 + 0.5* dot));
 }
}

inline void getF(double *F,double rho, double u1, double u2, double u3, double theta)
{
 double dot(0.0),term1(0.0),cs2(1.0),cSqr(0.0),Ki(0.0),index(0.0),wt(0.0);
 u1 = u1/theta;      
 u2 = u2/theta;      
 u3 = u3/theta;  
 double oneMinusTheta = (1.0-theta);
 for(int dv=0;dv<27;dv++)
 {
  cSqr = cx[dv]*cx[dv] + cy[dv]*cy[dv] + cz[dv]*cz[dv];
  index = cSqr/cs2;
  dot = (u1*cx[dv] + u2*cy[dv] + u3*cz[dv]);
  wt = oneMinusTheta*oneMinusTheta*oneMinusTheta* pow((theta/(2.0*oneMinusTheta)),index);
  Ki = (6.0*theta*theta + (cSqr/cs2)*(1.0-3.0*theta)) / (3.0*(1.0-theta));  
  term1 = 0.5*dot*dot - (0.5*Ki*(u1*u1+u2*u2+u3*u3));
  F[dv] = rho*wt*(1.0 + dot + term1);    
 }    
}



///////////////////
//  Get Moments  //
///////////////////
template <int numfield, typename dataType> 
void getMomentsNode(gridBCC3D< numfield,dataType>  &lbGridR27, int i1, int i2, int i3, double &rho, double &u1, double &u2, double &u3, double &theta)
{
 double csq,rhoInv;
 rho=0.0; u1=0.0; u2=0.0; u3=0.0;theta=0.0;
 for(int dv=0;dv<27;dv++)
 {
  rho += lbGridR27.Node(i1,i2,i3,dv);
  u1  += lbGridR27.Node(i1,i2,i3,dv)*cx[dv];
  u2  += lbGridR27.Node(i1,i2,i3,dv)*cy[dv];
  u3  += lbGridR27.Node(i1,i2,i3,dv)*cz[dv]; 
  csq = cx[dv]*cx[dv] + cy[dv]*cy[dv] + cz[dv]*cz[dv];
  theta += lbGridR27.Node(i1,i2,i3,dv)*csq;
 }
 rhoInv = 1.0/rho;
 u1 *= rhoInv; u2 *= rhoInv;  u3 *= rhoInv;  
 theta -= rho*(u1*u1 + u2*u2 + u3*u3);
 theta /= (3.0*rho);
}

void getMoments(double *fTemp, double &rho, double &u1, double &u2, double &u3, double &theta)
{
 double csq,rhoInv;
 rho=0.0; u1=0.0; u2=0.0; u3=0.0;theta=0.0;
 for(int dv=0;dv<27;dv++)
 {
  rho += fTemp[dv];
  u1  += fTemp[dv]*cx[dv];
  u2  += fTemp[dv]*cy[dv];
  u3  += fTemp[dv]*cz[dv]; 
  csq = cx[dv]*cx[dv] + cy[dv]*cy[dv] + cz[dv]*cz[dv];
  theta += fTemp[dv]*csq;
 }
 rhoInv = 1.0/rho;
 u1 *= rhoInv; u2 *= rhoInv;  u3 *= rhoInv;  
 theta -= rho*(u1*u1 + u2*u2 + u3*u3);
 theta /= (3.0*rho);
}

template <typename dataType> 
void getPabMoments(dataType* fTemp, double &rho, double &ux, double &uy, double &uz, double &pxx, double &pyy ,  double &pzz, double &pxy, double &pyz, double &pxz)
{ 
  double csq,rhoInv;
  rho=0.0; ux=0.0; uy=0.0; uz=0.0;pxx=0.0;pyy=0.0;pzz=0.0;pxy=0.0;pxz=0.0;
  pyz=0.0;
  for(int dv=0;dv<27;dv++)
  {
    rho += fTemp[dv];
    ux  += fTemp[dv]*cx[dv];
    uy  += fTemp[dv]*cy[dv];
    uz  += fTemp[dv]*cz[dv]; 
    pxx += fTemp[dv]*cx[dv]*cx[dv];
    pyy += fTemp[dv]*cy[dv]*cy[dv];
    pzz += fTemp[dv]*cz[dv]*cz[dv];
    pxy += fTemp[dv]*cx[dv]*cy[dv];
    pxz += fTemp[dv]*cx[dv]*cz[dv];
    pyz += fTemp[dv]*cy[dv]*cz[dv];  
  }
  rhoInv = 1.0/rho;
  ux *= rhoInv; uy *= rhoInv; uz *= rhoInv;
}

/////////////////////
//  Get Enstrophy  //
/////////////////////
template <int numfield, typename dataType> 
void getEnstrophyMomentsNode(gridBCC3D< numfield,dataType>  &lbGridR27, int i1, int i2, int i3, double &rho, double &ux, double &uy, double &uz, double &theta, double &pxx, double &pyy ,  double &pzz, double &pxy, double &pxz, double &pyz, double &sXX, double &sYY, double &sZZ, double &sXY, double &sYZ, double &sXZ, double beta, double dt)
{ 
  double csq,rhoInv;
  rho=0.0; ux=0.0; uy=0.0; uz=0.0;theta=0.0;pxx=0.0;pyy=0.0;pzz=0.0;pxy=0.0;pxz=0.0;
  pyz=0.0,sXX=0.0;sYY=0.0;sZZ=0.0;sXY=0.0;sXZ=0.0;sYZ=0.0;
  for(int dv=0;dv<27;dv++)
  {
    rho += lbGridR27.Node(i1,i2,i3,dv);
    ux += lbGridR27.Node(i1,i2,i3,dv)*cx[dv];
    uy += lbGridR27.Node(i1,i2,i3,dv)*cy[dv];
    uz += lbGridR27.Node(i1,i2,i3,dv)*cz[dv]; 
    csq = cx[dv]*cx[dv] + cy[dv]*cy[dv] + cz[dv]*cz[dv];
    theta += lbGridR27.Node(i1,i2,i3,dv)*csq;
    pxx += lbGridR27.Node(i1,i2,i3,dv)*cx[dv]*cx[dv];
    pyy += lbGridR27.Node(i1,i2,i3,dv)*cy[dv]*cy[dv];
    pzz += lbGridR27.Node(i1,i2,i3,dv)*cz[dv]*cz[dv];
    pxy += lbGridR27.Node(i1,i2,i3,dv)*cx[dv]*cy[dv];
    pxz += lbGridR27.Node(i1,i2,i3,dv)*cx[dv]*cz[dv];
    pyz += lbGridR27.Node(i1,i2,i3,dv)*cy[dv]*cz[dv];  
  }
  rhoInv = 1.0/rho;
  ux *= rhoInv; uy *= rhoInv; uz *= rhoInv;
  theta -= rho*(ux*ux + uy*uy + uz*uz);
  theta /= 3.0*rho;   
  double factor = beta/(dt*rho*theta);
  sXX = (rho*ux*ux+rho*theta-pxx)*factor;
  sYY = (rho*uy*uy+rho*theta-pyy)*factor;
  sZZ = (rho*uz*uz+rho*theta-pzz)*factor;
  sXY = (rho*ux*uy- pxy)*factor;
  sXZ = (rho*ux*uz- pxz)*factor;
  sYZ = (rho*uy*uz- pyz)*factor;
}

template <typename dataType> 
void getEnstrophyMoments(dataType* fTemp, dataType &rho, dataType &ux, dataType &uy, dataType &uz
, dataType &theta, dataType &pxx, dataType &pyy ,  dataType &pzz, dataType &pxy, dataType &pxz, dataType &pyz, dataType &sXX, dataType &sYY, dataType &sZZ
, dataType &sXY, dataType &sYZ, dataType &sXZ, dataType beta, dataType dt)
{ 
  double csq,rhoInv;
  rho=0.0; ux=0.0; uy=0.0; uz=0.0;theta=0.0;pxx=0.0;pyy=0.0;pzz=0.0;pxy=0.0;pxz=0.0;
  pyz=0.0,sXX=0.0;sYY=0.0;sZZ=0.0;sXY=0.0;sXZ=0.0;sYZ=0.0;

  for(int dv=0;dv<27;dv++)
  {
    rho += fTemp[dv];
    ux  += fTemp[dv]*cx[dv];
    uy  += fTemp[dv]*cy[dv];
    uz  += fTemp[dv]*cz[dv]; 
    csq = cx[dv]*cx[dv] + cy[dv]*cy[dv] + cz[dv]*cz[dv];
    theta += fTemp[dv]*csq;
    pxx += fTemp[dv]*cx[dv]*cx[dv];
    pyy += fTemp[dv]*cy[dv]*cy[dv];
    pzz += fTemp[dv]*cz[dv]*cz[dv];
    pxy += fTemp[dv]*cx[dv]*cy[dv];
    pxz += fTemp[dv]*cx[dv]*cz[dv];
    pyz += fTemp[dv]*cy[dv]*cz[dv];  
  }
  rhoInv = 1.0/rho;
  ux *= rhoInv; uy *= rhoInv; uz *= rhoInv;
  theta -= rho*(ux*ux + uy*uy + uz*uz);
  theta /= 3.0*rho;   
  double factor = beta/(dt*rho*theta);
  sXX = (rho*ux*ux+rho*theta-pxx)*factor;
  sYY = (rho*uy*uy+rho*theta-pyy)*factor;
  sZZ = (rho*uz*uz+rho*theta-pzz)*factor;
  sXY = (rho*ux*uy- pxy)*factor;
  sXZ = (rho*ux*uz- pxz)*factor;
  sYZ = (rho*uy*uz- pyz)*factor;
}

template<typename dataType>
void getGrad(dataType *fGrad, dataType rho, dataType ux, dataType uy, dataType uz, dataType Pxx, dataType Pyy, dataType Pzz, dataType Pxy, dataType Pyz, dataType Pxz) 
{
  dataType theta0 = 1./3.;
  for (int dv=0;dv<27;dv++) {
    double term1 = rho;
    double term2 = rho*(ux*cx[dv] + uy*cy[dv] + uz*cz[dv])/theta0;
    Pxx = Pxx - rho*theta0;
    Pyy = Pyy - rho*theta0;
    Pzz = Pzz - rho*theta0;
    double term3 = (1./(2.*theta0*theta0))*(Pxx*cx[dv]*cx[dv] + Pyy*cy[dv]*cy[dv] + Pzz*cz[dv]*cz[dv] + 2.*Pxy*cx[dv]*cy[dv] + 2.*Pyz*cy[dv]*cz[dv] + 2.*Pxz*cx[dv]*cz[dv]);
    double term4 = (1./(2.*theta0))*(Pxx + Pyy + Pzz);
    fGrad[dv] = wt[dv]*(term1 + term2 + term3 - term4);
  }
}

//  double rho(0.97834), u1(-0.006234), u2(0.023742), u3(0.01234523), theta(0.33333), feq[27];
 
//  std::cout << "Before " << rho << " " << u1 << " " << u2 << " " << u3 << " " << theta << std::endl;
//  getfEqIsothermal(feq, rho, u1, u2, u3, theta);
//  getMoments(feq, rho, u1, u2, u3, theta);
//  std::cout << "After " << rho << " " << u1 << " " << u2 << " " << u3 << " " << theta << std::endl;

// std::cout << "Before " << rho << " " << ux << " " << uy << " " << uz << " " << Pxx << " " << Pyy << " " << Pzz << " " << Pxy << " " << Pyz << " " << Pxz << std::endl;
//   getGrad(fGrad, rho, ux, uy,  uz, Pxx,  Pyy,  Pzz,  Pxy, Pyz,  Pxz);
//   getPabMoments(fGrad, rho, ux, uy, uz, Pxx,  Pyy,  Pzz,  Pxy, Pyz,  Pxz);
//   std::cout << "After " << rho << " " << ux << " " << uy << " " << uz << " " << Pxx << " " << Pyy << " " << Pzz << " " << Pxy << " " << Pyz << " " << Pxz << std::endl;