#include <iostream>
#include <math.h>

int main(){
  int n=64;
  double kL = 2.*M_PI/n;
  double U0 = 0.25;
  double nu = 1./6.;
  double ke = 0.;

  for (int j=0;j<=n;j++){
    for (int i=0;i<=n;i++){

      double x = i;
      double y = j;

      double rho = 1. - (3./4.*U0*U0*(cos(2.*kL*x) + cos(2.*kL*y)));
      double u1 = U0*cos(kL*x)*sin(kL*y);
      double u2 = -U0*sin(kL*x)*cos(kL*y);
      double u3 = 0.;

      std::cout << "At " << kL*x << " " << kL*y << " " << rho << " " << u1 << " " << u2 << std::endl;
      ke += 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
    }
  }

  std::cout << "KE " << ke << " " << U0*U0*M_PI*M_PI << std::endl;
  return 0;
}