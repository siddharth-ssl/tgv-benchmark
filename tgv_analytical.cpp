// The analytical is taken from the following article
// Advanced Automatic Code Generation for Multiple Relaxation-Time Lattice Boltzmann Methods
// Authors: Frederik Hennig, Markus Holzer, and Ulrich Rüde
// SIAM Journal on Scientific Computing 2023 45:4, C233-C254 

template <int numfield, typename dataType>
void l2norm(gridBCC3D<numfield,dataType> &gridLB, dataType U0, dataType nu, int t)
{
  dataType error_ke(0.), error_u1(0.), error_u2(0.), error_rho(0.);
  dataType fTemp[numfield], rho, ux, uy, uz, theta, ke(0.);
  dataType x, y, kL=2.*M_PI/gridLB.m1;
  dataType tgv_ux, tgv_uy, tgv_rho, tgv_ke(0.);
  dataType e_u1, e_u2, e_rho, e_ke;

  for (int k=gridLB.nB3; k<= gridLB.nE3;  k++) { 
    for (int j=gridLB.nB2; j<= gridLB.nE2;  j++) { 
		  for (int i=gridLB.nB1; i<= gridLB.nE1;  i++) {

        copyFromGrid(gridLB, fTemp, i, j, k);
        getMoments(fTemp, rho, ux, uy, uz, theta);

		    dataType ke_local = 0.5*rho*(ux*ux + uy*uy + uz*uz);
        ke += ke_local;

        x = (i);
        y = (j);
			
        tgv_ux = U0*cos(kL*x)*sin(kL*y)*exp(-2.*kL*kL*nu*t); // ux
			  tgv_uy = -U0*sin(kL*x)*cos(kL*y)*exp(-2.*kL*kL*nu*t); // uy
			  tgv_rho = 1. - 0.75*U0*U0*(cos(2*kL*x) + cos(2*kL*y))*exp(-4*nu*kL*kL*t); // rho

        dataType tgv_ke_local = 0.5*rho*(tgv_ux*tgv_ux + tgv_uy*tgv_uy);// KE
        tgv_ke += tgv_ke_local;

        e_u1 = ux - tgv_ux;
        e_u2 = uy - tgv_uy;
        e_rho = rho - tgv_rho;
        e_ke = ke_local - tgv_ke_local;

        error_u1 += e_u1*e_u1;
        error_u2 += e_u2*e_u2;
        error_rho += e_rho*e_rho;
        error_ke += e_ke*e_ke;
      }
    }
  }

	// dataType tgv_ke = U0*U0*M_PI*M_PI*exp(-4*nu*kL*kL*t); // KE
  std::cout << "TGV KE at time t=" << t << " is " << tgv_ke << std::endl;
  std::cout << "LB KE at time t=" << t << " is " << ke << std::endl;

  error_u1 = std::sqrt(error_u1/gridLB.m1*gridLB.m2*gridLB.m3);
  error_u2 = std::sqrt(error_u2/gridLB.m1*gridLB.m2*gridLB.m3);
  error_rho = std::sqrt(error_rho/gridLB.m1*gridLB.m2*gridLB.m3);
  error_ke = std::sqrt(error_ke/gridLB.m1*gridLB.m2*gridLB.m3);

  std::cout << "L2 Norm in u1 = " << error_u1 << std::endl;
  std::cout << "L2 Norm in u2 = " << error_u2 << std::endl;
  std::cout << "L2 Norm in rho = " << error_rho << std::endl;
  std::cout << "L2 Norm in ke = " << error_ke << std::endl;
}
