This code computes the density and velocity in 3D space for taylor-green vortex initialization (TGV).
The details of the TGV (the analytical expressions) can be found in the following article \
**Advanced Automatic Code Generation for Multiple Relaxation-Time Lattice Boltzmann Methods
Authors: Frederik Hennig, Markus Holzer, and Ulrich Rüde
SIAM Journal on Scientific Computing 2023 45:4, C233-C254**

Compile the code using the following steps
```bash
g++ -I . -O3 main.cpp -O tgv.out

or 

bash compile.sh
```

Execute the code using the following steps
```bash
./tgv.out 64 64 64 0.05 100 5 0.1

or 

bash run.sh
```

The inputs in the above line are in the following order
```bash
./tgv.out n1 n2 n3 Mach Reynolds iterations print_freq
```
where,
- n1 = Number of grid points in x-direction
- n2 = Number of grid points in y-direction
- n3 = Number of grid points in z-direction
- Mach = Mach number of the flow (used to calculate the initial kinetic energy)
- Reynolds = Reynolds number of the flow (used to calculate Viscosity of the fluid, with charLen=1)
- iterations = The number specified here is in terms of convection times. Convection time is defined as `$n1/(2\pi U0)$`. So if iterations = 5, it means the simulation would run for 5 convection times
- print_freq = The frequency at which the L2 norm w.r.t analytical expression for TGV is printed and the field data consisting velocities is printed. This is also expressed in terms of convection times. If print_freq = 0.1, it means data will be printed every 0.1 convection times.
