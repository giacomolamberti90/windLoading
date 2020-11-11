## Inflow condition
This folder contains the OpenFOAM code that generates inflow conditions for LES (see [1,2] for more details):
  - **XCDF.C**: pisoFoam solver + divergence free inflow condition (analogous to pisoFoam.C)
  - **Ugen.H**: digital filter procedure
  - **pEqn.H**: introduces the turbulent flow field in a plane next to the inlet inside the domain, then solve the standard Poisson equation
  - **UEqn.H**: standard momentum equation
  - **inflowProperties.H**: read dicti
  - **createFields.H**: define variables 
  - **Make**: contains files to compile the solver
  - **tutorial**: contains link to a sample directory which contains all the files needed to run a simulation using the inflow solver
 
[1] Xie, Zheng-Tong, and Ian P. Castro. "Efficient generation of inflow conditions for large eddy simulation of street-scale flows." Flow, turbulence and combustion 81.3 (2008): 449-470.

[2] Kim, Yusik, Ian P. Castro, and Zheng-Tong Xie. "Divergence-free turbulence inflow conditions for large-eddy simulations with incompressible flow solvers." Computers & Fluids 84 (2013): 56-68.
