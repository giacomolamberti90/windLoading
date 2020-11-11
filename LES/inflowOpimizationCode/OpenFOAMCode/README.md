## Inflow condition
This folder contains the OpenFOAM code that generates inflow conditions for LES (see [1,2] for more details):
  - **XCDF.C**: standard pisoFoam solver
  - **Ugen.H**: digital filter procedure
  - **pEqn.H**: introduces the turbulent flow field in a plane inside the domain next to the inlet boundary, then solve the standard Poisson equation
  - **UEqn.H**: standard momentum equation
  - **inflowProperties.H**: define inflow properties to be read from OpenFOAM dictionary 
  - **createFields.H**: define variables needed by the digital filter
  - **Make**: contains files to compile the solver
  - **tutorial**: contains link to a sample directory and instructions on how to run the simulation
 
[1] Xie, Zheng-Tong, and Ian P. Castro. "Efficient generation of inflow conditions for large eddy simulation of street-scale flows." Flow, turbulence and combustion 81.3 (2008): 449-470.

[2] Kim, Yusik, Ian P. Castro, and Zheng-Tong Xie. "Divergence-free turbulence inflow conditions for large-eddy simulations with incompressible flow solvers." Computers & Fluids 84 (2013): 56-68.
