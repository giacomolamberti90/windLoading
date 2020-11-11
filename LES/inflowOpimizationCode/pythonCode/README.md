## Inflow Optimization 
This directory contains the python code to perform gradient-based optimization of Reynolds stresses and integral time-scales (refer to [1] for more details):
  - **optXCDF.py**: main file used to perform the inflow optimization; the file first initializes the quantities of interest (i.e. Reynolds stresses and time-scales at inflow plane and building locations), then read the results of the two initial simulations (test0, test1) and finally runs the optimization to produce the input files for test2. The code automatically generates the dictionaries necessary to run the OpenFOAM simulations.
  - **multiObjGrad.py**: contains functions to perform one step of optimization for Reynolds stresses and time-scales
  - **runOpenFOAM.py**: (optional) automatically run OpenFOAM with the new files and check for simulation convergence; the function should be called from optXCDF.py.
  - **inputData.py**: defines input data, i.e. experimental values, Bezier points, hyper-parameters of optimization, folder paths, etc.
  - **bezier.py**: function that generates Bezier curve
  - **timeScale.py**: function to compute the integral time-scale
  - **postProcessing.py**: plot Reynolds stresses and time-scales throughout optimization steps
  - **exp**: folder containing experimental data

[1] Lamberti, Giacomo, et al. "Optimizing turbulent inflow conditions for large-eddy simulations of the atmospheric boundary layer." Journal of Wind Engineering and Industrial Aerodynamics 177 (2018): 32-44.
