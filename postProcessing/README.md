This directory contains MATLAB files to perform post-processing of:

**LES sensitivity analysis**:
  - **pressureLES.m**: load LES data
  - **accuracyCFD.m**: compute accuracy of LES sensitivity analysis
  - **convergenceLES.m**: check convergence of LES simulation
  - **inflowConditionsLES.m**: compare velocity statistics outcome of LES sensitivity analysis
  - **inflowAveregeSensitivity.m**: compare pressure statistics outcome of LES sensitivity analysis

**RANS UQ**:
  - **accuracyCFD.m**: compute accuracy of UQ analysis
  - **pressureUQ.m**: load aleatoric UQ data
  - **pressureEpistemicUQ.m**: load epistemic UQ data 
  - **pressureCombinedUQ.m**: combine the two UQ studies

**PoliMi**:
  - **pressurePoliMi**: load PoliMi data

**WoW**:
  - **pressureWoW.m**: load WoW data
  - **pressureWoW_HR.m**: load high-speed WoW data
  - **applyTF.m**: apply transfer function
  - **importScanivalveData.m**: import WoW pressure data
  - **ReadTHFile.m**: import WoW velocity data
  - **transferFunctionWoW.m**: apply transfer function
    
**Main files**:
  - **highRisePressure.m**: main file to perform post-processing and validation 
  - **areaAveragedPressure.m**: compute area-averaged pressure across panel
  - **designPressure.m**: compute design pressure of panel
  - **computeTimeScale.m**: compute turbulence integral time-scale
  - **gumbel.m**: apply gumbel method
  - **pressureContour.m**: plot contour of pressure statistics
  - **pressureProfile.m**: plot profile of pressure statistics
  - **pressureInstantaneousContour.m**: plot contour of instantaneous pressure
  - **pressurePDF.m**: plot pressure probability density function
  - **pressurePSD.m**: plot pressure power spectra
  - **pressureTimeHistory.m**: plot pressure time series
    
    
