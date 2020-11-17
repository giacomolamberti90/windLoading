This directory contains MATLAB files to perform post-processing of:

  1. **LES sensitivity analysis**:
    - **pressureLES.m**: load LES data
    - **accuracyCFD.m**: compute accuracy of LES sensitivity analysis
    - **convergenceLES.m**: check convergence of LES simulation
    - **inflowConditionsLES.m**: compare velocity statistics outcome of LES sensitivity analysis
    - **inflowAveregeSensitivity.m**: compare pressure statistics outcome of LES sensitivity analysis

  2. **RANS UQ**:
    - **accuracyCFD.m**: compute accuracy of UQ analysis
    - **pressureUQ.m**: load aleatoric UQ data
    - **pressureEpistemicUQ.m**: load epistemic UQ data 
    - **pressureCombinedUQ.m**: combine the two UQ studies

  3. **PoliMi**:
    - **pressurePoliMi**: load PoliMi data

  4. **WoW**:
    - **pressureWoW.m**: load WoW data
    - **pressureWoW_HR.m**: load high-speed WoW data
    - **applyTF.m**: apply transfer function
    - **importScanivalveData.m**: import WoW pressure data
    - **ReadTHFile.m**: import WoW velocity data
    - **transferFunctionWoW.m**: apply transfer function
    
  5. **Main files**:
    - **highRisePressure.m**: main file to perform 
    - **areaAveragedPressure.m**:
    - **designPressure.m**:
    - **computeTimeScale.m**:
    - **gumbel.m**:
    - **pressureContour.m**:
    - **pressureProfile.m**:
    - **pressureInstantaneousContour.m**:
    - **pressurePDF.m**:
    - **pressurePSD.m**:
    - **pressureTimeHistory.m**:
    
    
