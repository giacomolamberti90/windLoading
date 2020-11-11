## Tutorial: generate inflow conditions for LES in empty domain
Download tutorial at: https://drive.google.com/drive/folders/1P5VEuYgWpWsbwFKwINNdYuxq2M4uXwbe?usp=sharing

The OpenFOAM case contains the following folders:
  - **0**: contains the initial and boundary conditions (no difference from standard OpenFOAM case)
  - **constant**: in addition to the standard dictionaries, the folder should contain:
      - **inflowProperties**: the input needed by the inflow solver 
      - **UInlet**: velocity profile
      - **uuBarInlet**: streamwise Reynolds stress profile
      - **vvBarInlet**: vertical Reynolds stress profile (y needs to be vertical direction for inflow solver!)
      - **wwBarInlet**: spanwise Reynolds stress profile
      - **uvBarInlet**: shear stress profile
  - **system**: contains the standard OpenFOAM files
  
To run a different test case follow these steps:
  1. Update mesh, i.e. domain dimension and mesh resolution
  2. Update the dictionary constant/inflowProperties containing the following entries:
      - **genZone**: name of the inflow generation plane
      - **refD**: reference height
      - **LY, LZ**: vertical and spanwise dimensions of the inflow generation plane
      - **NY, NZ**: vertical and spanwise number of cells of the virtual mesh used by the digital filter
      - **lagT_u, lagT_v, lagT_w**: streamwise, vertical and spanwise time-scales (average over height)
      - **Uref**: reference velocity
  3. Update the files containing inflow profiles (*Inlet); the files should contain two columns, i.e. vertical coordinate and value of the quantity of interest
  4. Update system/topoSetDict and run topoSet
  5. Run setsToZones -noFlipMap
  6. Run XCDF
     
