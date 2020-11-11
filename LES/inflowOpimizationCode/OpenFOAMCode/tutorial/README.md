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
