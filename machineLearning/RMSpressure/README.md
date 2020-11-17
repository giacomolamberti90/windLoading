This directory contains the code relative to [1]:
  - **train.py**: training of machine learning models
  - **evaluate.py**: evaluation of machine learning models
  - **model.py**: deep learning model
  - **loader.py**: load data into dataframes
  - **utils.py**: functions to perform post-processing
  - **features.py**: load and compute features from RANS
  - **labels.py**: load rms pressure from LES and rotate coordinates

The weigths of the models at different wind directions can be found here: https://drive.google.com/drive/folders/1URstkFLY5JE210xiOb-Zypoxix9D6boW?usp=sharing

The data used to train/test the models can be found here: https://drive.google.com/drive/folders/1j_1n68Zv536TzPycXgckk86icpchRozQ?usp=sharing

The LES simulations used to obtain the rms presssure coefficient can be found here: https://drive.google.com/drive/folders/1nLW-Y6TaWN6vO9kB0KKBfXEzjeQ50tTg?usp=sharing

The RANS simulations used to obtain the features can be found here: https://drive.google.com/drive/folders/1UB2d9xUhWXW1iDNxZHWBlpwpQPx4GSbF?usp=sharing

[1] Lamberti, Giacomo, and Catherine Gorlé. "A multi-fidelity machine learning framework to predict wind loads on buildings" Journal of Wind Engineering and Industrial Aerodynamics (under review).
