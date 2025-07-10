# SM_QAQC_Results

Small and basic repo for generating performance plots from data acquired by the Sensor Module test stand used for BTL production. Code for analyzing raw waveform-level data, integrated waveform values, and analysis summary plots in root files.

## Setup
Code in run in a conda environment on Caltech Tier-II machines with el8 linux distribution. After cloning the repo, basic environment with ROOT can be created, activated, and have necessary packages install via:

```
conda env create --file qaqc_env.yml
conda activate qaqc
conda install -c conda-forge hdf5 numpy scipy psycopg2 root h5py matplotlib seaborn
``` 
