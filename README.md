# scCAL
## Description of files
- [fluorescent_compensation.py](https://github.com/genetics-dulab/scCAL/blob/main/fluorescent_compensation.py) source code for compensation for depth-dependent attenuation of fluorescence intensity
- [data_stats.py](https://github.com/genetics-dulab/scCAL/blob/main/data_stats.py) source code for data distribution, CV, infomation content, variability
- [CAL_lineage.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_lineage.py) source code related to lineage analysis
- [CAL_tissue.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_tissue.py) source code related to tissue analysis
- [CAL_symmetry.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_symmetry.py) source code related to symmetry analysis
- [CDCA.py](https://github.com/genetics-dulab/scCAL/blob/main/CDCA.py) source code related to CDAC analysis

## Prerequisites
All the dependencies are listed in requirements.txt. 

Most analyses were performed under python 3.7.1

Additional Python packages required for analysis include:

- pandas 1.0.3
- numpy 1.19.1
- scipy 1.2.1
- statsmodels 0.11.1

Additional Python packages required for plotting the figures include:

- seaborn 0.11.0
- matplotlib 3.0.2