# Single-cell Chromatin Acitivity Landscape Analysis

This repository contains the source code for reproducing the analyses of our study titled "Multidimensional Chromatin-based Regulation of *C. elegans* Embryogenesis" by Zhiguang Zhao# , Rong Fan# , Weina Xu, Yangyang Wang, Xuehua Ma, and Zhuo Du* 

## Description of files
- [fluorescent_compensation.py](https://github.com/genetics-dulab/scCAL/blob/main/fluorescent_compensation.py)  Compensation for depth-dependent attenuation of fluorescence intensity
- [data_stats.py](https://github.com/genetics-dulab/scCAL/blob/main/data_stats.py) Calculation of data distribution, CV, infomation content and variability.
- [CAL_lineage.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_lineage.py) Performing analyses related to lineage fate patterning and anterior-posterior fate asymmetry.
- [CAL_tissue.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_tissue.py) Performing analyses related to tissue convergence and tissue heterogeneity.
- [CAL_symmetry.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_symmetry.py) Performing analyses related to symmetry predetermination.
- [CDCA.py](https://github.com/genetics-dulab/scCAL/blob/main/CDCA.py) Performing analyses related to chromatin co-dynamic region.
- [data](https://github.com/genetics-dulab/scCAL/blob/main/data) folder contains the raw data needed for downstream analysis.

## Prerequisites
All the dependencies are listed in [requirements.txt](https://github.com/genetics-dulab/scCAL/blob/main/requirements.txt). 

Most analyses were performed under python 3.7.1

Additional Python packages required for analysis include:

- pandas (>= 1.0.3)
- numpy (>= 1.19.1)
- scipy (>= 1.2.1)
- statsmodels (>= 0.11.1)

Additional Python packages required for plotting the figures include:

- seaborn (>= 0.11.0)
- matplotlib (>= 3.0.2)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](https://github.com/genetics-dulab/scCAL/blob/main/LICENSE) file for details.

