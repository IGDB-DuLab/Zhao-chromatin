# Mapping Chromatin activity landscape during *C. elegans* embryogenesis at Single-cell resolution

This repository provides the source code for reproducing the results of the study titled " Chromatin dynamics during *C. elegans* embryogenesis at the single-cell level" by Zhiguang Zhao# , Rong Fan# , Weina Xu, Yangyang Wang, Xuehua Ma, and Zhuo Du*.

## Description of files
- [fluorescent_compensation.py](https://github.com/genetics-dulab/scCAL/blob/main/fluorescent_compensation.py)  Compensation for depth-dependent attenuation of fluorescence intensity during 3D live-cell imaging.
- [data_stats.py](https://github.com/genetics-dulab/scCAL/blob/main/data_stats.py) Calculating the genomic and cellular distribution, coefficient of variation, the information content of GFP expression integrated into different genomic positions.
- [CAL_lineage.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_lineage.py) Performing analyses related to lineage-dependent changes in chromatin activity landscape.
- [CAL_tissue.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_tissue.py) Performing analyses related to tissue-based convergence and lineage-dependent heterogeneity in chromatin landscape within the same tissue type.
- [CAL_symmetry.py](https://github.com/genetics-dulab/scCAL/blob/main/CAL_symmetry.py) Performing analyses related to chromatin landscape divergences between left-right symmetry cells and their progenitor cells. 
- [CDCA.py](https://github.com/genetics-dulab/scCAL/blob/main/CDCA.py) Performing analyses related to chromatin co-dynamic region.
- [data](https://github.com/genetics-dulab/scCAL/blob/main/data) folder contains the raw data needed for various analyses.

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

