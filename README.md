# WINDEX: Scanning whole-genome aligned haplotype data for signatures of positive selective sweeps using an HHMM
_Authors:_ Hannah Snell, Scott McCallum, Dhruv Raghavan, Ritambhara Singh, Sohini Ramachandran, Lauren Sugden

**Will be formally packaged upon publication.** 

## Getting started
To use WINDEX, please set up the following environment based on the `requirements.txt` file in this repository: 

1. Make a new python environment called `windex_env`:  
  ```sh
  python3 -m venv ~/.venv/windex_env
  ```  

2. Activate said environment:  
  ```sh
  source ~/.venv/windex_env/bin/activate
  ```  

3. Install all needed dependencies with the `requirements.txt` file:  
  ```sh
  pip install -r requirements.txt
  ```  

## Input data

### Site- and window-based statistics

WINDEX requires aligned haplotype data from human genome samples as input, specifically in the form of calculated summary statistics. The easiest way to achieve this is to calculate your choice of selection statistics from a VCF file containing a number of individuals from a known population. We have provided a site-based statistics file and a window based-statistics file as an example in `examples/emission_training`. Both files must follow the same tab-delimited input format as SWIF(r), which is:  

```sh
SNP_name \t Physical_Distance \t Map_Distance \t Stat1Name \t Stat2Name \t ...
``` 

The `SNP_name` column can contain either variant positions, such as for site-based statistical calculations, or window start positions, such as for window-based statistical calculations. **Important: do not mix site- and window-based statistics in the same input file; they must stay separate.** WINDEX does not use the `Physical_Distance` or `Map_Distance` columns directly, so these can contain placeholder numbers to satisfy the file compatibility between SWIF(r) and WINDEX. 

### Emission training

WINDEX also requires a set of trained Gaussian mixture models that represent the marginal distributions of each input statistic over each class (neutral/linked/sweep). Although it is not required, the easiest way to obtain these files is by training these marginal distributions using SWIF(r), since the statistic files are compatible. We also provide an easy step-by-step training demonstration with SWIF(r) as an example in `examples/emission_training` to demystify this process.

## Running WINDEX 

WINDEX can be called using this python command:

```sh
python hierarchical_hmm.py path_to_trained_sites path_to_trained_windows site_file_output_name window_file_output_name site_statistics_file window_statistics_file window_size=100000
```

Please see FAQ page for any warnings/errors that arise during WINDEX testing.

## Simple example 

## Whole-genome scans in data from the 1000 Genomes Project
