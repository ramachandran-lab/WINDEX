# WINDEX: Scanning whole-genome aligned haplotype data for signatures of positive selective sweeps using an HHMM
_Authors:_ Hannah Snell, Scott McCallum, Dhruv Raghavan, Ritambhara Singh, Sohini Ramachandran, Lauren Sugden

**Will be formally packaged upon publication.** 

Processed VCFs are being uploaded to a Zenodo repository. Will update shortly. 

## Getting started
To use WINDEX, please start by cloning this repository: 

```sh
git clone https://github.com/hmsnell/WINDEX.git
```

Next, set up the following environment based on the `requirements.txt` file in this repository. 
WINDEX is compatible with Python versions 3.10 through 3.12: 

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
python hierarchical_hmm.py path_to_trained_sites path_to_trained_windows site_file_output_name window_file_output_name site_statistics_file window_statistics_file window_size
```

Arguments: 
- `path_to_trained_sites`: path to dir where the trained emissions are for your site statistics
- `path_to_trained_windows`: path to dir where the trained emissions are for your window statistics
- `site_file_output`: name of desired output for site-level observations
- `window_file_output`: name of desired output for window-level observations
- `site_statistics`: path to site statistic file for the simulation/genome in question
- `window_statistics`: path to window statistic file for the simulation/genome in question
- `window_size`: window size used to calculate window-based statistics

Please see FAQ page for any warnings/errors that arise during WINDEX testing.

## Simple example 

Here we outline a simple and intuitive example for how to preprocess the data, run WINDEX, and do a post-hoc analysis. We will use a subset of simulations generated for WINDEX validation in the paper. These simulations are specifically representing 1Mbp generic contigs in the Yoruba Nigerian (YRI) population under the Three Population out of Africa demographic model (Gutenkunst _et al.,_ 2009) with the introduction of a positive selective sweep 300 generations before present, with a selection coefficient of 0.05, at position 500001. We already calculated all of the site- and window-based statistics (in 40kb-sized windows) for you - they are the same as the ones reported in the WINDEX paper. You can explore these files in the following directories: 

```sh
~/example/emission_training/site_stats/stat_files/
# and 
~/example/emission_training/windowstats/stat_files/
```
### Step One: WINDEX Emission Training 

We will start by splitting these files into training and testing sets so that we can use these files to train the emissions needed for WINDEX. The easiest way to do this training is by using SWIF(r) since all files are compatible as input. First, use the script `split_train_test.sh` to generate a single training and testing file for each class (neutral, linked, sweep) that contains all concatenated summary statistics from the simulations. **Important: directories `site_stats/` and `window_stats/` are configured for SWIF(r) training already by their file structure. Please note this file organization for your own experiments and do not alter this structure for this example.**

You will know that you ran `split_train_test.sh` correctly when you can go into the `simulations/` directory in either `site_stats/` or `window_stats/` and see a single file in each class folder with the extension `.train`. You can also check that the testing data was deposited in the `testing_data/` folders, too. 

Next, we can run the SWIF(r) training by using the script `train_distributions.sh`. For each statistic level (site and window), you will use this script to launch the training. Once the training is done, you can optionally (but still recommended) retrain SWIF(r) so that the marginal distributions used for the emissions have only 1 component each. Do this by going into the `AODE/` directories of both `site_stats/` and `window_stats/` and manually editing both the `marginal_component_nums` and the `joint_component_nums` files so that all columns contain 1's. Then, you can use the same `train_distributions.sh` script to launch SWIF(r) retraining with the appropriate command. 

Once training (and retraining) are complete at both the site- and window-levels, you are now ready to run WINDEX! 

## Step Two: Running WINDEX 

WINDEX is easy to run using the `run_windex.sh` script found in `~/examples/run_windex` folder. The required inputs should already be configured for this example. You will know that you ran this correctly when the `output/` directory contains files with the extensions `.site_classified` and `.window_classified`. 

Now that you've run WINDEX, you can explore the post-hoc analyses we used in our paper by loading the output files in the the `post_hoc_analysis.Rmd` file. 

## WINDEX Transition Probabilities 

## Whole-genome scans in data from the 1000 Genomes Project
