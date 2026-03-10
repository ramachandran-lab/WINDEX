#!/bin/bash

### SLURM suggested params

#SBATCH -J run_windex 	   
#SBATCH -N 1 						
#SBATCH -n 1                        
#SBATCH -t 2-00:00 			        	
#SBATCH --mem 10G 				    

<< read_me
This script can be used to execute WINDEX after training the GMMs on your statistic files. 
Input params (in order): 
path_to_trained_sites = path to directory where trained distributions for site statistics were generated.
path_to_trained_windows = path to directory where trained distributions for window statistics were generated.
outputname_sites = desired output name for WINDEX site classification file.
outputname_windows = desired output name for WINDEX window classification file.
site_stat_file = path to site statistic file for one simulation/genome.
window_stat_file = path to window statistic file for one simulation/genome.
window_size = size of window in which the window statistics were calculated. 
read_me

module load python
source ~/venv/windex_env/bin/activate

### arguments in order ### 
    # path to training dir for sites
    # path to training dir for windows 
    # output file name for sites
    # output file name for windows
    # path to site statisitc file for single simulation/genome
    # path to window statistic file for single simulation/genome
    # window size used to calculate window statistics

# run WINDEX

for i in {81..100}; do # run WINDEX on all testing files (sims 1-80 were used for training!)

    python hmm_scripts/hierarchical_hmm.py \
        ~/example/emission_training/site_stats/ \
        ~/example/emission_training/window_stats/ \
        example_$i.site_classified \
        example_$i.window_classified \
        ~/example/emission_training/site_stats/stat_files/YRI.300.0.05.$i.stats \
        ~/example/emission_training/window_stats/stat_files/YRI.300.0.05.$i.stats \
        40000; 

done
        
