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

# use hierarchical_hmm.py -h

# run WINDEX

for i in {81..100}; do # run WINDEX on all testing files (sims 1-80 were used for training!)

    python ~/scripts/hierarchical_hmm.py \
        --path2trained_sites ~/example/emission_training/site_stats/ \
        --path2trained_windows ~/example/emission_training/window_stats/ \
        --sites_out example_$i.site_classified \
        --windows_out example_$i.window_classified \
        --datafile_sites ~/example/emission_training/site_stats/stat_files/YRI.300.0.05.$i.stats \
        --datafile_windows ~/example/emission_training/window_stats/stat_files/YRI.300.0.05.$i.stats \
        --window_size 40000; 

done
        
