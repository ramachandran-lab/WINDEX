#!/bin/sh

### SLURM suggested params

#SBATCH -J retrain_distributions_w_swifr	        
#SBATCH -n 4                        
#SBATCH -t 1-0 					    
#SBATCH --mem 50G 

<< read_me
This script can be used to train the emissions for WINDEX using its predecessor tool SWIF(r). 
SWIF(r) forms Gaussian mixture models for each input class, and model components are determined by BIC. 
For simplicity in this example, we reduce all GMM components to 1 for all input statistics after the first training step, then retrain. 
Optional code is provided to run testing data using NB-SWIF(r) as a comparison to WINDEX.
This training should be run for your choice of site and window statistics separately.
read_me

module load python # load python 3.13 into your env
source ~/venv/windex_env/bin/activate # activate windex python env

swifr_path=~/example/emission_training/site_stats # path to site or window statistic dir (needs to follow SWIF(r) dir setup)

# train SWIF(r) to obtain marginal statistic distributions
swifr_train --path $swifr_path

# optional, after examining distributions, retrain if Gaussian components need to be reduced
#swifr_train --path $swifr_path --retrain

# optional, but NB-SWIF(r) test comparison with same priors as used for WINDEX (priors ordered N, S, L)
#swifr_test --path2trained $swifr_path --nb --pi 0.88 0.04 0.08 --file $swifr_path/testing_data/sweep/YRI.300.0.05.test # testing sweep class
#swifr_test --path2trained $swifr_path --nb --pi 0.88 0.04 0.08 --file $swifr_path/testing_data/linked/YRI.300.0.05.test # testing linked class
#swifr_test --path2trained $swifr_path --nb --pi 0.88 0.04 0.08 --file $swifr_path/testing_data/neutral/YRI.300.0.05.test # testing neutral class
