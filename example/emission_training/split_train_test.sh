#!/bin/sh

### SLURM suggested params

#SBATCH -J generate_train_test_data        
#SBATCH -n 1    
#SBATCH --time 1:00:00                    				    
#SBATCH --mem 5G 

<< read_me
This script can be used to prepare your statistic files for input to training. 
This script splits the statistic files at an 80/20 ratio for training and testing, respectively. 
As an example, normalized statistics generated from 3POOA simulations in SLiM for YRI are used (80 training, 20 testing). 
It is not required to use SWIF(r) to generate the trained emissions, but WINDEX integrates with the SWIF(r) output to make it easier for users, so it is recommended. 
If you choose to use SWIF(r) for WINDEX emission training, please follow the SWIF(r) directory and file format available here: https://github.com/ramachandran-lab/SWIFr
read_me

module load python
source ~/venv/windex_env/bin/activate

proj_dir=~/example/emission_training
window_header=$'SNP_name\tPhysical_Distance\tMap_Distance\tTheta_Pi\tTheta_W\tTajima_D\tFay_Wu_H\tZeng_E\tGarud_H\tPBS\tNSS' # change these headers to match the statistics you are using
site_header=$'SNP_name\tPhysical_Distance\tMap_Distance\tXP-EHH\tDDAF\tnSL\tiHS\tFST'

### for window statistic files 

# make training files
for i in {1..80}; do
    tail -n+2 $proj_dir/window_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1==480001 {print}' | (echo "$window_header"; cat -) | sed 's/\s/\t/g' > $proj_dir/window_stats/simulations/sweep/YRI.300.0.05.sweep.train
for i in {1..80}; do 
    tail -n+2 $proj_dir/window_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1==440001 || $1==520001 {print}' | (echo "$window_header"; cat - ) | sed 's/\s/\t/g' > $proj_dir/window_stats/simulations/linked/YRI.300.0.05.linked.train
for i in {1..80}; do 
    tail -n+2 $proj_dir/window_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1!=440001 && $1!=520001 && $1!=480001 {print}' | (echo "$window_header"; cat - )| sed 's/\s/\t/g' > $proj_dir/window_stats/simulations/neutral/YRI.300.0.05.neutral.train

# make testing files
for i in {81..100}; do
    tail -n+2 $proj_dir/window_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1==480001 {print}' | (echo "$window_header"; cat - ) | sed 's/\s/\t/g' > $proj_dir/window_stats/testing_data/sweep/YRI.300.0.05.sweep.test
for i in {81..100}; do
    tail -n+2 $proj_dir/window_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1==440001 || $1==520001 {print}' | (echo "$window_header"; cat - )| sed 's/\s/\t/g' > $proj_dir/window_stats/testing_data/linked/YRI.300.0.05.linked.test
for i in {81..100}; do
    tail -n+2 $proj_dir/window_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1!=480001 && $1!=440001 && $1!=520001 {print}' | (echo "$window_header"; cat - ) | sed 's/\s/\t/g' > $proj_dir/window_stats/testing_data/neutral/YRI.300.0.05.neutral.test

### for site statistic files 

# make training files
for i in {1..80}; do
    tail -n+2 $proj_dir/site_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1==500001 {print}' | (echo "$site_header"; cat - ) | sed 's/\s/\t/g' > $proj_dir/site_stats/simulations/sweep/YRI.300.0.05.sweep.train
for i in {1..80}; do 
    tail -n+2 $proj_dir/site_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1>=440001 && $1<520001 && $1!=500001 {print}' | (echo "$site_header"; cat - ) | sed 's/\s/\t/g' > $proj_dir/site_stats/simulations/linked/YRI.300.0.05.linked.train
for i in {1..80}; do 
    tail -n+2 $proj_dir/site_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1<=440001 || $1>=520001 {print}' | (echo "$site_header"; cat - ) | sed 's/\s/\t/g' > $proj_dir/site_stats/simulations/neutral/YRI.300.0.05.neutral.train

# # make testing files
for i in {81..100}; do
    tail -n+2 $proj_dir/site_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1==480001 {print}' | (echo "$site_header"; cat - ) | sed 's/\s/\t/g' > $proj_dir/site_stats/testing_data/sweep/YRI.300.0.05.sweep.test
for i in {81..100}; do
    tail -n+2 $proj_dir/site_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1>=440001 && $1<520001 && $1!=500001 {print}' | (echo "$site_header"; cat - ) | sed 's/\s/\t/g' > $proj_dir/site_stats/testing_data/linked/YRI.300.0.05.linked.test
for i in {81..100}; do
    tail -n+2 $proj_dir/site_stats/stat_files/YRI.300.0.05.$i.stats; done | sort -k1,1n | awk '$1<=440001 || $1>=520001 {print}' | (echo "$site_header"; cat - ) | sed 's/\s/\t/g' > $proj_dir/site_stats/testing_data/neutral/YRI.300.0.05.neutral.test
