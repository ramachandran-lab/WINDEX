#!/bin/bash

#SBATCH -J 1kgp_window_stats     
#SBATCH -N 1 						
#SBATCH -n 4                        
#SBATCH -t 1-00 			        
#SBATCH --mem 10G 				    

module load python
source ~/.venv/windex_env/bin/activate

windowsize=100000

vcf_path=/path/to/vcfs
vcf_pop_path=/path/to/population/definitions
out_path=/out/path
script_path=/script/location

for i in {1..22}; do
    for pop in "yri" "ceu" "chb"; do

    outputname=$out_path/$pop/${pop}_filtered_chr${i}.windowstats
    
    # convert vcf to freq file for lauren's stats (accounting for ancestral vs derived alleles)
    cat $vcf_pop_path/$pop/${pop}_filtered_chr${i}.recode.vcf | awk -F'\t' 'BEGIN{print "CHROM\tPOS\tREF\tALT\tAA\tDERIVED_COUNT\tTOTAL_ALLELES\tDAF"} !/^#/ {ref=$4; alt=$5; if (match($8,/AA=([A-Za-z]+)/,a)) {aa=toupper(a[1]); total=0; derived=0; for(i=10;i<=NF;i++){if($i~/\.\/\./)continue; split($i,g,/[:\/|]/); for(j=1;j<=2;j++){if(g[j]~/[01]/){total++; if(aa==ref && g[j]==1)derived++; else if(aa==alt && g[j]==0)derived++;}}} if(total>0) daf=derived/total; else daf="NA"; print $1,$2,ref,alt,aa,derived,total,daf}}' OFS="\t" > $out_path/$pop/${pop}_filtered_${i}.frq_derived
    awk 'BEGIN{OFS="\t"; print "SNP","CHROM","CHROM_POS","ALLELE1","FREQ1","ALLELE2","FREQ2"} NR>1{daf=$8; aa=$5; ref=$3; alt=$4; if(daf=="NA"){aaf="NA"} else {aaf=1-daf;} if(aa==ref){der=alt} else if(aa==alt){der=ref} else next; print $2, $1, $2, aa, aaf, der, daf}' $out_path/$pop/${pop}_filtered_${i}.frq_derived > $out_path/$pop/${pop}_filtered_${i}.ancestral_derived_freq

    # run stats
    if [ $pop == "yri" ]; then

    chr_size=`cat $vcf_pop_path/$pop/yri_filtered_chr${i}.recode.vcf | tail -1 | awk '{print $2}'`
    python $script_path/window_helper.py $vcf_pop_path/$pop/yri_filtered_chr${i}.recode.vcf $vcf_pop_path/ceu/ceu_filtered_chr${i}.recode.vcf $vcf_pop_path/chb/chb_filtered_chr${i}.recode.vcf $out_path/$pop/${pop}_filtered_${i}.ancestral_derived_freq $outputname $windowsize $chr_size

    elif [ $pop == "ceu" ]; then

    chr_size=`cat $vcf_pop_path/$pop/yri_filtered_chr${i}.recode.vcf | tail -1 | awk '{print $2}'`
    python $script_path/window_helper.py $vcf_pop_path/$pop/ceu_filtered_chr${i}.recode.vcf $vcf_pop_path/chb/chb_filtered_chr${i}.recode.vcf $vcf_pop_path/yri/yri_filtered_chr${i}.recode.vcf $out_path/$pop/${pop}_filtered_${i}.ancestral_derived_freq $outputname $windowsize $chr_size

    else 
    
    chr_size=`cat $vcf_pop_path/$pop/yri_filtered_chr${i}.recode.vcf | tail -1 | awk '{print $2}'`
    python $script_path/window_helper.py $vcf_pop_path/$pop/chb_filtered_chr${i}.recode.vcf $vcf_pop_path/yri/yri_filtered_chr${i}.recode.vcf $vcf_pop_path/ceu/ceu_filtered_chr${i}.recode.vcf $out_path/$pop/${pop}_filtered_${i}.ancestral_derived_freq $outputname $windowsize $chr_size

    fi

    done
done

