#!/bin/bash

#SBATCH -J 1kgp_site_stats     	    
#SBATCH -N 1 						
#SBATCH -n 1                        
#SBATCH -t 2-00:00 			        
#SBATCH --mem 25G 				    

module load python 
source ~/.venv/windex_env/bin/activate

for i in {1..22}; do

    for j in "yri" "ceu" "chb"; do
    echo $i, $j

    chrom=chr$i
    focal_pop=$j

    file_path=/path/to/vcfs
    pop_path=/path/to/population/definitions
    out_path=/output/path

#################### generate single-population VCF files ######################

# cat $file_path/ceu_yri_fin_chb_stu_${chrom}.vcf.gz | vcftools --gzvcf - --recode --recode-INFO-all --keep $pop_path/threepops_samples.txt --out $file_path/threepops_${chrom}
# cat $file_path/threepops_${chrom}.recode.vcf | bcftools view -i 'INFO/AC>0 & INFO/AC<INFO/AN' - | vcftools --vcf - --maf 0.05 --min-alleles 2 --max-alleles 2 --remove-indels --recode --recode-INFO-all --out $file_path/threepops_${chrom}_filtered
 
    for pop in "yri" "ceu" "chb"; do

        # split into single pop files
        cat $file_path/threepops_${chrom}_filtered.recode.vcf | vcftools --vcf - --recode --recode-INFO-all --keep $pop_path/${pop}_samples.txt --out $out_path/${pop}_filtered_${chrom}

        # make mapfiles from single-pop vcfs
        grep -v '#' $out_path/${pop}_filtered_${chrom}.recode.vcf | awk '{print $1"\t"$3"\t"$2"\t"$2 }' > $out_path/${pop}_filtered_${chrom}.map; done

    #################### calculate individual statistics ###########################

    # XP-EHH 
    echo 'XP-EHH'

    for pop in "yri" "ceu" "chb"; do

        if [ $pop == $focal_pop ]; then 
            continue 
        else
            
            if [[ $focal_pop == "yri" && $pop == "ceu" ]]; then
            # run xpehh 
            /users/hsnell/data_sramacha/hsnell/swifr/programs/selscan/bin/linux/./selscan --xpehh --vcf $out_path/${focal_pop}_filtered_${chrom}.recode.vcf --vcf-ref $out_path/${pop}_filtered_${chrom}.recode.vcf --map $out_path/${focal_pop}_filtered_${chrom}.map --maf 0.001 --threads 12 --trunc-ok --keep-low-freq --wagh --out $out_path/${focal_pop}_filtered_${chrom}

            elif [[ $focal_pop == "ceu" && $pop == "chb" ]]; then
            # run xpehh 
            /users/hsnell/data_sramacha/hsnell/swifr/programs/selscan/bin/linux/./selscan --xpehh --vcf $out_path/${focal_pop}_filtered_${chrom}.recode.vcf --vcf-ref $out_path/${pop}_filtered_${chrom}.recode.vcf --map $out_path/${focal_pop}_filtered_${chrom}.map --maf 0.001 --threads 12 --trunc-ok --keep-low-freq --wagh --out $out_path/${focal_pop}_filtered_${chrom}

            elif [[ $focal_pop == "chb" && $pop == "yri" ]]; then
            # run xpehh 
            /users/hsnell/data_sramacha/hsnell/swifr/programs/selscan/bin/linux/./selscan --xpehh --vcf $out_path/${focal_pop}_filtered_${chrom}.recode.vcf --vcf-ref $out_path/${pop}_filtered_${chrom}.recode.vcf --map $out_path/${focal_pop}_filtered_${chrom}.map --maf 0.001 --threads 12 --trunc-ok --keep-low-freq --wagh --out $out_path/${focal_pop}_filtered_${chrom}

            fi; 
        fi; 
    done

    # DDAF 
    echo 'DDAF'

    for pop in "yri" "ceu" "chb"; do

        # get derived allele frequencies for each population
        cat $out_path/${pop}_filtered_${chrom}.recode.vcf | awk -F'\t' 'BEGIN{print "CHROM\tPOS\tREF\tALT\tAA\tDERIVED_COUNT\tTOTAL_ALLELES\tDAF"} !/^#/ {ref=$4; alt=$5; if (match($8,/AA=([A-Za-z]+)/,a)) {aa=toupper(a[1]); total=0; derived=0; for(i=10;i<=NF;i++){if($i~/\.\/\./)continue; split($i,g,/[:\/|]/); for(j=1;j<=2;j++){if(g[j]~/[01]/){total++; if(aa==ref && g[j]==1)derived++; else if(aa==alt && g[j]==0)derived++;}}} if(total>0) daf=derived/total; else daf="NA"; print $1,$2,ref,alt,aa,derived,total,daf}}' OFS="\t" > $out_path/${pop}_filtered_${chrom}.frq_derived 
        cat $out_path/${pop}_filtered_${chrom}.frq_derived | awk '$5 != N {print}' | awk '{print $8}' | tail -n+2 > $out_path/${pop}_filtered_${chrom}.DAF; done

    cat $out_path/${focal_pop}_filtered_${chrom}.frq_derived | awk '$5 != N {print}' | awk '{print $2}' > $out_path/${focal_pop}_filtered_${chrom}.DAF.pos
    paste $out_path/${focal_pop}_filtered_${chrom}.DAF $out_path/yri_filtered_${chrom}.DAF $out_path/ceu_filtered_${chrom}.DAF | awk 'BEGIN { OFS = "\t" } {print ($1-(($2+$3)/2))}' | sed '1s/^/DDAF\n/' | paste $out_path/${focal_pop}_filtered_${chrom}.DAF.pos - > $out_path/${focal_pop}_filtered_${chrom}.DDAF.out

    # nSL
    echo 'nSL'

    /users/hsnell/data_sramacha/hsnell/swifr/programs/selscan/bin/linux/./selscan --nsl --vcf $out_path/${focal_pop}_filtered_${chrom}.recode.vcf --map  $out_path/${focal_pop}_filtered_${chrom}.map --maf 0.001 --threads 12 --trunc-ok --keep-low-freq --out $out_path/${focal_pop}_filtered_${chrom} 

    # iHS 
    echo 'iHS' 

    /users/hsnell/data_sramacha/hsnell/swifr/programs/selscan/bin/linux/./selscan --ihs --vcf $out_path/${focal_pop}_filtered_${chrom}.recode.vcf --map  $out_path/${focal_pop}_filtered_${chrom}.map --maf 0.001 --threads 12 --trunc-ok --max-extend -1 --out $out_path/${focal_pop}_filtered_${chrom}


    # Fst
    echo 'Fst'

    for pop in "yri" "ceu" "chb"; do

        grep "^#CHROM" $out_path/${pop}_filtered_${chrom}.recode.vcf | sed 's/\t/\n/g' | tail -n+10 > $out_path/${pop}_filtered_${chrom}.txt; done

    if [ $pop == "yri" ]; then
        vcftools --vcf $file_path/threepops_${chrom}_filtered.recode.vcf --weir-fst-pop $out_path/${focal_pop}_filtered_${chrom}.txt --weir-fst-pop $out_path/ceu_filtered_${chrom}.txt --weir-fst-pop $out_path/chb_filtered_${chrom}.txt --stdout | cut -f2,3 | sed "s/POS/SNP_name/g" | sed "s/WEIR_AND_COCKERHAM_FST/FST/g" > $out_path/${focal_pop}_filtered_${chrom}.fst.out

    elif [ $pop == "ceu" ]; then
        vcftools --vcf $file_path/threepops_${chrom}_filtered.recode.vcf --weir-fst-pop $out_path/${focal_pop}_filtered_${chrom}.txt --weir-fst-pop $out_path/chb_filtered_${chrom}.txt --weir-fst-pop $out_path/yri_filtered_${chrom}.txt --stdout | cut -f2,3 | sed "s/POS/SNP_name/g" | sed "s/WEIR_AND_COCKERHAM_FST/FST/g" > $out_path/${focal_pop}_filtered_${chrom}.fst.out

    elif [ $pop == "chb" ]; then
        vcftools --vcf $file_path/threepops_${chrom}_filtered.recode.vcf --weir-fst-pop $out_path/${focal_pop}_filtered_${chrom}.txt --weir-fst-pop $out_path/yri_filtered_${chrom}.txt --weir-fst-pop $out_path/ceu_filtered_${chrom}.txt --stdout | cut -f2,3 | sed "s/POS/SNP_name/g" | sed "s/WEIR_AND_COCKERHAM_FST/FST/g" > $out_path/${focal_pop}_filtered_${chrom}.fst.out

    fi; 

    done

done