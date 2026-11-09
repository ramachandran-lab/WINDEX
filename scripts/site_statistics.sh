#!/bin/bash

#SBATCH -J site_statistics	    
#SBATCH -p gpu --gres=gpu:1         
#SBATCH -N 1 						
#SBATCH -n 1                       
#SBATCH -t 2-00:00 			       
#SBATCH --mem 5G 				    

module load python
source ~/.venv/windex_env/bin/activate

### for 1000 Genomes Data 

pol_dir=polarity_files

for pop in 'YRI' 'CEU' 'CHB'; do
    for chrom in (1..22); do

        echo -e "array task id: $chrom || params = pop: $pop"

        vcf_file="threepops_chr${chrom}.recode.vcf"

        # set pops and refs for stat calculations
        if [ $pop == "YRI" ]; then 
            focal_pop="YRI"
            ref_pop="CEU"
            alt_pop="CHB"
        elif [ $pop == "CEU" ]; then
            focal_pop="CEU"
            ref_pop="CHB"
            alt_pop="YRI"
        elif [ $pop == "CHB" ]; then
            focal_pop="CHB"
            ref_pop="YRI"
            alt_pop="CEU"
        fi

        # perform stat calculations
        module load bcftools

        # filter for biallelic sites and MAF >= 0.01 
        if [[ ! -s ${vcf_file}.filtered ]]; then
            bcftools norm -m- ${vcf_file} -Ov | \
            bcftools +fill-tags -- -t AC,AN | \
            bcftools view -i 'AC/AN>=0.01' -Ov -o ${vcf_file}.filtered

            bcftools +fill-tags ${vcf_file} -Ou -- -t AC,AN | bcftools view -m2 -M2 -e 'AC==0 || AC==AN' -Ou | bcftools annotate -x INFO -Oz -o ${base}.filtered.vcf.gz

        fi

        # per-site ancestral allele polarity (assisted w Claude Opus v5)
        if [[ ! -s ${pol_dir}/chr${chrom}.alt_is_derived ]]; then
            bcftools query -f '%POS\t%REF\t%ALT\t%INFO/AA\n' ${vcf_file}.filtered \
              | awk -F'\t' -v c=$chrom 'BEGIN{OFS="\t"}
                  {aa = toupper(substr($4,1,1))
                   if      (aa == toupper($2)) { d = 1;  r++ }
                   else if (aa == toupper($3)) { d = 0;  a++ }
                   else                        { d = -1; u++ }
                   if (!($1 in seen)) { seen[$1]; print $1, d }}
                  END {printf("chr%-3s REF_anc=%-9d ALT_anc=%-8d unknown=%-8d frac_flipped=%.4f frac_unknown=%.4f\n",
                              c, r, a, u, a/(r+a+u), u/(r+a+u)) > "/dev/stderr"}' \
              > ${pol_dir}/chr${chrom}.alt_is_derived 2>> ${pol_dir}/summary.txt
        fi

        module unload bcftools

        source ~/venv/windex_env/bin/activate
        module load vcftools

        # calculate 
        ~/venv/windex_env/bin/python site_helper.py ${vcf_file}.filtered $focal_pop $ref_pop $alt_pop 1kgp_${pop}_${chrom}.site_stats ${pol_dir}/chr${chrom}.alt_is_derived
        
        # cleanup
        # rm ${vcf_file}.filtered
        # rm ${vcf_file}.filtered.dedup
        # rm ${vcf_file}.filtered.maf

    done
done


# ------------------ running selscan to obtain XP-EHH (need to run on CPUs) ------------------------- #

module load bcftools

focal_pop='YRI'

for i in (1..22); do 
    
    vcf_file="threepops_chr${i}.recode.vcf"

    if [[ ! -f "${vcf_file}.filtered" ]]; then

        # perform stat calculations
        # filter for biallelic sites and MAF >= 0.01 
        bcftools norm -m- ${vcf_file} -Ov | \
        bcftools +fill-tags -- -t AC,AN | \
        bcftools view -i 'AC/AN>=0.01' -Ov -o ${vcf_file}.filtered 

    fi

    # filter biallelic SNPs, polymorphic in the pooled sample
    bcftools +fill-tags ${vcf_file} -Ou -- -t AC,AN | bcftools view -m2 -M2 -e 'AC==0 || AC==AN' -Ou | bcftools annotate -x INFO -Oz -o ${base}.filtered.vcf.gz
    bcftools index ${base}.filtered.vcf.gz

    # split by population (same sites in every file)
    for pop in 'YRI' 'CEU' 'CHB'; do
        bcftools view -S pops_1kgp.${pop} ${base}.filtered.vcf.gz -Ov -o ${base}_${pop}.filtered.vcf
    done

    # one map file for all three (identical loci)
    bcftools view -H ${base}.filtered.vcf.gz | awk '{print $1"\t"$1"_"$2"\t"$2"\t"$2}' > ${base}.map

    # XP-EHH 
    echo 'XP-EHH'

    for pop in 'YRI' 'CEU' 'CHB'; do

        if [ $pop == $focal_pop ]; then 
            continue 
        else
            if [[ $focal_pop == "YRI" && $pop == "CEU" ]]; then
                # run xpehh 
                ./selscan --xpehh --vcf ${base}_${focal_pop}.filtered.vcf --vcf-ref ${base}_${pop}.filtered.vcf --map ${base}.map --maf 0.001 --threads 4 --trunc-ok --keep-low-freq --wagh --out ${base}_${focal_pop}

            elif [[ $focal_pop == "CEU" && $pop == "CHB" ]]; then
                ./selscan --xpehh --vcf ${base}_${focal_pop}.filtered.vcf --vcf-ref ${base}_${pop}.filtered.vcf --map ${base}.map --maf 0.001 --threads 4 --trunc-ok --keep-low-freq --wagh --out ${base}_${focal_pop}

            elif [[ $focal_pop == "CHB" && $pop == "YRI" ]]; then
                ./selscan --xpehh --vcf ${base}_${focal_pop}.filtered.vcf --vcf-ref ${base}_${pop}.filtered.vcf --map ${base}.map --maf 0.001 --threads 4 --trunc-ok --keep-low-freq --wagh --out ${base}_${focal_pop}

            fi
        fi

    done
done