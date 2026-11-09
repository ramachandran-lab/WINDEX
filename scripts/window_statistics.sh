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

### for 1000 Genomes Data 

for pop in 'YRI' 'CEU' 'CHB'; do
    for chrom in {1..22}; do

        echo -e "params = pop: $pop | chrom: $chrom"

        vcf_file="threepops_chr${chrom}.recode.vcf"

        chr_size=`cat $vcf_file | tail -1 | awk '{print $2}'`

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
        bcftools view -m2 -M2 -i 'AC/AN>=0.01' -Ov -o ${vcf_file}.filtered ${vcf_file} 

        # get population-specific file
        for p in "YRI" "CEU" "CHB"; do
            bcftools view -S ${p}_1kgp.samples -o ${vcf_file}.filtered.${p} ${vcf_file}.filtered 
        done 

        # calculate derived frequency from VCF and output in lauren's file format (assisted w Claude Sonnet v5)
        if bcftools view -h "${vcf_file}.filtered.${pop}" | grep -q '^##INFO=<ID=AA,'; then
            query_fmt='%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA[\t%GT]\n'
        else
            query_fmt='%CHROM\t%POS\t%REF\t%ALT\t.[\t%GT]\n'
        fi

        bcftools query -f "$query_fmt" \
            "${vcf_file}.filtered.${pop}" \
        | awk -F'\t' '
            BEGIN { OFS="\t"; print "CHROM","POS","REF","ALT","AA","DERIVED_COUNT","TOTAL_ALLELES","DAF" }
            {
                chrom=$1; pos=$2; ref=$3; alt=$4
                aa_raw=$5

                # for 1kgp data specifically 
                sub(/\|.*/, "", aa_raw)

                if (aa_raw == "." || aa_raw == "" || aa_raw == "N") { # grab Ns too
                    aa = ref
                } else {
                    aa = toupper(aa_raw)
                }

                if (aa != ref && aa != alt) { print chrom,pos,ref,alt,aa,"NA","NA","NA"; next }

                total=0; derived=0
                for (col=6; col<=NF; col++) {
                    gt=$col
                    gsub(/\|/,"/",gt)
                    n=split(gt,alleles,"/")
                    for (k=1; k<=n; k++) {
                        a=alleles[k]
                        if (a=="0" || a=="1") {
                            total++
                            if ((aa==ref && a=="1") || (aa==alt && a=="0")) derived++ # for 1kgp 
                        }
                    }
                }
                daf = (total>0) ? derived/total : "NA"
                print chrom,pos,ref,alt,aa,derived,total,daf
            }' > "${vcf_file}.filtered.${pop}.frq_derived"

        awk 'BEGIN{OFS="\t"; print "SNP","CHROM","CHROM_POS","ALLELE1","FREQ1","ALLELE2","FREQ2"}
        NR>1{daf=$8; aa=$5; ref=$3; alt=$4;
            if(daf=="NA"){aaf="NA"} else {aaf=1-daf;}
            if(aa==ref){der=alt} else if(aa==alt){der=ref} else next;
            print $2, $1, $2, der, daf, aa, aaf}' \
        "${vcf_file}.filtered.${pop}.frq_derived" > "${vcf_file}.filtered.${pop}.ancestral_derived_freq"
        
        # setup env for stat calculation
        module unload bcftools
        source ~/venv/windex_env/bin/activate

        outputname="1kgp_${pop}_${chrom}.window_stats"
        ~/venv/windex_env/bin/python window_helper.py ${vcf_file}.filtered.${focal_pop} ${vcf_file}.filtered.${ref_pop} ${vcf_file}.filtered.${alt_pop} ${vcf_file}.filtered.${focal_pop}.ancestral_derived_freq $outputname $window_size $chr_size
        
        # cleanup
        # rm ${vcf_file}.filtered
    
    done
done
