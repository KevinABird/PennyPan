#!/bin/bash
#$ -N callSNPs_unfiltered
#$ -M cmclaughlin@hudsonalpha.org
#$ -m besa
#$ -cwd
#$ -pe smp 32
#$ -q mecat.q
#$ -p -50
#$ -o callSNPs_unfiltered.out
#$ -e callSNPs_unfiltered.err

source ~/.bashrc
conda activate /home/jrifkin/.conda/envs/snp_calling/

##############################################################
###--Index alignment files. Call and filter variant sites--###
##############################################################

# Index all the sorted bam files
for bam in sorted_bam/*.bam; do
    samtools index "$bam"
done

# Generate raw VCF file with bcftools (Compute all genotypes | From all genotypes, call only the variants)
bcftools mpileup -a FMT/AD,FMT/DP -Ou -f refs/Tarvensevar_MN106_872_v4.0.fa -b bam_list.txt | bcftools call -mv -Oz -o pennycressGBS_unfiltered_raw.vcf.gz

# Rename the samples 
bcftools reheader -s samples.txt pennycressGBS_unfiltered_raw.vcf.gz -o pennycressGBS_unfiltered.vcf.gz