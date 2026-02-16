#!/bin/bash
#$ -N filter_pennycress
#$ -M cmclaughlin@hudsonalpha.org
#$ -m besa
#$ -cwd
#$ -pe smp 32
#$ -q mecat.q
#$ -p -50
#$ -o filter_pennycress.out
#$ -e filter_pennycress.err

source ~/.bashrc
conda activate vcftools

###----------Relaxed filter for overlap with WGS libraries----------###
#Result is 637,506 SNPs 
VCF_IN=./pennycressGBS_unfiltered.vcf.gz
VCF_OUT=./pennycressGBS_filtered_20241125.vcf.gz

# set filters
MISS=0.70
QUAL=30
MIN_DEPTH=5
MAX_DEPTH=50

# filtering with vcftools
vcftools --gzvcf $VCF_IN \
--remove-indels --max-missing $MISS --minQ $QUAL \
--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
$VCF_OUT

###----------Indels removed, stricter filter----------###
#Result is 36,770 SNPs  
#VCF_IN=./pennycressGBS_unfiltered.vcf.gz
#VCF_OUT=./pennycressGBS_filtered_indelsremoved_20241121.vcf.gz

# set filters
#MISS=0.90
#QUAL=30
#MIN_DEPTH=10
#MAX_DEPTH=50

# filtering with vcftools
#vcftools --gzvcf $VCF_IN \
#--remove-indels --max-missing $MISS --minQ $QUAL \
#--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
#$VCF_OUT

###----------Relax filter for GWAS----------###
#Result is 228,994 sites
#VCF_IN=./pennycressGBS_unfiltered.vcf.gz
#VCF_OUT=./pennycressGBS_filtered_20241121.vcf.gz

# set filters
#MISS=0.80
#QUAL=30
#MIN_DEPTH=10
#MAX_DEPTH=50

# filtering with vcftools
#vcftools --gzvcf $VCF_IN \
#--max-missing $MISS --minQ $QUAL \
#--minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > \
#$VCF_OUT


