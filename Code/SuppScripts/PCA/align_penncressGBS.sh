#!/bin/bash
#$ -N align_pennycress
#$ -M cmclaughlin@hudsonalpha.org
#$ -m besa
#$ -cwd
#$ -pe smp 32
#$ -q mecat.q
#$ -p -50
#$ -o align_pennycress.out
#$ -e align_pennycress.err

source ~/.bashrc
conda activate /home/jrifkin/.conda/envs/bwa

####################################################
###--Index reference, align reads, check scores--###
####################################################

# Refence genome
REF=refs/Tarvensevar_MN106_872_v4.0.fa

# Index reference genome for bwa aligner
bwa-mem2 index $REF

# Loop through all samples
for row in {A..H}; do
  for col in {01..12}; do

  SAMPLE="${row}${col}"

  echo ${SAMPLE}

  # Read names for demultiplexed paired-end data
  R1=trimmed-${SAMPLE}.1.fastq.gz
  R2=trimmed-${SAMPLE}.2.fastq.gz

  # Align the paired-end reads, make tmp SAM, and sorted BAM file
  bwa-mem2 mem -t 32 -M $REF trimmed/${R1} trimmed/${R2} > tmp_sam/tmp_${SAMPLE}.sam 
  samtools sort -O BAM -@ 32 tmp_sam/tmp_${SAMPLE}.sam > \
  sorted_bam/sorted_${SAMPLE}.bam
  
  samtools flagstat sorted_bam/sorted_${SAMPLE}.bam > \
  flagstat/flagstat_sorted_${SAMPLE}.txt
  
  done
done
