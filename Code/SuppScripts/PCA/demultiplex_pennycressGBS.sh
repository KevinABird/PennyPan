#!/bin/bash
#$ -N demultiplex_pennycress
#$ -M cmclaughlin@hudsonalpha.org
#$ -m besa
#$ -cwd
#$ -pe smp 48
#$ -q mecat.q
#$ -p -50
#$ -S /bin/bash
#$ -o demultiplex_pennycress.out
#$ -e demultiplex_pennycress.err

##################################################
###----------Demultiplex the GBS data----------###
##################################################

/home/jrifkin/.conda/envs/CutAdapt/bin/cutadapt -e 1 -g ^file:barcodes.fasta \
-o ./trimmed/"trimmed-{name}.1.fastq.gz" -p ./trimmed/"trimmed-{name}.2.fastq.gz" \
 /home/jrifkin/my_scratch/Pennycress/GBS/FP-GBS-96_S8_L004_R1_001.fastq-002.gz \
 /home/jrifkin/my_scratch/Pennycress/GBS/FP-GBS-96_S8_L004_R2_001.fastq-001.gz
