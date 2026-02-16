#!/bin/bash
#$ -N Pennycress_polishing_lib_align
#$ -cwd
#$ -M jrifkin@hudsonalpha.org
#$ -m aes
#$ -V
#$ -q mecat.q
#$ -pe smp 64
#$ -p -100
#$ -o Pennycress_polishing_lib_align_$JOB_ID.out
#$ -e Pennycress_polishing_lib_align_$JOB_ID.err



source ~/.bashrc

conda activate bwa

# -----------------------------------------------------------------------------
# Align reads, and check stats
# -----------------------------------------------------------------------------
WORKDIR=/home/jrifkin/my_scratch/Pennycress/polishing_library_SNP_calls
IN_DIR=${WORKDIR}/raw_fq
OUT_DIR=${WORKDIR}/sorted_bam
STAT_DIR=${WORKDIR}/stats
REFDIR=${WORKDIR}

for SAMPLE in Tibet33 AK34W MN106ref MN134 Lovettomn PI650286 Ames32873


do
    echo $SAMPLE
    bwa-mem2 mem -t 64 -M \
            ${REFDIR}/Tarvensevar_MN106_872_v4.0.fa \
            ${IN_DIR}/${SAMPLE}*R1.fastq ${IN_DIR}/${SAMPLE}*R2.fastq > tmp_${SAMPLE}.sam 
    samtools sort -O BAM -@ 64 tmp_${SAMPLE}.sam > ${OUT_DIR}/sorted_${SAMPLE}.bam
    samtools flagstat ${OUT_DIR}/sorted_${SAMPLE}.bam > \
    ${STAT_DIR}/flagstat_sorted_${SAMPLE}.txt
done





#for SAMPLE in Tibet33 AK34W MN106ref MN134 Lovettomn PI650286 Ames32873
#do
#    echo $SAMPLE
#    nohup bzip2 -dc *${SAMPLE}*R1.fastq.bz2 > ${SAMPLE}_R1.fastq &
#    nohup bzip2 -dc *${SAMPLE}*R2.fastq.bz2 > ${SAMPLE}_R2.fastq &
#done 