#!/bin/bash
#$ -N pennycress_polishing_pileup_and_calls
#$ -cwd
#$ -M jrifkin@hudsonalpha.org
#$ -m aes
#$ -V
#$ -q mecat.q
#$ -pe smp 48
#$ -p -100
#$ -o Pennycress_polishing_call_$JOB_ID.out
#$ -e Pennycress_polishing_call_$JOB_ID.err



source ~/.bashrc
# -----------------------------------------------------------------------------
# Activate env if necessary, set variable names
# -----------------------------------------------------------------------------

conda activate snp_calling

IN_DIR=/home/jrifkin/my_scratch/Pennycress/polishing_library_SNP_calls/sorted_bam
OUT_DIR=/home/jrifkin/my_scratch/Pennycress/polishing_library_SNP_calls/variant_calls
REF_FA=Tarvensevar_MN106_872_v4.0.fa

LOGFILE="pennycress_pileup_$(date +"%Y_%m_%d_%I_%M_%p").log"
touch ${LOGFILE}



# -----------------------------------------------------------------------------
# Filter and deduplicate sorted BAMs
# -----------------------------------------------------------------------------

for SAMPLE in AK34W Ames32873 Lovettomn MN106ref MN134 PI650286 Tibet33

do
    echo $SAMPLE

    picard FixMateInformation \
    I=${IN_DIR}/sorted_${SAMPLE}.bam \
    O=${IN_DIR}/sorted_${SAMPLE}_fixmate.bam \
    ADD_MATE_CIGAR=true

    picard AddOrReplaceReadGroups \
    I=${IN_DIR}/sorted_${SAMPLE}_fixmate.bam \
    O=${IN_DIR}/sorted_${SAMPLE}_fixmate_RG.bam \
    RGID=1 \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${SAMPLE}

    picard MarkDuplicates \
    I=${IN_DIR}/sorted_${SAMPLE}_fixmate_RG.bam \
    O=${IN_DIR}/ready_${SAMPLE}.bam \
    M=stats/${SAMPLE}_picard_stats.txt \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true \
    SORTING_COLLECTION_SIZE_RATIO=0.05 \
    REMOVE_DUPLICATES=true \
    TMP_DIR=${TMP_DIR} \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    MAX_RECORDS_IN_RAM=2750000 \
    READ_NAME_REGEX='[a-zA-Z0-9]+[:_][0-9]+[:_][a-zA-Z0-9]+[:_][0-9]+[:_]([0-9]+)[:_]([0-9]+)[:_]([0-9]+)' \
    CREATE_INDEX=true
done 


    #samtools fixmate --threads 48  -O bam ${IN_DIR}/sorted_${SAMPLE}.bam ${IN_DIR}/sorted_${SAMPLE}_fixmate.bam

    #samtools markdup  --threads 48  -O bam ${IN_DIR}/sorted_${SAMPLE}_fixmate.bam ${IN_DIR}/ready_${SAMPLE}.bam


# -----------------------------------------------------------------------------
# Call and filter SNPs
# -----------------------------------------------------------------------------

bcftools mpileup -a FMT/AD,FMT/DP -d 500 -q 20 -Q 20 --threads 48 -Ou -f ${REF_FA} ${IN_DIR}/ready_*.bam | bcftools call --threads 48 -mv -Ob -o ${OUT_DIR}/raw_calls_q20Q20.bcf 

## filter to keep sites with >= 8 reads
LO_LIM=8
filt="(FORMAT/DP)>=${LO_LIM}"
bcftools filter -Oz -i "$filt" -S . ${OUT_DIR}/raw_calls_q20Q20.bcf > ${OUT_DIR}/d8_set_to_missing_Q20.vcf.gz 

conda deactivate