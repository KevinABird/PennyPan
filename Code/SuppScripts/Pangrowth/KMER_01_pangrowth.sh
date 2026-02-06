#! /bin/bash

source ~/.bashrc
source 99_init_script_vars.sh

## create temp directory and echo path to .o file
TMP_DIR=`/bin/mktemp -d -p /mnt/data1/tmp`
cd ${TMP_DIR}
echo "${TMP_DIR}"

mamba activate ${MAMBA}/bioinfo_tools

# -----------------------------------------------------------------------------
# Subset sequences to genes, NLRs, collapse BED ranges such that each position 
# only represented once in each subset
# -----------------------------------------------------------------------------

OUTDIR="${WORKDIR}/pangrowth_output"
mkdir ${OUTDIR}
OUTFILE="${OUTDIR}/pangrowth_alpha_values.txt"
touch ${OUTFILE}

PGDIR="~/bin/pangrowth"

HV_NLRS="${INFO_DIR}/hvNLRsForKmer20250527.txt"

while read -a line
do 
	REF_FA=${line[0]}
	GENE_GFF=${line[1]}
	NLR_NAME=${line[2]}
	NAME=${NLR_NAME/Tarvensevar/}
	
	rsync -avuP ${REF_FA} .
	REF_FA=$(basename ${REF_FA})
	
	for c in ${chroms[@]}
	do
	
		rsync -avuP ${REF_DIR}/${NAME}.${c}.fa.gz .
		gunzip ${NAME}.${c}.fa.gz
		cat ${NAME}.${c}.fa >> ${NAME}_chr_only.fa
	
	done
	
	gzip -cd ${GENE_GFF} | \
	grep -v "^#" | \
	awk '{ if ($3 == "gene") print $1"\t"($4 - 1)"\t"$5 }' > ${NLR_NAME}_genes.bed
	sort -k1,1 -k2,2n ${NLR_NAME}_genes.bed > ${NLR_NAME}_genes.sorted.bed
	bedtools merge -i ${NLR_NAME}_genes.sorted.bed > ${NLR_NAME}_genes.merged.bed
	
	grep -v "^ID" ${INFO_DIR}/NLR_block_genes_20250430.txt | \
	grep "${NLR_NAME}" | awk '{ print $10"\t"($13 - 1)"\t"$14 }' > ${NLR_NAME}_NLR.bed
	sort -k1,1 -k2,2n ${NLR_NAME}_NLR.bed > ${NLR_NAME}_NLR.sorted.bed
	bedtools merge -i ${NLR_NAME}_NLR.sorted.bed > ${NLR_NAME}_NLR.merged.bed
	
	## get lists of specific hvNLRs of interest
	grep "${NLR_NAME}" ${HV_NLRS} | cut -f3 > ${NAME}_hvNLRs.txt
	gzip -cd ${GENE_GFF} | \
	grep -f ${NAME}_hvNLRs.txt | \
	awk '{ if ($3 == "gene") print $1"\t"($4 - 1)"\t"$5 }' > ${NAME}_hvNLRs.bed
	sort -k1,1 -k2,2n ${NAME}_hvNLRs.bed > ${NAME}_hvNLRs.sorted.bed
	bedtools merge -i ${NAME}_hvNLRs.sorted.bed > ${NAME}_hvNLRs.merged.bed
	
	seqtk subseq ${REF_FA} ${NLR_NAME}_genes.merged.bed > ${NLR_NAME}_genes.fa
	seqtk subseq ${REF_FA} ${NLR_NAME}_NLR.merged.bed > ${NLR_NAME}_NLR.fa
	seqtk subseq ${REF_FA} ${NAME}_hvNLRs.merged.bed > ${NAME}_hvNLRs.fa	

done < ${INFO_DIR}/nlr_and_gene_pangrowth_files.txt


# -----------------------------------------------------------------------------
# Run pangrowth for multiple k-mer lengths
# -----------------------------------------------------------------------------

for k in 25 31 37 43 49 55
do
	SET="pennycress_n7_k${k}"
	${PGDIR}/pangrowth hist -t 8 -k ${k} ./*_chr_only.fa > ${SET}_hist.txt
	${PGDIR}/pangrowth growth -h ${SET}_hist.txt > ${SET}_growth.txt
	ALPHA=$(python ${PGDIR}/scripts/plot_growth_AMH.py ${SET}_growth.txt ${SET}_growth.pdf)

	printf "%s\t%s\t%s\n" \
	"$k" "whole_genome" "$ALPHA" >> ${OUTFILE}
	
	SET="pennycress_n7_k${k}"
	${PGDIR}/pangrowth hist -t 8 -k ${k} ./*genes.fa > ${SET}_genes_hist.txt
	${PGDIR}/pangrowth growth -h ${SET}_genes_hist.txt > ${SET}_genes_growth.txt
	ALPHA=$(python ${PGDIR}/scripts/plot_growth_AMH.py ${SET}_genes_growth.txt ${SET}_genes_growth.pdf)
	
	printf "%s\t%s\t%s\n" \
	"$k" "genes" "$ALPHA" >> ${OUTFILE}
	
	SET="pennycress_n7_k${k}"
	${PGDIR}/pangrowth hist -t 8 -k ${k} ./*NLR.fa > ${SET}_NLR_hist.txt
	${PGDIR}/pangrowth growth -h ${SET}_NLR_hist.txt > ${SET}_NLR_growth.txt
	ALPHA=$(python ${PGDIR}/scripts/plot_growth_AMH.py ${SET}_NLR_growth.txt ${SET}_NLR_growth.pdf)
	
	printf "%s\t%s\t%s\n" \
	"$k" "NLR" "$ALPHA" >> ${OUTFILE}
	
	SET="pennycress_n7_k${k}"
	${PGDIR}/pangrowth hist -t 8 -k ${k} ./*_hvNLRs.fa > ${SET}_hvNLRs_hist.txt
	${PGDIR}/pangrowth growth -h ${SET}_hvNLRs_hist.txt > ${SET}_hvNLRs_growth.txt
	ALPHA=$(python ${PGDIR}/scripts/plot_growth_AMH.py ${SET}_hvNLRs_growth.txt ${SET}_hvNLRs_growth.pdf)
	
	printf "%s\t%s\t%s\n" \
	"$k" "hvNLRs" "$ALPHA" >> ${OUTFILE}
	
done

rsync -avuP *_chr_only.fa ${REF_DIR}
rsync -avuP *.txt ${OUTDIR}
rsync -avuP *.pdf ${OUTDIR}

# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------

rsync -avuP ${TMP_DIR}/* ${OUTDIR}

cd ${OUTDIR}
rm -rf ${TMP_DIR}

mamba deactivate
