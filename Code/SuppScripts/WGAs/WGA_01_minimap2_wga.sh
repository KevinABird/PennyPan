#!/bin/bash

source ~/.bashrc
source 99_init_script_vars.sh

# -----------------------------------------------------------------------------
# Run minimap2
# -----------------------------------------------------------------------------

REF_FA="Tarvensevar_MN106_872_v4.0.fa"
OUTDIR="${WORKDIR}/mm2_output"

rsync -avuP ${REF_DIR}/${REF_FA} .

mkdir ${OUTDIR}

## loop over genomes to be compared against MN106, subset FASTAs to chromosomes only, and 
## run minimap2
while read -a line
do

	rsync -avuP ${REF_DIR}/${line[1]} .

	### ---------------	
	mamba activate ${MAMBA}/bioinfo_tools
	## keep chromosomes only
	for c in ${chroms[@]}
	do
		echo ${c} > tmp.list
		seqtk subseq ${line[1]} tmp.list >> ${line[0]}_chr_only.fa
	done
	mamba deactivate
	### ---------------
	
	OUT_05PAF="MN106_x_${line[0]}_asm05-f_0.paf"
	COMP_FA="${line[0]}_chr_only.fa"
	
	### ---------------
	mamba activate ${MAMBA}/minimap2
	
	minimap2 -cx asm5 --eqx -t 32 -f 0 ${REF_FA} ${COMP_FA} > ${OUT_05PAF}
	
	mamba deactivate
	### ---------------	

done < ${INFO_DIR}/mm2_file_list.txt

mamba deactivate

rsync -avuP ./*.paf ${OUTDIR}