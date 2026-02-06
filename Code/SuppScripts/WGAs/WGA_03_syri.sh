#!/bin/bash

source ~/.bashrc
source 99_init_script_vars.sh

TMP_DIR=`/bin/mktemp -d -p /mnt/data1/tmp`
cd ${TMP_DIR}
echo "${TMP_DIR}"


# -----------------------------------------------------------------------------
# Run SyRI on minimap2 output (.paf format)
# -----------------------------------------------------------------------------

mamba activate ${MAMBA}/syri

REF_FA="Tarvensevar_MN106_872_v4.0.fa"
rsync -avuP ${REF_DIR}/${REF_FA} .

INDIR="${WORKDIR}/filtered_mm2_output"
OUTDIR="${WORKDIR}/syri_output"

mkdir ${OUTDIR}

while read -a line
do

	SAMP="${line[0]}"
	COMP_FA="${line[1]}"
	IN_PAF="filtered_MN106_x_${SAMP}_asm05-f_0_1000_0.8_2.paf"
	
	rsync -avuP ${REF_DIR}/${COMP_FA} .
	rsync -avuP ${INDIR}/${IN_PAF} .

	syri \
		-c ${IN_PAF} \
		-r ${REF_FA} \
		-q ${COMP_FA} \
		-k -F P \
		--prefix ${SAMP}_asm05-f_0_1000_0.8_2.
		
	rm ${IN_PAF}

done < ${INFO_DIR}/mm2_file_list.txt

mamba deactivate


# -----------------------------------------------------------------------------
# Clean up temp dir
# -----------------------------------------------------------------------------

rsync -avuP *asm05-f_0_1000_0.8_2* ${OUTDIR}
cd ${OUTDIR}
rm -rf ${TMP_DIR}