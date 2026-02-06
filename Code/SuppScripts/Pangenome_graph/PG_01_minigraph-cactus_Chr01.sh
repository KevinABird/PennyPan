#! /bin/bash

source ~/.bashrc
source 99_init_script_vars.sh

## create temp directory and echo path to .o file
TMP_DIR=`/bin/mktemp -d -p /mnt/data1/tmp`
cd ${TMP_DIR}
echo "${TMP_DIR}"


# -----------------------------------------------------------------------------
# Declare variables and prep files
# -----------------------------------------------------------------------------

readonly CACTUS_IMAGE=docker://quay.io/comparative-genomics-toolkit/cactus:v2.9.0

RUN_ID="${TAXON}_Chr01"

readonly JOBSTORE_IMAGE=jobStore_${RUN_ID}.img
readonly CACTUS_SCRATCH=${TMP_DIR}/cactus-${RUN_ID}
OUTDIR="${MC_OUT_DIR}/${RUN_ID}"
mkdir -p ${OUTDIR}

LOGFILE="${LOG_DIR}/01_minigraph-cactus_Chr01_$(date +"%Y_%m_%d_%I_%M_%p").log"
touch ${LOGFILE}

echo "Chr01 - prepping for construction - " $(date -u) >> ${LOGFILE}

# creating the local working directories
mkdir -p -m 777 ${CACTUS_SCRATCH}/upper ${CACTUS_SCRATCH}/work
truncate -s 300M "${JOBSTORE_IMAGE}"
apptainer exec ${CACTUS_IMAGE} mkfs.ext3 -d ${CACTUS_SCRATCH} "${JOBSTORE_IMAGE}"

mkdir -m 700 -p ${CACTUS_SCRATCH}/tmp
mkdir cactus_wd

# copy input files or cactus can't find them
cp ${SEQFILE} ./Chr01_seqfile.txt
sed -i 's/xCHROMx/Chr01/g' Chr01_seqfile.txt

cp ${CPFILE} ./Chr01_file_copy.txt
sed -i 's/xCHROMx/Chr01/g' Chr01_file_copy.txt

while read -a line
do
	cp ${line[0]} .
done < Chr01_file_copy.txt

# -----------------------------------------------------------------------------
# Run Minigraph-Cactus pangenome
# -----------------------------------------------------------------------------

echo "Chr01 - running minigraph-cactus graph construction - " $(date -u) >> ${LOGFILE}
apptainer exec --cleanenv \
  --overlay ${JOBSTORE_IMAGE} \
  --bind ${CACTUS_SCRATCH}/tmp:/tmp \
  --env PYTHONNOUSERSITE=1 \
  ${CACTUS_IMAGE} \
  cactus-pangenome \
  	--workDir=cactus_wd \
  	--binariesMode local \
  	--reference ${REFSAMP} \
  	--chop \
  	--vcf clip \
  	--giraffe clip \
  	--gfa clip \
  	--gbz full clip \
  	--xg clip \
  	--chrom-vg clip \
  	--chrom-og full clip \
  	--outDir ${RUN_ID} \
  	--outName ${RUN_ID} \
  	js \
  	Chr01_seqfile.txt

echo "Chr01 - graph construction complete - " $(date -u) >> ${LOGFILE}

# -----------------------------------------------------------------------------
# Viz with odgi
# -----------------------------------------------------------------------------

mamba activate ${MAMBA}/pg_tools

cut -f 1 Chr01_seqfile.txt > samp.list
sed -e 's/$/#/' -i samp.list

echo "Chr01 - starting odgi - " $(date -u) >> ${LOGFILE}
OG=$(find . -name "*.og")
for OGFN in ${OG}
do

	base=$(basename ${OGFN} .og)
	SORTOG="${base}.sorted.og"
	PNG="${TAXON}_${base}.sorted.png"
	INV_PNG="${TAXON}_${base}.sorted.inversions.png"
	COL_PNG="${TAXON}_${base}.sorted.collapsed.png"

	odgi sort \
		-i ${OGFN} \
		--threads ${MC_NTHREADS} \
		-P \
		-O \
		-Y \
		-o ${SORTOG}

	odgi viz \
		-i ${SORTOG} \
		--threads ${MC_NTHREADS} \
		-o ${PNG}
		
	odgi viz \
		-i ${SORTOG} \
		--threads ${MC_NTHREADS} \
		-o ${INV_PNG} \
		-M samp.list \
		-z
	
	odgi viz \
		-i ${SORTOG} \
		--threads ${MC_NTHREADS} \
		-o ${COL_PNG} \
		-M samp.list

done

echo "Chr01 - odgi complete - " $(date -u) >> ${LOGFILE}
mamba deactivate

# -----------------------------------------------------------------------------
# Clean up tmp dir
# -----------------------------------------------------------------------------
rm *.list
rm *.fa
rm *.fa.gz
rm *.img
rm Chr01_seqfile.txt
rm Chr01_file_copy.txt
rsync -avuP *.png ${ODGI_OUT}
rm *.png
rsync -avuP $TMP_DIR/* ${OUTDIR}

cd ${OUTDIR}
rm -rf ${TMP_DIR}