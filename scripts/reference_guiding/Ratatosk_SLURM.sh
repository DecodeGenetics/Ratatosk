#!/bin/bash

MIN_NB_THREADS=""
MAX_NB_THREADS=""
MAX_NB_NODES=""
REF_GENOME=""
SHORT_READS_BAM=""
LONG_READS_BAM=""
SHORT_READS_PHASING=""
LONG_READS_PHASING=""
OUT_PREFIX=""
PARTITION="nomosix"

PRINT_HELP=0

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--min-cores) MIN_NB_THREADS="$2"; shift ;;
        -C|--max-cores) MAX_NB_THREADS="$2"; shift ;;
        -n|--max-nodes) MAX_NB_NODES="$2"; shift ;;
        -r|--reference) REF_GENOME="$2" ; shift ;;
	-s|--in-short-bam) SHORT_READS_BAM="$2" ; shift ;;
	-l|--in-long-bam) LONG_READS_BAM="$2" ; shift ;;
	-p|--in-short-phase) SHORT_READS_PHASING="$2" ; shift ;;
	-P|--in-long-phase) LONG_READS_PHASING="$2" ; shift ;;
	-o|--out-pref) OUT_PREFIX="$2" ; shift ;;
	-m|--slurm_part) PARTITION="$2" ; shift ;;
	-h|--help) PRINT_HELP=1 ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ ${PRINT_HELP} -eq 1 ]
then

	echo "Usage:"
	echo ""
	echo "bash $0 -c <MIN_NB_THREADS> -C <MAX_NB_THREADS> -n <MAX_NB_NODES> -r <REFRENCE_GENOME> -s <SHORT_READS_BAM> -l <LONG_READS_BAM> -o <OUTPUT_PATH> -p [SHORT_READS_PHASING] -P [LONG_READS_PHASING] -m [SLURM_PARTITION]"
	echo ""
	echo "<Mandatory>: "
	echo "-c <MIN_NB_THREADS>: Minimum number of threads that can be given to a job."
	echo "-C <MAX_NB_THREADS>: Maximum number of threads that can be given to a job."
	echo "-n <MAX_NB_NODES>: Maximum number of nodes that can be used at a given time."
	echo "-r <REFRENCE_GENOME>: Reference genome in FASTA"
	echo "-s <SHORT_READS_BAM>: Input BAM file of the short reads. Short reads from the same pair must have the same name."
	echo "-l <LONG_READS_BAM>: Input BAM file of the long reads."
	echo "-o <OUTPUT_PATH>: Output path."
	echo ""
	echo "[Optional]: "
	echo "-m [SLURM_PARTITION]: Slurm partition (default: nomosix)."
	echo "-p [SHORT_READS_PHASING]: Short reads phasing (default: None)"
	echo "-P [LONG_READS_PHASING]: Long reads phasing (default: None)"
	echo ""

	exit 0
fi

# Check min number of threads is provided and greater than 0
if [ "${MIN_NB_THREADS}" == "" ]
then
	1>&2 echo "Missing minimum number of threads for a job (-c). Abort."
	exit 1

elif [ "${MIN_NB_THREADS}" -eq 0 ]
then
	1>&2 echo "Minimum number of threads cannot be 0. Abort."
	exit 2
fi

# Check max number of threads is provided and greater than 0
if [ "${MAX_NB_THREADS}" == "" ]
then
	1>&2 echo "Missing maximum number of threads for a job (-C). Abort."
	exit 3

elif [ "${MAX_NB_THREADS}" -eq 0 ]
then
	1>&2 echo "Maximum number of threads cannot be 0. Abort."
	exit 4
fi

# Check min number of threads is not greater than max number of threads
if [ "${MIN_NB_THREADS}" -gt "${MAX_NB_THREADS}" ]
then
	1>&2 echo "Minimum number of threads cannot be greater than maximum number of threads. Abort."
	exit 5
fi

# Check short read BAM file is provided and exists
if [ "${SHORT_READS_BAM}" == "" ]
then
	1>&2 echo "Missing input short read BAMs (-s). Abort."
	exit 6

elif [ ! -f "${SHORT_READS_BAM}" ]
then
	1>&2 echo "Short read BAM file not found. Abort."
	exit 7
fi

# Check long read BAM file is provided and exists
if [ "${LONG_READS_BAM}" == "" ]
then
	1>&2 echo "Missing input long read BAMs (-l). Abort."
	exit 8

elif [ ! -f "${LONG_READS_BAM}" ]
then
	1>&2 echo "Long read BAM file not found. Abort."
	exit 9
fi

# Check reference genome file is provided and exists
if [ "${REF_GENOME}" == "" ]
then
	1>&2 echo "Missing input reference genome (-r). Abort."
	exit 10

elif [ ! -f "${REF_GENOME}" ]
then
	1>&2 echo "Reference genome file not found. Abort."
	exit 11
fi

# Check output prefix is provided
if [ "${OUT_PREFIX}" == "" ]
then
	1>&2 echo "Missing output path (-o). Abort."
	exit 12
fi

# Check short reads phasing exists if provided
if [ "${SHORT_READS_PHASING}" != "" ] && [ ! -f "${SHORT_READS_PHASING}" ]
then
	1>&2 echo "Short reads phasing file not found. Abort."
	exit 13

# Check long reads phasing exists if provided
elif [ "${LONG_READS_PHASING}" != "" ] && [ ! -f "${LONG_READS_PHASING}" ]
then
	1>&2 echo "Long reads phasing file not found. Abort."
	exit 14
fi

# Check short and long reads phasing are provided together
if [ "${SHORT_READS_PHASING}" != "" ] && [ "${LONG_READS_PHASING}" == "" ]
then
	1>&2 echo "Short reads phasing cannot be used without the long reads phasing. Abort."
	exit 15

elif [ "${SHORT_READS_PHASING}" == "" ] && [ "${LONG_READS_PHASING}" != "" ]
then
	1>&2 echo "Long reads phasing cannot be used without the short reads phasing. Abort."
	exit 16
fi

# Check number of nodes is not 0
if [ "${MAX_NB_NODES}" -eq 0 ]
then
	1>&2 echo "Maximum number of nodes cannot be 0. Abort."
	exit 17
fi

OUT_PREFIX=$(mkdir -p ${OUT_PREFIX}; cd ${OUT_PREFIX}; pwd)

echo "Minimum number of threads: ${MIN_NB_THREADS}"
echo "Maximum number of threads: ${MAX_NB_THREADS}"
echo "Input reference genome file: ${REF_GENOME}"
echo "Input short read file: ${SHORT_READS_BAM}"
echo "Input long read file: ${LONG_READS_BAM}"
echo "Output prefix: ${OUT_PREFIX}"
echo "Slurm partition: ${PARTITION}"

if [ "${SHORT_READS_PHASING}" != "" ]
then
	echo "Short reads phasing: ${SHORT_READS_PHASING}"
fi

if [ "${LONG_READS_PHASING}" != "" ]
then
	echo "Long reads phasing: ${LONG_READS_PHASING}"
fi

PREFIX_PATH_SEG="${OUT_PREFIX}/segments"
PREFIX_PATH_CORRECTED="${OUT_PREFIX}/ratatosk"
PREFIX_PATH_LOG="${OUT_PREFIX}/slurm_log"

PREFIX_FILE_BIN="${PREFIX_PATH_LOG}/binCorrectionarray.sbatch"
NAME_SR_UNMAPPED_IN_FILE="${PREFIX_PATH_SEG}/sample_sr_unmapped"

BIN_SZ=5000000

# Creating folders
mkdir -p "${OUT_PREFIX}"
mkdir -p "${PREFIX_PATH_SEG}"
mkdir -p "${PREFIX_PATH_CORRECTED}"
mkdir -p "${PREFIX_PATH_LOG}"

# Load chromosomes names and lengths in array
mapfile -t CHR_NAMES < <(samtools view -H "${SHORT_READS_BAM}" | grep "@SQ" | awk '{print substr($2,4,length($2))}')
mapfile -t CHR_LENGTHS < <(samtools view -H "${SHORT_READS_BAM}" | grep "@SQ" | awk '{print substr($3,4,length($3))}')

# Extract bins
CMD="python3 segmentBAM.py -t ${MAX_NB_THREADS} -s ${SHORT_READS_BAM} -l ${LONG_READS_BAM} -o ${PREFIX_PATH_SEG}/sample > ${PREFIX_PATH_SEG}/sample.bin;"

SLURM_OUT=$(sbatch -p ${PARTITION} -J Ratatosk_extractBins --mem=24G --cpus-per-task=${MAX_NB_THREADS} -t 2-0:0 -o ${PREFIX_PATH_LOG}/extractBins.slurm.out -e -${PREFIX_PATH_LOG}/extractBins.slurm.err --wrap="${CMD}")

# Extract job ID from SLURM output.
if ! echo ${SLURM_OUT} | grep -q "[1-9][0-9]*$"; then
	echo "Job(s) submission failed."
	echo ${SLURM_OUT}
	exit 18
else
	JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
fi

# Bin correction
echo "#!/bin/bash" > ${PREFIX_FILE_BIN}
echo "" >> ${PREFIX_FILE_BIN}
echo "#SBATCH --array=1-${MAX_NB_NODES}" >> ${PREFIX_FILE_BIN}
echo "#SBATCH --partition=${PARTITION}" >> ${PREFIX_FILE_BIN}
echo "#SBATCH --job-name=Ratatosk_correctSegONT1" >> ${PREFIX_FILE_BIN}
echo "#SBATCH --mem=16G" >> ${PREFIX_FILE_BIN}
echo "#SBATCH --time=1-0:0" >> ${PREFIX_FILE_BIN}
echo "#SBATCH --ntasks=${MIN_NB_THREADS}" >> ${PREFIX_FILE_BIN}
echo "#SBATCH --output=${PREFIX_PATH_LOG}/Ratatosk_correctSegONT1_%A_%a.out">> ${PREFIX_FILE_BIN}
echo "#SBATCH --error=${PREFIX_PATH_LOG}/Ratatosk_correctSegONT1_%A_%a.err" >> ${PREFIX_FILE_BIN}
echo "#SBATCH --dependency=afterok:${JOB_ID}" >> ${PREFIX_FILE_BIN}

# For all ONT reads with a good enough MAPQ, correct bin by bin
for j in $(seq 0 $((${#CHR_NAMES[@]})))
do
	chr_name="${CHR_NAMES[$j]}"
	chr_len="${CHR_LENGTHS[$j]}"

	# If chromosome name is not empty
	if [ ! -z ${chr_name} ]
	then

		for pos in $(seq 0 $((${BIN_SZ})) $((${chr_len})))
		do
			NAME_LR_IN_FILE="${PREFIX_PATH_SEG}/sample_lr_${chr_name}_${pos}"
			NAME_SR_IN_FILE="${PREFIX_PATH_SEG}/sample_sr_${chr_name}_${pos}"

			NAME_LR_OUT_FILE="${NAME_LR_IN_FILE}_corrected"

			RATATOSK_CMD="Ratatosk -v -c ${MIN_NB_THREADS} -s ${NAME_SR_IN_FILE}.fa -l ${NAME_LR_IN_FILE}.fq -u ${NAME_SR_UNMAPPED_IN_FILE}.fa -o ${NAME_LR_OUT_FILE}"

			if [ "${SHORT_READS_PHASING}" != "" ] && [ "${LONG_READS_PHASING}" != "" ]
			then
				RATATOSK_CMD="${RATATOSK_CMD} -p ${SHORT_READS_PHASING} -P ${LONG_READS_PHASING}"
			fi

			CMD="if [ -f ${NAME_LR_IN_FILE}.fq ] && [ -s ${NAME_LR_IN_FILE}.fq ]; then if [ -f ${NAME_SR_IN_FILE}.fa ] && [ -s ${NAME_SR_IN_FILE}.fa ]; then" # Check input file exists as they should
			CMD="${CMD} ${RATATOSK_CMD};"  # Ratatosk correction
			CMD="${CMD} else cp ${NAME_LR_IN_FILE}.fq ${NAME_LR_OUT_FILE}.fastq; fi; fi;"

			echo "${CMD}" >> ${PREFIX_FILE_BIN}
		done
	fi
done

SLURM_OUT=$(sbatch "${PREFIX_FILE_BIN}")

# Extract job ID from SLURM output.
if ! echo ${SLURM_OUT} | grep -q "[1-9][0-9]*$"; then
	echo "Job(s) submission failed."
	echo ${SLURM_OUT}
	exit 19
else
	JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
fi

CMD="cat \$(ls -t ${PREFIX_PATH_SEG}/sample_lr_*_corrected.fastq) > ${PREFIX_PATH_SEG}/sample_lr_map.fastq;" # Concat corrected bins
CMD="${CMD} samtools bam2fq -@ ${MAX_NB_THREADS} -n ${SHORT_READS_BAM} | gzip >${PREFIX_PATH_SEG}/sample_sr.fastq.gz;" # Extract all short reads
CMD="${CMD} Ratatosk -v -c ${MAX_NB_THREADS} -s ${PREFIX_PATH_SEG}/sample_sr.fastq.gz -l ${PREFIX_PATH_SEG}/sample_lr_unknown.fq -a ${PREFIX_PATH_SEG}/sample_lr_map.fastq -o ${PREFIX_PATH_SEG}/sample_lr_unknown_corrected;" # Ratatosk
CMD="${CMD} cat ${PREFIX_PATH_SEG}/sample_lr_map.fastq ${PREFIX_PATH_SEG}/sample_lr_unknown_corrected.fastq > ${PREFIX_PATH_CORRECTED}/sample_corrected.fastq;" # Concat corrected bins
CMD="${CMD} rm -rf ${PREFIX_PATH_SEG};"

SLURM_OUT=$(sbatch -p ${PARTITION} -J Ratatosk_correctSegONT2 --dependency=afterok:${JOB_ID} --mem=350G --cpus-per-task=${MAX_NB_THREADS} -t 7-0:0 -o ${PREFIX_PATH_LOG}/correctSegONT2.slurm.out -e ${PREFIX_PATH_LOG}/correctSegONT2.slurm.err --wrap="${CMD}")

# Extract job ID from SLURM output.
if ! echo ${SLURM_OUT} | grep -q "[1-9][0-9]*$"; then
	echo "Job(s) submission failed."
	echo ${SLURM_OUT}
	exit 20
else
	JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
fi
