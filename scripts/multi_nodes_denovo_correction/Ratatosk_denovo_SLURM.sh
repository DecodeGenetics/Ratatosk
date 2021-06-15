#!/bin/bash

MIN_NB_THREADS=""
MAX_NB_THREADS=""
MAX_NB_NODES=""
SHORT_READS_FQ=""
LONG_READS_FQ=""
OUT_PATH=""
PARTITION="nomosix"
NB_LINES_SEGMENT=1000000

PRINT_HELP=0

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--min-cores) MIN_NB_THREADS="$2"; shift ;;
        -C|--max-cores) MAX_NB_THREADS="$2"; shift ;;
        -n|--max-nodes) MAX_NB_NODES="$2"; shift ;;
	-s|--in-short-fq) SHORT_READS_FQ="$2" ; shift ;;
	-l|--in-long-fq) LONG_READS_FQ="$2" ; shift ;;
	-o|--out-pref) OUT_PATH="$2" ; shift ;;
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
	echo "bash $0 -c <MIN_NB_THREADS> -C <MAX_NB_THREADS> -n <MAX_NB_NODES> -s <SHORT_READS_FQ> -l <LONG_READS_FQ> -o <OUTPUT_PATH> -m [SLURM_PARTITION]"
	echo ""
	echo "<Mandatory>: "
	echo "-c <MIN_NB_THREADS>: Minimum number of threads that can be given to a job."
	echo "-C <MAX_NB_THREADS>: Maximum number of threads that can be given to a job."
	echo "-n <MAX_NB_NODES>: Maximum number of nodes that can be used at a given time."
	echo "-s <SHORT_READS_FQ>: Input FASTQ file of the short reads. Short reads from the same pair must have the same name."
	echo "-l <LONG_READS_FQ>: Input FASTQ file of the long reads."
	echo "-o <OUTPUT_PATH>: Output path."
	echo ""
	echo "[Optional]: "
	echo "-m [SLURM_PARTITION]: Slurm partition (default: nomosix)."
	echo ""

	exit 0
fi

# Check min number of threads is provided and greater than 0
if [ "${MIN_NB_THREADS}" == "" ]
then
	1>&2 echo "Missing minimum number of threads for a job (-c). Abort."
	exit 1

elif [ "${MIN_NB_THREADS}" -le 0 ]
then
	1>&2 echo "Minimum number of threads cannot be 0. Abort."
	exit 2
fi

# Check max number of threads is provided and greater than 0
if [ "${MAX_NB_THREADS}" == "" ]
then
	1>&2 echo "Missing maximum number of threads for a job (-C). Abort."
	exit 3

elif [ "${MAX_NB_THREADS}" -le 0 ]
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

# Check short read FASTQ file is provided and exists
if [ "${SHORT_READS_FQ}" == "" ]
then
	1>&2 echo "Missing input short read FASTQs (-s). Abort."
	exit 6

elif [ ! -f "${SHORT_READS_FQ}" ]
then
	1>&2 echo "Short read FASTQ file not found. Abort."
	exit 7
fi

# Check long read FASTQ file is provided and exists
if [ "${LONG_READS_FQ}" == "" ]
then
	1>&2 echo "Missing input long read FASTQs (-l). Abort."
	exit 8

elif [ ! -f "${LONG_READS_FQ}" ]
then
	1>&2 echo "Long read FASTQ file not found. Abort."
	exit 9
fi

# Check output prefix is provided
if [ "${OUT_PATH}" == "" ]
then
	1>&2 echo "Missing output path (-o). Abort."
	exit 12
fi

# Check number of nodes is not 0
if [ "${MAX_NB_NODES}" -le 0 ]
then
	1>&2 echo "Maximum number of nodes cannot be 0. Abort."
	exit 17
fi

OUT_PATH=$(mkdir -p ${OUT_PATH}; cd ${OUT_PATH}; pwd)

echo "Minimum number of threads: ${MIN_NB_THREADS}"
echo "Maximum number of threads: ${MAX_NB_THREADS}"
echo "Maximum number of nodes: ${MAX_NB_NODES}"
echo "Input short read file: ${SHORT_READS_FQ}"
echo "Input long read file: ${LONG_READS_FQ}"
echo "Output prefix: ${OUT_PATH}"
echo "Slurm partition: ${PARTITION}"

PREFIX_PATH_SEG="${OUT_PATH}/segments"
PREFIX_PATH_CORRECTED="${OUT_PATH}/ratatosk"
PREFIX_PATH_LOG="${OUT_PATH}/slurm_log"

PREFIX_FILE_BIN1="${PREFIX_PATH_LOG}/correctReads1.sbatch"
PREFIX_FILE_BIN2="${PREFIX_PATH_LOG}/correctReads2.sbatch"

# Creating folders
mkdir -p "${OUT_PATH}"
mkdir -p "${PREFIX_PATH_SEG}"
mkdir -p "${PREFIX_PATH_CORRECTED}"
mkdir -p "${PREFIX_PATH_LOG}"

# Compute number of lines in the long read file
NB_LINES_LR=$(zcat ${LONG_READS_FQ} | wc -l | awk '{print $1}')

# Submit first "large" job: Split input long reads into batches and computing correction index for 1st pass 
CMD="zcat ${LONG_READS_FQ} | split -l ${NB_LINES_SEGMENT} -a 6 -d - ${PREFIX_PATH_SEG}/sample.fastq.part_;" # Split ONT reads into chunks
CMD="${CMD} Ratatosk index -v -1 -c ${MAX_NB_THREADS} -s ${SHORT_READS_FQ} -l ${LONG_READS_FQ} -o ${PREFIX_PATH_SEG}/sample;" # Build Ratatosk index for 1st correction pass

SLURM_OUT=$(sbatch -p ${PARTITION} -J Ratatosk_buildIndex1 --mem=250G --cpus-per-task=${MAX_NB_THREADS} -t 2-0:0 -o ${PREFIX_PATH_LOG}/buildIndex1.slurm.out -e -${PREFIX_PATH_LOG}/buildIndex1.slurm.err --wrap="${CMD}")

# Extract job ID from SLURM output.
if ! echo ${SLURM_OUT} | grep -q "[1-9][0-9]*$"; then
	echo "Job(s) submission failed."
	echo ${SLURM_OUT}
	exit 18
else
	JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
fi

# Submit many "small" jobs (as a job array). Each job corrects a batch of long reads. 
echo "#!/bin/bash" > ${PREFIX_FILE_BIN1}
echo "" >> ${PREFIX_FILE_BIN1}
echo "#SBATCH --array=1-${MAX_NB_NODES}" >> ${PREFIX_FILE_BIN1}
echo "#SBATCH --partition=${PARTITION}" >> ${PREFIX_FILE_BIN1}
echo "#SBATCH --job-name=Ratatosk_correctReads1" >> ${PREFIX_FILE_BIN1}
echo "#SBATCH --mem=100G" >> ${PREFIX_FILE_BIN1}
echo "#SBATCH --time=1-0:0" >> ${PREFIX_FILE_BIN1}
echo "#SBATCH --ntasks=${MIN_NB_THREADS}" >> ${PREFIX_FILE_BIN1}
echo "#SBATCH --output=${PREFIX_PATH_LOG}/Ratatosk_correctReads1_%A_%a.out" >> ${PREFIX_FILE_BIN1}
echo "#SBATCH --error=${PREFIX_PATH_LOG}/Ratatosk_correctReads1_%A_%a.err" >> ${PREFIX_FILE_BIN1}
echo "#SBATCH --dependency=afterok:${JOB_ID}" >> ${PREFIX_FILE_BIN1}

i=0 # Counter for the batches to process

for j in {000000..999999} # Assumes no more than 1M batches of 1M lines have been created
do

	if [ $i -lt ${NB_LINES_LR} ]
	then
		echo "Ratatosk correct -v -1 -c ${MIN_NB_THREADS} -g ${PREFIX_PATH_SEG}/sample.index.k31.fasta -d ${PREFIX_PATH_SEG}/sample.index.k31.rtsk -l ${PREFIX_PATH_SEG}/sample.fastq.part_${j} -o ${PREFIX_PATH_SEG}/sample.corrected.part_${j}; rm -rf ${PREFIX_PATH_SEG}/sample.fastq.part_${j};" >> ${PREFIX_FILE_BIN1}
		i=$((i+NB_LINES_SEGMENT))
	else
		break
	fi
done

SLURM_OUT=$(sbatch "${PREFIX_FILE_BIN1}")

# Extract job ID from SLURM output.
if ! echo ${SLURM_OUT} | grep -q "[1-9][0-9]*$"; then
	echo "Job(s) submission failed."
	echo ${SLURM_OUT}
	exit 19
else
	JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
fi

# Submit second "large" job: List corrected batches and computing correction index for 2nd pass 
CMD="ls -lh ${PREFIX_PATH_SEG}/sample.corrected.part_*.2.fastq | awk '{print \$9}' > ${PREFIX_PATH_SEG}/parts.txt;" # List all 1st pass corrected batches
CMD="${CMD} Ratatosk index -v -2 -c ${MAX_NB_THREADS} -g ${PREFIX_PATH_SEG}/sample.index.k63.fasta -l ${PREFIX_PATH_SEG}/parts.txt -o ${PREFIX_PATH_SEG}/sample;" # Build Ratatosk index for 2nd correction pass

SLURM_OUT=$(sbatch -p ${PARTITION} -J Ratatosk_buildIndex2 --dependency=afterok:${JOB_ID} --mem=150G --cpus-per-task=${MAX_NB_THREADS} -t 2-0:0 -o ${PREFIX_PATH_LOG}/buildIndex2.slurm.out -e ${PREFIX_PATH_LOG}/buildIndex2.slurm.err --wrap="${CMD}")

# Extract job ID from SLURM output.
if ! echo ${SLURM_OUT} | grep -q "[1-9][0-9]*$"; then
	echo "Job(s) submission failed."
	echo ${SLURM_OUT}
	exit 20
else
	JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
fi

# Submit many "small" jobs (as a job array). Each job corrects a batch of corrected long reads. 
echo "#!/bin/bash" > ${PREFIX_FILE_BIN2}
echo "" >> ${PREFIX_FILE_BIN2}
echo "#SBATCH --array=1-${MAX_NB_NODES}" >> ${PREFIX_FILE_BIN2}
echo "#SBATCH --partition=${PARTITION}" >> ${PREFIX_FILE_BIN2}
echo "#SBATCH --job-name=Ratatosk_correctReads2" >> ${PREFIX_FILE_BIN2}
echo "#SBATCH --mem=100G" >> ${PREFIX_FILE_BIN2}
echo "#SBATCH --time=1-0:0" >> ${PREFIX_FILE_BIN2}
echo "#SBATCH --ntasks=${MIN_NB_THREADS}" >> ${PREFIX_FILE_BIN2}
echo "#SBATCH --output=${PREFIX_PATH_LOG}/Ratatosk_correctReads2_%A_%a.out">> ${PREFIX_FILE_BIN2}
echo "#SBATCH --error=${PREFIX_PATH_LOG}/Ratatosk_correctReads2_%A_%a.err" >> ${PREFIX_FILE_BIN2}
echo "#SBATCH --dependency=afterok:${JOB_ID}" >> ${PREFIX_FILE_BIN2}

i=0

for j in {000000..999999} # Assumes no more than 1M batches of 1M lines have been created
do

	if [ $i -lt ${NB_LINES_LR} ]
	then
		echo "Ratatosk correct -v -2 -c ${MIN_NB_THREADS} -g ${PREFIX_PATH_SEG}/sample.index.k63.fasta -d ${PREFIX_PATH_SEG}/sample.index.k63.rtsk -l ${PREFIX_PATH_SEG}/sample.corrected.part_${j}.2.fastq -o ${PREFIX_PATH_SEG}/sample.corrected.part_${j}; rm -rf ${PREFIX_PATH_SEG}/sample.corrected.part_${j}.2.fastq;" >> ${PREFIX_FILE_BIN2}
		i=$((i+NB_LINES_SEGMENT))
	else
		break
	fi
done

SLURM_OUT=$(sbatch "${PREFIX_FILE_BIN2}")

# Extract job ID from SLURM output.
if ! echo ${SLURM_OUT} | grep -q "[1-9][0-9]*$"; then
	echo "Job(s) submission failed."
	echo ${SLURM_OUT}
	exit 21
else
	JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
fi

# Final job: Patch all the corrected batches together, clean up tmp files
CMD="cat ${PREFIX_PATH_SEG}/sample.corrected.part_*.fastq > ${PREFIX_PATH_CORRECTED}/sample.corrected.fastq; rm -rf ${PREFIX_PATH_SEG};"

SLURM_OUT=$(sbatch -p ${PARTITION} -J Ratatosk_catCorrectedReads --dependency=afterok:${JOB_ID} --mem=2G --cpus-per-task=1 -t 1-0:0 -o ${PREFIX_PATH_LOG}/catCorrectedReads.slurm.out -e ${PREFIX_PATH_LOG}/catCorrectedReads.slurm.err --wrap="${CMD}")

# Extract job ID from SLURM output.
if ! echo ${SLURM_OUT} | grep -q "[1-9][0-9]*$"; then
	echo "Job(s) submission failed."
	echo ${SLURM_OUT}
	exit 22
else
	JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
fi
