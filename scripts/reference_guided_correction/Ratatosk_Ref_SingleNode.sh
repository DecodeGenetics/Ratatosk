#!/bin/bash

NB_THREADS=1
REF_GENOME=""
SHORT_READS_BAM=""
LONG_READS_BAM=""
SHORT_READS_PHASING=""
LONG_READS_PHASING=""
OUT_PREFIX=""

PRINT_HELP=0

MAX_NB_CORES=$(nproc --all)

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--cores) NB_THREADS="$2"; shift ;;
        -r|--reference) REF_GENOME="$2" ; shift ;;
		-s|--in-short-bam) SHORT_READS_BAM="$2" ; shift ;;
		-l|--in-long-bam) LONG_READS_BAM="$2" ; shift ;;
		-p|--in-short-phase) SHORT_READS_PHASING="$2" ; shift ;;
		-P|--in-long-phase) LONG_READS_PHASING="$2" ; shift ;;
		-o|--out-pref) OUT_PREFIX="$2" ; shift ;;
		-h|--help) PRINT_HELP=1 ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ ${PRINT_HELP} -eq 1 ]
then

	echo "Usage:"
	echo ""
	echo "bash $0 -r <REFRENCE_GENOME> -s <SHORT_READS_BAM> -l <LONG_READS_BAM> -o <OUTPUT_PATH> -c [NB_THREADS] -p [SHORT_READS_PHASING] -P [LONG_READS_PHASING]"
	echo ""
	echo "<Mandatory>: "
	echo "-r <REFRENCE_GENOME>: Reference genome in FASTA"
	echo "-s <SHORT_READS_BAM>: Input BAM file of the short reads. Short reads from the same pair must have the same name."
	echo "-l <LONG_READS_BAM>: Input BAM file of the long reads."
	echo "-o <OUTPUT_PATH>: Output path."
	echo ""
	echo "[Optional]: "
	echo "-c [NB_THREADS]: Number of threads to use (default: 1)"
	echo "-p [SHORT_READS_PHASING]: Short reads phasing (default: None)"
	echo "-P [LONG_READS_PHASING]: Long reads phasing (default: None)"
	echo ""

	exit 0
fi

# Check if short read BAM file is provided and exists
if [ "${SHORT_READS_BAM}" == "" ]
then
	1>&2 echo "Missing input short read BAMs (-s). Abort."
	exit 1

elif [ ! -f "${SHORT_READS_BAM}" ]
then
	1>&2 echo "Short read BAM file not found. Abort."
	exit 2
fi

# Check long read BAM file is provided and exists
if [ "${LONG_READS_BAM}" == "" ]
then
	1>&2 echo "Missing input long read BAMs (-l). Abort."
	exit 3

elif [ ! -f "${LONG_READS_BAM}" ]
then
	1>&2 echo "Long read BAM file not found. Abort."
	exit 4
fi

# Check reference genome file is provided and exists
if [ "${REF_GENOME}" == "" ]
then
	1>&2 echo "Missing input reference genome (-r). Abort."
	exit 5

elif [ ! -f "${REF_GENOME}" ]
then
	1>&2 echo "Reference genome file not found. Abort."
	exit 6
fi

# Check number of threads doesn't exceed number of cores on machine
if [ "${NB_THREADS}" -gt "${MAX_NB_CORES}" ]
then
	1>&2 echo "Number of threads required exceed what is available. Abort."
	exit 7
fi

# Check output prefix is provided
if [ "${OUT_PREFIX}" == "" ]
then
	1>&2 echo "Missing output path (-o). Abort."
	exit 8
fi

# Check short reads phasing exists if provided
if [ "${SHORT_READS_PHASING}" != "" ] && [ ! -f "${SHORT_READS_PHASING}" ]
then
	1>&2 echo "Short reads phasing file not found. Abort."
	exit 9

# Check long reads phasing exists if provided
elif [ "${LONG_READS_PHASING}" != "" ] && [ ! -f "${LONG_READS_PHASING}" ]
then
	1>&2 echo "Long reads phasing file not found. Abort."
	exit 10
fi

# Check short and long reads phasing are provided together
if [ "${SHORT_READS_PHASING}" != "" ] && [ "${LONG_READS_PHASING}" == "" ]
then
	1>&2 echo "Short reads phasing cannot be used without the long reads phasing. Abort."
	exit 11

elif [ "${SHORT_READS_PHASING}" == "" ] && [ "${LONG_READS_PHASING}" != "" ]
then
	1>&2 echo "Long reads phasing cannot be used without the short reads phasing. Abort."
	exit 12
fi

OUT_PREFIX=$(mkdir -p ${OUT_PREFIX}; cd ${OUT_PREFIX}; pwd)

echo "Number of threads: ${NB_THREADS}"
echo "Input reference genome file: ${REF_GENOME}"
echo "Input short read file: ${SHORT_READS_BAM}"
echo "Input long read file: ${LONG_READS_BAM}"
echo "Output path: ${OUT_PREFIX}"

if [ "${SHORT_READS_PHASING}" != "" ]
then
	echo "Short reads phasing: ${SHORT_READS_PHASING}"
fi

if [ "${LONG_READS_PHASING}" != "" ]
then
	echo "Long reads phasing: ${LONG_READS_PHASING}"
fi

RatatoskBin() {

	if [ "${1}" != "unmapped_sr" ] && [ "${1}" != "unknown_lr" ] # Bin is not unmapped short reads or unmapped/low mapq long reads
	then
		if [ ! "${4}" -eq 0 ] && [ ! "${5}" -eq 0 ] # The bin contains at least one long and one short read to perform the correction
		then
			if [ "${9}" != "" ] && [ "${10}" != "" ] # Phasing is provided
			then
				if [ -f "${6}" ] # Unmapped short reads available to assist correction
				then
					Ratatosk correct -c "${8}" -s "${2}" -l "${3}" -u "${6}" -p "${9}" -P "${10}" -o "${7}/sample_lr_${1}_corrected"
				else
					Ratatosk correct -c "${8}" -s "${2}" -l "${3}" -p "${9}" -P "${10}" -o "${7}/sample_lr_${1}_corrected"
				fi
			else
				if [ -f "${6}" ] # Unmapped short reads available to assist correction
				then
					Ratatosk correct -c "${8}" -s "${2}" -l "${3}" -u "${6}" -o "${7}/sample_lr_${1}_corrected"
				else
					Ratatosk correct -c "${8}" -s "${2}" -l "${3}" -o "${7}/sample_lr_${1}_corrected"
				fi
			fi
		else
			cp "${3}" "${7}/sample_lr_${1}_corrected.fastq" # Copy uncorrected long reads to bin output
		fi

		if [ ! -f "${7}/sample_lr_${1}_corrected.fastq" ] # If no bin output, something went wrong with the correction. Return error code.
		then

			return 1
		fi
	fi

	return 0
}

PREFIX_PATH_SEG="${OUT_PREFIX}/segments"
PREFIX_PATH_CORRECTED="${OUT_PREFIX}/ratatosk"
NAME_SR_UNMAPPED_IN_FILE="${PREFIX_PATH_SEG}/sample_sr_unmapped.fa"

# Creating folders
mkdir -p "${OUT_PREFIX}"
mkdir -p "${PREFIX_PATH_SEG}"
mkdir -p "${PREFIX_PATH_CORRECTED}"

export -f RatatoskBin

# Extract bins
echo "- Creating bins of long and short reads"

python3 segmentBAM.py -t ${NB_THREADS} -s ${SHORT_READS_BAM} -l ${LONG_READS_BAM} -o "${PREFIX_PATH_SEG}/sample" > "${PREFIX_PATH_SEG}/sample.bins"

if [ $? -eq 0 ] && [ -f "${PREFIX_PATH_SEG}/sample.bins" ]
then

	# Correct each bin in parallel
	echo "- Correcting bins"

	#NB_BINS=$(wc -l "${PREFIX_PATH_SEG}/sample.bins" | cut -d ' ' -f1) # Get number of bins
	#NB_BINS=$((NB_BINS-2)) # Remove 2 bins from count (low mapq and unmapped short/long reads)
	#NB_THREADS_PER_JOB=$((NB_THREADS/NB_BINS)) # Get min number of threads per bin
	#if [ "${NB_THREADS_PER_JOB}" -eq 0 ]; then NB_THREADS_PER_JOB=1; fi
	#NB_JOBS=$((NB_THREADS/NB_THREADS_PER_JOB))
	#if [ "${NB_JOBS}" -eq 0 ]; then NB_JOBS=1; fi
	#parallel -j ${NB_JOBS} -C '\t' RatatoskBin {1} {2} {3} {4} {5} ${NAME_SR_UNMAPPED_IN_FILE} ${PREFIX_PATH_SEG} ${NB_THREADS_PER_JOB} :::: "${PREFIX_PATH_SEG}/sample.bins" && echo OK

	while read BIN
	do
		BIN_NAME=$(echo ${BIN} | awk '{print $1}')
		BIN_SR_FILEN=$(echo ${BIN} | awk '{print $2}')
		BIN_LR_FILEN=$(echo ${BIN} | awk '{print $3}')
		BIN_SR_SZ=$(echo ${BIN} | awk '{print $4}')
		BIN_LR_SZ=$(echo ${BIN} | awk '{print $5}')

		RatatoskBin ${BIN_NAME} ${BIN_SR_FILEN} ${BIN_LR_FILEN} ${BIN_SR_SZ} ${BIN_LR_SZ} \
					${NAME_SR_UNMAPPED_IN_FILE} ${PREFIX_PATH_SEG} \
					${NB_THREADS} ${SHORT_READS_PHASING} ${LONG_READS_PHASING}

		if [ $? -eq 0 ]
		then
			if [ "${BIN_NAME}" != "unmapped_sr" ] && [ "${BIN_NAME}" != "unknown_lr" ]
			then
				rm -rf ${BIN_SR_FILEN} ${BIN_LR_FILEN} # Bin correction finished, can delete tmp files
			fi
		else

			1>&2 echo "Bin correction failed. Abort."
			exit 13 # Something happened, exit
		fi

	done < "${PREFIX_PATH_SEG}/sample.bins"

	if [ -f "${PREFIX_PATH_SEG}/sample_lr_unknown.fq" ]
	then

		echo "- Correcting bin of low MAPQ and unmapped long reads"
		# Concat corrected bins
		cat ${PREFIX_PATH_SEG}/sample_lr_*_corrected.fastq > "${PREFIX_PATH_SEG}/sample_lr_map.fastq"; rm -rf ${PREFIX_PATH_SEG}/sample_lr_*_corrected.fastq;
		# Extract short reads
		samtools bam2fq -@ ${NB_THREADS} -n ${SHORT_READS_BAM} | gzip > "${PREFIX_PATH_SEG}/sample_sr.fastq.gz"
		# Final Ratatosk correction
		Ratatosk -c ${NB_THREADS} -s "${PREFIX_PATH_SEG}/sample_sr.fastq.gz" -l "${PREFIX_PATH_SEG}/sample_lr_unknown.fq" -a "${PREFIX_PATH_SEG}/sample_lr_map.fastq" -o "${PREFIX_PATH_SEG}/sample_lr_unknown_corrected"

		if [ $? -eq 0 ] && [ -f "${PREFIX_PATH_SEG}/sample_lr_unknown_corrected.fastq" ]
		then

			echo "- Concatenating corrected bins and cleaning up"
			# Concat corrected bins
			cat "${PREFIX_PATH_SEG}/sample_lr_map.fastq" "${PREFIX_PATH_SEG}/sample_lr_unknown_corrected.fastq" > "${PREFIX_PATH_CORRECTED}/sample_corrected.fastq"
			# Delete tmp files
			rm -rf ${PREFIX_PATH_SEG}
		else

			1>&2 echo "Failed to correct bin of unmapped and low mapq long reads. Abort."
			exit 14
		fi

	else
		# Concat corrected bins
		cat ${PREFIX_PATH_SEG}/sample_lr_*_corrected.fastq > "${PREFIX_PATH_CORRECTED}/sample_corrected.fastq";
		# Delete tmp files
		rm -rf ${PREFIX_PATH_SEG}
	fi
else

	1>&2 echo "Failed to segment input BAM into bins. Abort."
	exit 15
fi