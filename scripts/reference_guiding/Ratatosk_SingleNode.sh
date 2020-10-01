#!/bin/bash

NB_THREADS=$1
REF_GENOME=$2
SHORT_READS_BAM=$3
LONG_READS_BAM=$4
OUT_PREFIX=$(mkdir -p $5; cd $5; pwd)

MAX_NB_CORES=$(nproc --all)

PREP_FILES=0

# Check if short read BAM file exists
if [ ! -f "${SHORT_READS_BAM}" ]
then
	1>&2 echo "Short read BAM file not found. Abort."
	exit 1
fi

# Check if long read BAM file exists
if [ ! -f "${LONG_READS_BAM}" ]
then
	1>&2 echo "Long read BAM file not found. Abort."
	exit 2
fi

# Check if reference genome file exists
if [ ! -f "${REF_GENOME}" ]
then
	1>&2 echo "Reference genome file not found. Abort."
	exit 3
fi

if [ "${NB_THREADS}" -gt "${MAX_NB_CORES}" ]
then
	1>&2 echo "Number of threads required exceed what is available. Abort."
	exit 4
fi

RatatoskBin() {

	if [ "${1}" != "unmapped_sr" ] && [ "${1}" != "unknown_lr" ] # Bin is not unmapped short reads or unmapped/low mapq long reads
	then
		if [ ! "${4}" -eq 0 ] && [ ! "${5}" -eq 0 ] # The bin contains at least one long and one short read to perform the correction
		then
			if [ -f "${6}" ] # Unmapped short reads available to assist correction
			then
				Ratatosk -c "${8}" -q 13 -s "${2}" -l "${3}" -u "${6}" -o "${7}/sample_lr_${1}_corrected"
			else
				Ratatosk -c "${8}" -q 13 -s "${2}" -l "${3}" -o "${7}/sample_lr_${1}_corrected"
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

echo "Number of threads: ${NB_THREADS}"
echo "Input reference genome file: ${REF_GENOME}"
echo "Input short read file: ${SHORT_READS_BAM}"
echo "Input long read file: ${LONG_READS_BAM}"
echo "Output prefix: ${OUT_PREFIX}"

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

		RatatoskBin ${BIN_NAME} ${BIN_SR_FILEN} ${BIN_LR_FILEN} ${BIN_SR_SZ} ${BIN_LR_SZ} ${NAME_SR_UNMAPPED_IN_FILE} ${PREFIX_PATH_SEG} ${NB_THREADS}

		if [ $? -eq 0 ]
		then
			if [ "${BIN_NAME}" != "unmapped_sr" ] && [ "${BIN_NAME}" != "unknown_lr" ]
			then
				rm -rf ${BIN_SR_FILEN} ${BIN_LR_FILEN} # Bin correction finished, can delete tmp files
			fi
		else

			1>&2 echo "Bin correction failed. Abort."
			exit 5 # Something happened, exit with error code 5
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
		Ratatosk -c ${NB_THREADS} -q 13 -s "${PREFIX_PATH_SEG}/sample_sr.fastq.gz" -l "${PREFIX_PATH_SEG}/sample_lr_unknown.fq" -a "${PREFIX_PATH_SEG}/sample_lr_map.fastq" -o "${PREFIX_PATH_SEG}/sample_lr_unknown_corrected"

		if [ $? -eq 0 ] && [ -f "${PREFIX_PATH_SEG}/sample_lr_unknown_corrected.fastq" ]
		then

			echo "- Concatenating corrected bins and cleaning up"
			# Concat corrected bins
			cat "${PREFIX_PATH_SEG}/sample_lr_map.fastq" "${PREFIX_PATH_SEG}/sample_lr_unknown_corrected.fastq" > "${PREFIX_PATH_CORRECTED}/sample_corrected.fastq"
			# Delete tmp files
			rm -rf ${PREFIX_PATH_SEG}
		else

			1>&2 echo "Failed to correct bin of unmapped and low mapq long reads. Abort."
			exit 6
		fi

	else
		# Concat corrected bins
		cat ${PREFIX_PATH_SEG}/sample_lr_*_corrected.fastq > "${PREFIX_PATH_CORRECTED}/sample_corrected.fastq";
		# Delete tmp files
		rm -rf ${PREFIX_PATH_SEG}
	fi
else

	1>&2 echo "Failed to segment input BAM into bins. Abort."
	exit 7
fi