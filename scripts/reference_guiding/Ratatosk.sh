#!/bin/bash

REF_GENOME=$1
SHORT_READS_BAM=$2
LONG_READS_BAM=$3
OUT_PREFIX=$4
PARTITION=$5

BIN_SZ="5000000" # Length of segments in bp
MIN_MAPQ_SR="0" # Minimum MAPQ of short reads (0 to have secondary alignments)
MIN_MAPQ_LR="30" # Minimum MAQ of long reads
MIN_QS_SR="28" # Minimum quality score for a base in the short reads. Replaced by a 'N' in bin if below that threshold

PREP_FILES=0

# Check if short read BAM file exists
if [ ! -f "${SHORT_READS_BAM}" ]
then
	1>&2 echo "Short read BAM file not found. Abort."
	PREP_FILES=1
fi

# Check if long read BAM file exists
if [ ! -f "${LONG_READS_BAM}" ]
then
	1>&2 echo "Long read BAM file not found. Abort."
	PREP_FILES=2
fi

# Check if reference genome file exists
if [ ! -f "${REF_GENOME}" ]
then
	1>&2 echo "Reference genome file not found. Abort."
	PREP_FILES=3
fi

# If all input file exists
if [ "$PREP_FILES" -eq 0 ]
then

	# If partition is empty, set default slurm partition
	if [ -z "${PARTITION}" ]
	then
		PARTITION="nomosix"
	fi

	echo "Input reference genome file: ${REF_GENOME}"
	echo "Input short read file: ${SHORT_READS_BAM}"
	echo "Input long read file: ${LONG_READS_BAM}"
	echo "Output prefix: ${OUT_PREFIX}"
	echo "Slurm partition: ${PARTITION}"

	mkdir "${OUT_PREFIX}"

	PREFIX_PATH_SEG="${OUT_PREFIX}/segments"
	PREFIX_PATH_CORRECTED="${OUT_PREFIX}/ratatosk"
	PREFIX_FILE_BIN="${OUT_PREFIX}/bin.sh"
	NAME_SR_UNMAPPED_IN_FILE="${PREFIX_PATH_SEG}/sample_sr_unmapped"

	# Load chromosomes names and lengths in array
	mapfile -t CHR_NAMES < <(samtools view -H "${SHORT_READS_BAM}" | grep "@SQ" | awk '{print substr($2,4,length($2))}')
	mapfile -t CHR_LENGTHS < <(samtools view -H "${SHORT_READS_BAM}" | grep "@SQ" | awk '{print substr($3,4,length($3))}')

	# Extract bins
	CMD="python3 segmentBAM.py -s ${SHORT_READS_BAM} -l ${LONG_READS_BAM} -t 24 -o ${OUT_PREFIX}/segments -b ${BIN_SZ} -m ${MIN_MAPQ_SR} -n ${MIN_MAPQ_LR} -q ${MIN_QS_SR};"
	CMD="${CMD} samtools bam2fq -@ 24 -n ${SHORT_READS_BAM} | gzip >${OUT_PREFIX}/segments/sample_sr.fastq.gz;"
	
	SLURM_OUT=$(sbatch -p ${PARTITION} -J extractBins --mem=24G --cpus-per-task=24 -t 2-0:0 -o extractBins.slurm.log --wrap="${CMD}")

	# Extract job ID from SLURM output.
	if ! echo ${SLURM_OUT} | grep -q "[1-9][0-9]*$"; then
		echo "Job(s) submission failed."
		echo ${SLURM_OUT}
		exit 1
	else
		JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
	fi

	# Extract graph of unmapped reads
	CMD="Bifrost build -v -t 24 -k 63 -s ${NAME_SR_UNMAPPED_IN_FILE}.fa -o ${NAME_SR_UNMAPPED_IN_FILE}_graph"

	SLURM_OUT=$(sbatch -p ${PARTITION} -J buildUnmapped --dependency=afterok:${JOB_ID} --mem=48G --cpus-per-task=24 -t 1-0:0 -o _buildUnmapped.slurm.log --wrap="${CMD}")

	# Extract job ID from SLURM output.
	if ! echo ${SLURM_OUT} | grep -q "[1-9][0-9]*$"; then
		echo "Job(s) submission failed."
		echo ${SLURM_OUT}
		exit 1
	else
		JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
	fi

	# Bin correction
	echo "#!/bin/bash" > ${PREFIX_FILE_BIN}
	echo "#SBATCH --partition=${PARTITION}" >> ${PREFIX_FILE_BIN}
	echo "#SBATCH --job-name=correctSegONT1" >> ${PREFIX_FILE_BIN}
	echo "#SBATCH --mem=16G" >> ${PREFIX_FILE_BIN}
	echo "#SBATCH --time=2-0:0" >> ${PREFIX_FILE_BIN}
	echo "#SBATCH --cpus-per-task=8" >> ${PREFIX_FILE_BIN}

	# For all ONT reads with mapq good enough, correct chromsome by chromosome, then bin by bin
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
	
				CMD="if [ -f ${NAME_LR_IN_FILE}.fq ] && [ -s ${NAME_LR_IN_FILE}.fq ]; then if [ -f ${NAME_SR_IN_FILE}.fa ] && [ -s ${NAME_SR_IN_FILE}.fa ]; then" # Check input file exists as they should
				CMD="${CMD} /usr/bin/time -v /odinn/users/guillaumeh/software/Ratatosk/build/src/Ratatosk -v -c 8 -q 13 -i ${NAME_SR_IN_FILE}.fa -u ${NAME_SR_UNMAPPED_IN_FILE}.fa -m ${NAME_SR_UNMAPPED_IN_FILE}_graph.gfa -l ${NAME_LR_IN_FILE}.fq -o ${NAME_LR_OUT_FILE};"  # Ratatosk correction
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
		exit 1
	else
		JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
	fi

	CMD="cat \$(ls -t ${PREFIX_PATH_SEG}/sample_lr_*_corrected.fastq) > ${PREFIX_PATH_SEG}/sample_lr_map.fastq;" # Concat corrected bins
	CMD="${CMD} Ratatosk -v -c 48 -q 13 -i ${OUT_PREFIX}/segments/sample_sr.fastq.gz -l ${PREFIX_PATH_SEG}/sample_lr_unknown.fq -p ${PREFIX_PATH_SEG}/sample_lr_map.fastq -o ${PREFIX_PATH_SEG}/sample_lr_unknown_corrected;" # Ratatosk
	CMD="${CMD} cat ${PREFIX_PATH_SEG}/sample_lr_map.fastq ${PREFIX_PATH_SEG}/sample_lr_unknown_corrected.fastq > ${PREFIX_PATH_SEG}/sample_corrected.fastq" # Concat corrected bins

	SLURM_OUT=$(sbatch -p ${PARTITION} -J correctSegONT2 --dependency=afterok:${JOB_ID} --mem=350G --cpus-per-task=48 -t 7-0:0 -o correctSegONT2.slurm.log --wrap="${CMD}")

	# Extract job ID from SLURM output.
	if ! echo ${SLURM_OUT} | grep -q "[1-9][0-9]*$"; then
		echo "Job(s) submission failed."
		echo ${SLURM_OUT}
		exit 1
	else
		JOB_ID=$(echo ${SLURM_OUT} | grep -oh "[1-9][0-9]*$")
	fi
fi
