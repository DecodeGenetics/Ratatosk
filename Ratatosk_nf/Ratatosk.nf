#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process chunkLR {

	label 'small_node'

	input:

		path bam_in
		path fastq_in

	output:

		path "lr.part_?????????"

	shell '/bin/bash', '-euo', 'pipefail'

	script:

		if (bam_in) {
			"""
			samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}

			NB_LINES=\$(\${samtools} bam2fq -@ ${task.cpus} -n ${bam_in} | wc -l | awk '{print \$1}')
			NB_READS_PER_JOB=\$((NB_LINES/${params.pipeline.nb_split_correction_jobs}))
			NB_FASTQ_LINES_PER_JOB=\$(((NB_READS_PER_JOB/4)*4))

			\${samtools} bam2fq -@ ${task.cpus} -n ${bam_in} | split -l \${NB_FASTQ_LINES_PER_JOB} -a 9 -d - lr.part_
			"""
		}
		else if (fastq_in) {
			"""
			pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

			SYMLINKED_FN=\$(readlink -f \"${fastq_in}\")
			IS_GZ_COMPRESSED=\$(file \${SYMLINKED_FN} | { grep \"gzip compressed\" || true; } | wc -l | awk '{print \$1}')
			
			if [ \${IS_GZ_COMPRESSED} -eq 1 ]
			then
			
				NB_LINES=\$(\${pigz} -p ${task.cpus} -c -d ${fastq_in} | wc -l | awk '{print \$1}')
				NB_READS_PER_JOB=\$((NB_LINES/${params.pipeline.nb_split_correction_jobs}))
				NB_FASTQ_LINES_PER_JOB=\$(((NB_READS_PER_JOB/4)*4))
			
				\${pigz} -p ${task.cpus} -c -d ${fastq_in} | split -l \${NB_FASTQ_LINES_PER_JOB} -a 9 -d - lr.part_
			else
			
				NB_LINES=\$(wc -l ${fastq_in} | awk '{print \$1}')
				NB_READS_PER_JOB=\$((NB_LINES/${params.pipeline.nb_split_correction_jobs}))
				NB_FASTQ_LINES_PER_JOB=\$(((NB_READS_PER_JOB/4)*4))
			
				split -l \${NB_FASTQ_LINES_PER_JOB} -a 9 -d ${fastq_in} lr.part_
			fi
			"""
		}
		else error("No reads in BAM or FASTQ format provided in input")
}

process getSuffixSampleSR {

	executor = 'local' 

    memory = '10.MB'
    cpus = 1
    time = '1.m'

    input:
    	path sr_fq

    output:
    	stdout

    //shell '/bin/bash', '-euo', 'pipefail'

    """
	SUFFIX_1ST_SR_REC=\$(zcat -f ${sr_fq} | head -n 1 | awk '{print substr(\$0,length(\$0)-1,2)}')

	echo \"\${SUFFIX_1ST_SR_REC}\"
    """
}

process extractSR {

	label 'small_node'

	input:

		path sr_bam_in

	output:

		path "sr.fq.gz"

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	samtools=\${SAMTOOLS:-${params.tools.samtools.bin}}
	pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

	\${samtools} bam2fq -@ ${task.cpus} -n ${sr_bam_in} | \${pigz} -p ${task.cpus} -c > sr.fq.gz
	"""
}

process buildIndex_1 {

	label 'large_node'

	input:

		path "lr_chunks_*"
		path sr_fq

	output:

		path('ratatosk.index.k31.rtsk'), emit: ratatosk_index1

		tuple path('ratatosk.index.k31.fasta.gz'), path('ratatosk.index.k31.bfi'), emit: bifrost_index1
		tuple path('ratatosk.index.k63.fasta.gz'), path('ratatosk.index.k63.bfi'), emit: bifrost_index2


	shell '/bin/bash', '-euo', 'pipefail'

	"""
	ratatosk=\${RATATOSK:-${params.tools.ratatosk.bin}}

	find . -name 'lr_chunks_*' > lr_chunks.list
	\${ratatosk} index -c ${task.cpus} -v -1 -Q ${params.max_lr_bq} -s ${sr_fq} -l lr_chunks.list -o ratatosk

	if [ ! -s \"ratatosk.index.k31.fasta.gz\" ] || [ ! -s \"ratatosk.index.k31.rtsk\" ] || [ ! -s \"ratatosk.index.k63.fasta.gz\" ]
	then
		echo \"Something went wrong during construction of the index 1: (some) output index files are either missing or empty\" 1>&2;
		exit 1
	fi
	"""
}

process correctChunk_1 {

	label 'medium_node'

	input:

		tuple val(id), path(lr_chunk_fq), path('index1.fasta.gz'), path('index1.bfi'), path('index1.rtsk')

	output:

		path "${id}.2.fastq"

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	ratatosk=\${RATATOSK:-${params.tools.ratatosk.bin}}

	\${ratatosk} correct -c ${task.cpus} -v -1 -Q ${params.max_lr_bq} -g index1.fasta.gz -d index1.rtsk -l ${lr_chunk_fq} -o ${id}

	if [ ! -s \"${id}.2.fastq\" ]
	then
		echo \"Something went wrong during correction with index 1\" 1>&2;
		exit 1
	fi
	"""	
}

process buildIndex_2 {

	label 'large_node'

	input:
		path "lr_chunks_*"
		tuple path('index2.fasta.gz'), path('index2.bfi')

	output:
		path('ratatosk.index.k63.rtsk'), emit: ratatosk_index2


	shell '/bin/bash', '-euo', 'pipefail'

	"""
	ratatosk=\${RATATOSK:-${params.tools.ratatosk.bin}}

	find . -name 'lr_chunks_*' > lr_chunks.list
	\${ratatosk} index -c ${task.cpus} -v -2 -Q ${params.max_lr_bq} -g index2.fasta.gz -l lr_chunks.list -o ratatosk

	if [ ! -s \"ratatosk.index.k63.rtsk\" ]
	then
		echo \"Something went wrong during construction of the index 2: (some) output index files are either missing or empty\" 1>&2;
		exit 1
	fi
	"""
}

process correctChunk_2 {

	label 'medium_node'

	input:

		tuple val(id), path(lr_chunk1_fq), path(lr_chunk_fq), path('index2.fasta.gz'), path('index2.bfi'), path('index2.rtsk') // lr_chunk_fq: Input raw chunk, lr_chunk1_fq: Corrected chunk with index 1

	output:

		path "${id}.fastq.gz"

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	ratatosk=\${RATATOSK:-${params.tools.ratatosk.bin}}
	pigz=\${PIGZ_BIN:-${params.tools.pigz.bin}}

	\${ratatosk} correct -c ${task.cpus} -v -2 -G -Q ${params.max_lr_bq} -g index2.fasta.gz -d index2.rtsk -l ${lr_chunk1_fq} -L ${lr_chunk_fq} -o ${id}

	if [ -s \"${id}.fastq.gz\" ]
	then

		\${pigz} -p ${task.cpus} -t ${id}.fastq.gz

		if [ ! \$? -eq 0 ]
		then
			echo \"Output file is not properly compressed\" 1>&2;
			exit 1
		fi
	else

		echo \"Something went wrong during correction with index 2\" 1>&2;
		exit 1
	fi
	"""	
}

process mergeCorrectedChunks {

	publishDir "${params.out_dir}"
	label 'small_node'

	input:
		path "lr_corrected_chunks_*"

	output:
		path "lr.corrected.fastq.gz"

	shell '/bin/bash', '-euo', 'pipefail'

	"""
	cat lr_corrected_chunks_* > lr.corrected.fastq.gz
	"""
}

workflow {

	if (!params.in_lr_bam && !params.in_lr_fq) error("No long reads in input.")
	if (params.in_lr_bam && params.in_lr_fq) error("Input long reads can be in BAM or FASTQ format but not both.")

	if (!params.in_sr_bam && !params.in_sr_fq) error("No short reads in input.")
	if (params.in_sr_bam && params.in_sr_fq) error("Input short reads can be in BAM or FASTQ format but not both.")

	lr_bam = params.in_lr_bam ? [params.in_lr_bam] : []
	lr_fq = params.in_lr_fq ? [params.in_lr_fq] : []

	// Split input long read fastq into smaller chunks
	lr_chunked = chunkLR(lr_bam, lr_fq).collect()

	// Extract short reads
	sr_fq = params.in_sr_fq ? Channel.fromPath(params.in_sr_fq) : extractSR(Channel.fromPath(params.in_sr_bam))

	getSuffixSampleSR(sr_fq)
		.subscribe {
			if ("$it"=="/1\n" || "$it"=="/2\n") println "\nWARNING: /1 or /2 suffix detected in the input short read names. " +
														"Ratatosk REQUIRES that input short reads from the same pair have EXACTLY the same FASTQ name. " +
														"Make sure the input FASTQ names are correctly formatted.\n"
		}

	// Build Ratatosk index for 1st correction pass
	index1 = buildIndex_1(lr_chunked, sr_fq)

	// Associate each long read chunk to an identifier
	lr_chunked.flatten().map{ file -> tuple(file.name.substring(file.name.length()-9), file) }.groupTuple().set { lr_id2chunk }
	// Associate each (id, chunk) to the Ratatosk index 1
	lr_id2chunk.combine(index1.bifrost_index1).combine(index1.ratatosk_index1).set { lr_chunked_split }

	// Correct each chunk (1st pass)
	lr_chunked_corr1 = correctChunk_1(lr_chunked_split).collect()

	// Build Ratatosk index for 2nd correction pass
	index2 = buildIndex_2(lr_chunked_corr1, index1.bifrost_index2)

	// Associate each corrected long read chunk to an identifier
	lr_chunked_corr1.flatten().map{ file -> tuple(file.simpleName, file) }.groupTuple().set { lr_corr1_id2chunk }
	// Associate each (id, corrected chunk) to the uncorrected chunk it is from
	lr_corr1_id2chunk.join(lr_id2chunk).set { lr_corr1_id2chunk2 }
	// Associate each (id, corrected chunk, uncorrected chunk) to the Ratatosk index 2
	lr_corr1_id2chunk2.combine(index1.bifrost_index2).combine(index2.ratatosk_index2).set { lr_chunked_corr1_split }

	// Correct each chunk (2nd pass)
	lr_chunked_corr2 = correctChunk_2(lr_chunked_corr1_split).collect()
	// Merge corrected chunks
	lr_chunked_corr = mergeCorrectedChunks(lr_chunked_corr2)
}