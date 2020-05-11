# Reference-guided correction

While Ratatosk is a reference-free method, a reference-guided preprocessing of the input data is proposed to:
- improve the correction
- reduce the running time
- distribute the workload over many nodes of an HPC

## Overview

1. Short and long reads are mapped to a reference genome.
2. Short and long reads with good mapping quality are distributed over bins corresponding to segments of the reference genome they map to.
3. Bins are corrected independently.
4. Unmapped long reads and long reads with low mapping quality are corrected *de novo* using all short reads. The correction is assisted by the previously corrected long read bins.

## Requirements

* [Slurm](https://slurm.schedmd.com)
* [python3](https://python.org)
* [samtools](https://htslib.org)
* [Ratatosk & Bifrost](https://github.com/GuillaumeHolley/Ratatosk)

The following python3 module is required:
* [pysam](https://pysam.readthedocs.io)

## Usage

1. As the input of the pipeline is BAM files, we recommend [bwa-mem](https://github.com/lh3/bwa) to map the short reads and [minimap2](https://github.com/lh3/minimap2) for the long reads. Input BAM files must be sorted and indexed (`samtools sort` and `samtools index`).

2. Script `scripts/reference_guiding/Ratatosk.sh` is a bash script submitting sbatch commands. It performs the read binning, correction and merging:
```
bash scripts/reference_guiding/Ratatosk.sh <reference_genome.fa> <short_reads.bam> <long_reads.bam> <output_prefix> <slurm_partition>
```

**Input**

- *reference_genome.fa*: FASTA file of the reference genome to which the input short and long reads are mapped to.
- *short_reads.bam*: BAM file of the input paired-end short reads
- *long_reads.bam*: BAM file of the input long reads
- *output_prefix*: Output path where temporary and output files are written. This path must contain enough space to run the binning. See notes below.
- *slurm_partition*: Slurm partition where the sbatch jobs will be submitted. If no partition is given, the default Slurm partition *nomosix* is used.

**Output**

Upon completion of the script, the corrected long reads are available in `<output_prefix>/ratatosk/sample_corrected.fastq`. Folder `<output_prefix>/segments` contains temporary files and can be deleted.

**Important**

Script `scripts/reference_guiding/Ratatosk.sh` is a work in progress. Here are a few considerations to run it that will be subsequently improved:

- This script was designed specifically to run on a human genome dataset. As such, the script assumes that a machine wih at least 48 cores and 350 GB of RAM is available on the Slum partition selected.
- Given a short reads FASTA/FASTQ file of size *X* GB and a long reads FASTA/FASTQ file of size *Y* GB, the output path *output_prefix* must have at least 2.5X + 2.5Y GB of disk available