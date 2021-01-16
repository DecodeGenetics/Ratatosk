# Reference-guided correction

An optional reference-guided preprocessing of the input data is proposed to:
- improve the correction
- reduce the running time
- distribute the workload over nodes of an HPC

**Important**: A reference-guided correction is **not** recommended if your reference genome is fragmented

## Overview

1. Short and long reads are mapped to a reference genome.
2. Short and long reads with good mapping quality are distributed over bins corresponding to segments of the reference genome they map to.
3. Bins are corrected independently.
4. Unmapped and low mapping quality long reads are corrected *de novo* using all short reads. The correction is assisted by the previously corrected long read bins.

## Requirements

* [python3](https://python.org) with the [pysam](https://pysam.readthedocs.io) module
* [samtools](https://htslib.org)
* [Ratatosk](https://github.com/GuillaumeHolley/Ratatosk)

To run on multiple machines:

* [Slurm](https://slurm.schedmd.com)

## Usage

As the input of the pipeline is BAM files, we recommend [bwa-mem](https://github.com/lh3/bwa) to map the short reads and [winnowmap2](https://github.com/marbl/Winnowmap) for the long reads. Input BAM files must be sorted and indexed (`samtools sort` and `samtools index`). Phasing files can also be provided to assist with the correction (see [Phasing](https://github.com/DecodeGenetics/Ratatosk/tree/master/phasing.md)).

### Single machine

Script `Ratatosk_SingleNode.sh` is a bash script performing the read binning, correction and merging in parallel:
```
Usage:

bash Ratatosk_SingleNode.sh -r <REFRENCE_GENOME> -s <SHORT_READS_BAM> -l <LONG_READS_BAM> -o <OUTPUT_PATH> -c [NB_THREADS] -p [SHORT_READS_PHASING] -P [LONG_READS_PHASING]

<Mandatory>: 
-r <REFRENCE_GENOME>: Reference genome in FASTA
-s <SHORT_READS_BAM>: Input BAM file of the short reads. Short reads from the same pair must have the same name.
-l <LONG_READS_BAM>: Input BAM file of the long reads.
-o <OUTPUT_PATH>: Output path.

[Optional]: 
-c [NB_THREADS]: Number of threads to use (default: 1)
-p [SHORT_READS_PHASING]: Short reads phasing (default: None)
-P [LONG_READS_PHASING]: Long reads phasing (default: None)
```

Corrected long reads are output to `<OUTPUT_PATH>/ratatosk/sample_corrected.fastq`.

**Important**

Script `Ratatosk_SingleNode.sh` is a work in progress. Here are a few considerations to run it that will be subsequently improved:

- Given a short reads FASTA/FASTQ file of size *X* GB and a long reads FASTA/FASTQ file of size *Y* GB, the output path *output_prefix* must have at least **2.5X + 2.5Y GB** of disk available.

### Multiple machines

Script `Ratatosk_SLURM.sh` is a bash script submitting sbatch commands. It performs the read binning, correction and merging:
```
Usage:

bash Ratatosk_SLURM.sh -c <MIN_NB_THREADS> -C <MAX_NB_THREADS> -n <MAX_NB_NODES> -r <REFRENCE_GENOME> -s <SHORT_READS_BAM> -l <LONG_READS_BAM> -o <OUTPUT_PATH> -p [SHORT_READS_PHASING] -P [LONG_READS_PHASING] -m [SLURM_PARTITION]

<Mandatory>: 
-c <MIN_NB_THREADS>: Minimum number of threads that can be given to a job.
-C <MAX_NB_THREADS>: Maximum number of threads that can be given to a job.
-n <MAX_NB_NODES>: Maximum number of nodes that can be used at a given time.
-r <REFRENCE_GENOME>: Reference genome in FASTA
-s <SHORT_READS_BAM>: Input BAM file of the short reads. Short reads from the same pair must have the same name.
-l <LONG_READS_BAM>: Input BAM file of the long reads.
-o <OUTPUT_PATH>: Output path.

[Optional]: 
-m [SLURM_PARTITION]: Slurm partition (default: nomosix).
-p [SHORT_READS_PHASING]: Short reads phasing (default: None)
-P [LONG_READS_PHASING]: Long reads phasing (default: None)
```

Corrected long reads are available in `<output_prefix>/ratatosk/sample_corrected.fastq` when all jobs are finished. You can check the status of the jobs with:
```
squeue -u <username>
```
After correction, folder `<output_prefix>/segments` contains temporary files that can be deleted.


**Important**

Script `Ratatosk_SLURM.sh` is a work in progress. Here are a few considerations to run it that will be subsequently improved:

- This script was designed specifically to run on a human genome dataset. As such, the script assumes that a machine with at least **48 cores** and **350 GB of RAM** is available on the Slum partition selected.
- Given a short reads FASTA/FASTQ file of size *X* GB and a long reads FASTA/FASTQ file of size *Y* GB, the output path *output_prefix* must have at least **2.5X + 2.5Y GB** of disk available.