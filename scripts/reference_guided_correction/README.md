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
* [Ratatosk](https://github.com/DecodeGenetics/Ratatosk/)

To run on multiple machines:

* [Slurm](https://slurm.schedmd.com)

## Usage

As the input of the pipeline is BAM files, we recommend [bwa-mem](https://github.com/lh3/bwa) to map the short reads and [winnowmap2](https://github.com/marbl/Winnowmap) for the long reads. Input BAM files must be sorted and indexed (`samtools sort` and `samtools index`). Phasing files can also be provided to assist with the correction (see [Phasing](phasing.md)).

### Single node

Script `Ratatosk_Ref_SingleNode.sh` is a bash script performing the read binning, correction and merging in parallel on a single compute node:
```
Usage:

bash Ratatosk_Ref_SingleNode.sh -r <REFRENCE_GENOME> -s <SHORT_READS_BAM> -l <LONG_READS_BAM> -o <OUTPUT_PATH> -c [NB_THREADS] -p [SHORT_READS_PHASING] -P [LONG_READS_PHASING]

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

The script `Ratatosk_Ref_SingleNode.sh` outputs temporary files in `<OUTPUT_PATH>`. Make sure `<OUTPUT_PATH>` has plenty of free space.

### Multiple nodes

Script `Ratatosk_Ref_SLURM.sh` is a bash script performing the read binning, correction and merging in parallel over multiple compute nodes. The script submits sbatch commands to the SLURM workload manager which must be available on your cluster:
```
Usage:

bash Ratatosk_Ref_SLURM.sh -c <MIN_NB_THREADS> -C <MAX_NB_THREADS> -n <MAX_NB_NODES> -r <REFRENCE_GENOME> -s <SHORT_READS_BAM> -l <LONG_READS_BAM> -o <OUTPUT_PATH> -p [SHORT_READS_PHASING] -P [LONG_READS_PHASING] -m [SLURM_PARTITION]

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

The script submits a batch of "small" SLURM jobs requiring 16GB of RAM and 24h wall-clock time. Once completed, a final "large" SLURM job requiring 250GB of RAM and 7 days wall-clock time is submitted. The memory and time requirements must be adapted in the script to fit the genome size and coverage of your input data set. The current requirements are for a 60x short reads and 47x long reads human genome data set.

Corrected long reads are available in `<OUTPUT_PATH>/ratatosk/sample_corrected.fastq` when all jobs are finished. You can check the status of the jobs with:
```
squeue -u <username>
```
After correction, folder `<OUTPUT_PATH>/segments` contains temporary files that can be deleted.

**Important**

The script `Ratatosk_Ref_SLURM.sh` outputs temporary files in `<OUTPUT_PATH>`. Make sure `<OUTPUT_PATH>` has plenty of free space.
