# Multi nodes *de novo* correction

The following describes how to perform *de novo* correction of a long read data set with Ratatosk over multiple compute nodes.

## Overview

1. The long read data set is split into batches (single node)
2. A Ratatosk index is built from all the input short and long reads for the **1st** correction pass (single node)
3. Each batch of long reads is corrected using the previously built index (multiple nodes)
4. A Ratatosk index is built from all the short and **corrected** long reads for the **2nd** correction pass (single node)
5. Each batch of corrected long reads is re-corrected using the previously built index (multiple nodes)
6. Batches of corrected long reads are put back into one file (single node)

Steps 2 and 4 are the most memory consuming but can only be performed on single nodes.

## Requirements

* [Slurm](https://slurm.schedmd.com)

## Usage

Script `Ratatosk_denovo_SLURM.sh` is a bash script performing the read splitting, correction and merging in parallel over multiple compute nodes. The script submits commands to the SLURM workload manager:
```
Usage:

bash Ratatosk_denovo_SLURM.sh -c <MIN_NB_THREADS> -C <MAX_NB_THREADS> -n <MAX_NB_NODES> -s <SHORT_READS_FQ> -l <LONG_READS_FQ> -o <OUTPUT_PATH> -m [SLURM_PARTITION]

<Mandatory>: 
-c <MIN_NB_THREADS>: Minimum number of threads that can be given to a job.
-C <MAX_NB_THREADS>: Maximum number of threads that can be given to a job.
-n <MAX_NB_NODES>: Maximum number of nodes that can be used at a given time.
-s <SHORT_READS_FQ>: Input FASTA/FASTQ file of the short reads. Short reads from the same pair must have the same name.
-l <LONG_READS_FQ>: Input FASTA/FASTQ file of the long reads.
-o <OUTPUT_PATH>: Output path.

[Optional]: 
-m [SLURM_PARTITION]: Slurm partition (default: nomosix).
```

The memory and time requirements must be adapted in the script to fit the genome size and coverage of your input data set. The current requirements are for a 60x short reads human genome data set.

Corrected long reads are available in `<OUTPUT_PATH>/ratatosk/sample_corrected.fastq` when all jobs are finished. You can check the status of the jobs with:
```
squeue -u <username>
```

**Important**

The script `Ratatosk_denovo_SLURM.sh` outputs temporary files in `<OUTPUT_PATH>`. Make sure `<OUTPUT_PATH>` has plenty of free space.
