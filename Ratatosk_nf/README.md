# Ratatosk.nf #

Ratatosk.nf is a Nextflow pipeline designed to run Ratatosk on a multi-node compute architecture (cluster or cloud system).

## Requirements ##

* [Nextflow](https://www.nextflow.io/)
* [Ratatosk](https://github.com/DecodeGenetics/Ratatosk)
* [Samtools](http://www.htslib.org/download/)
* [Pigz](https://zlib.net/pigz/)

The pipeline was extensively tested with Nextflow 22.01, Ratatosk 0.9, samtools 1.17 and pigz 2.7.

## Usage ##

**Before starting**
- Ratatosk works best with paired-end short reads in input (`--sr_fq_in`): **reads from the same pair must have the same FASTA/FASTQ name**.
- Ratatosk was designed primarily for correcting ONT R9.4 reads for which the maximum base quality is 40 (default value). When correcting a different type of long reads, adjust the maximum base quality accordingly with `--max_lr_bq`, e.g `--max_lr_bq 90` must be used with ONT R10.
- Several temporary files are written to disk by the pipeline so make sure to have enough free space.

### Running the pipeline ###

IMPORTANT: See [Cluster configuration](#cluster-configuration) before running `Ratatosk.nf` to configure the pipeline to your cluster system.

```bash
nextflow run -profile cluster Ratatosk.nf \
--in_lr_fq long in_long_reads.fastq.gz --in_sr_fq in_short_reads.fastq.gz --out_dir /my/output/directory/
```
This command runs the `Ratatosk.nf` pipeline with Nextflow: Ratatosk corrects the long read file (`--in_lr_fq long in_long_reads.fastq.gz`) using an index built from the short read file (`--in_sr_fq in_short_reads.fastq.gz`). The output corrected reads are written to `/my/output/directory/lr.corrected.fastq.gz` (`--out_dir /my/output/directory/`).

### Pipeline arguments
The following are mandatory:
- **--in_lr_fq** or **--in_lr_bam**: Long reads to correct, either in FASTQ or BAM format
- **--in_sr_fq** or **--in_sr_bam**: Short reads to use for correction, either in FASTQ or BAM format. If in FASTQ format, **reads from the same pair must have the same FASTA/FASTQ name**
- **--out_dir**: Output directory of the corrected long reads

The following are optional:
- **--max-lr-bq**: Maximum base quality of the input long reads to correct. Default is 40.

Alternatively, one can skip the arguments on the command line and instead, edit the parameter file `params.yaml` which uses the same argument names as the command line. Once the file is edited, the pipeline can be run with the following command:
```bash
nextflow run -profile cluster -params-file params.yaml Ratatosk.nf
```

### Cluster configuration

By default, `Ratatosk.nf` will run jobs on SLURM. You can instead use the workload manager of your choice by:
- adding `-process.executor=your_workload_manager` on the command line
- modifying the value of `cluster.process.executor` in `nextflow.config`

Nextflow supports a wide variety of workload managers and cloud systems: SLURM, SGE, LSF, AWS, Google Cloud, etc. See the [Nextflow executor documentation](https://www.nextflow.io/docs/latest/executor.html) for more information.

The pipeline uses 3 profiles of nodes:
- **small_node**: 16 cores, 4GB of RAM per core. Mostly used for decompressing and splitting FASTQ files, extracting reads from BAM files.
- **medium_node**: 32 cores, 4GB of RAM per core. Used for loading a Ratatosk index and correcting a long read FASTQ chunk.
- **large_node**: 64 cores, 6GB of RAM per core. Used for creating Ratatosk indexes.

These profiles can be edited in `nextflow.config` to fit your cluster configuration. Keep in mind:
- There is exactly one job with the **large_node** profile running at any given time. This job (creating a Ratatosk index) is very CPU, RAM and IO intensive so it is important to give **large_node** your "best" node specs.
- There are up to 50 jobs with the **medium_node** profile running in parallel at any given time. These jobs do the bulk of the correction. The 50 jobs limit can be changed by modifying the value of `pipeline.nb_split_correction_jobs` in `config/pipeline.config`

## Advanced settings

By default, the pipeline assumes that all required tool binaries (ratatosk, samtools and pigz) are available in your PATH. If not:
- the path to each tool binary can be set in `config/tools.config`
```json
tools.ratatosk.bin = '/local/path/to/ratatosk/binary'
tools.samtools.bin = '/local/path/to/samtools/binary'
tools.pigz.bin = '/local/path/to/pigz/binary'
```
- the path to each tool binary can be set in global variables in your `~/.bashrc`
```bash
export RATATOSK=/local/path/to/ratatosk/binary
export SAMTOOLS=/local/path/to/samtools/binary
export PIGZ_BIN=/local/path/to/pigz/binary
```






