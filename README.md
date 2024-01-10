# Ratatosk

### Hybrid error correction of long reads using colored de Bruijn graphs

Ratatosk is a *de novo* error correction tool for erroneous long reads designed for accurate variant calling and assembly. It is based on a compacted and colored de Bruijn graph built from accurate short reads. Reads color paths in the graph while vertices are annotated with candidate *de novo* SNPs and short repeats. We demonstrate that Ratatosk can reduce the raw error rate of ONT reads several fold on average with a mean error rate as low as 1.4%. Variant calling on Ratatosk corrected data shows 99.91% and 95.88% F1 for SNP and indels respectively. An assembly of the Ashkenazi individual HG002 created from Ratatosk corrected ONT reads yields a contig N50 of 39.7 Mbp and a quality value of 48.5.

## Table of Contents

* [Installation](#installation)
  * [Using conda](#using-conda)
  * [From source](#from-source)
* [Usage](#usage)
* [Variant calling](#variant-calling)
* [Interface](#interface)
* [FAQ](#faq)
* [Troubleshooting](#troubleshooting)
* [Citation](#citation)
* [Contact](#contact)
* [License](#license)

## Installation

### Using conda
Ratatosk is available on [bioconda](https://bioconda.github.io/recipes/ratatosk/README.html):

```sh
# Can also replace `conda` with `mamba`
conda create -n ratatosk -c conda-forge -c bioconda ratatosk
```
### From source

#### Requirements

* C++11 compiler:
    * [GCC](https://gcc.gnu.org/) >= 5.1.0
    * [Clang](http://clang.llvm.org/) >= 3.5
* [Cmake](https://cmake.org/) >= 2.8.12
* [Zlib](https://zlib.net/)

All are probably already installed on your computer as those are installed by default on most operating systems. They can be downloaded and installed by following the instructions on their respective websites. However, it is most likely they are all available via a package manager for your operating system: 

* **Ubuntu/Debian**:
```
sudo apt-get install build-essential cmake zlib1g-dev
```
* **MacOS** (with [Homebrew](https://brew.sh/)):
```
brew install --with-toolchain llvm
brew install cmake zlib
```
* **Windows 10**:

1. Open the Windows Store
2. Search and install the `Ubuntu` app (from `Canonical Group Limited`)
3. Open the Windows main menu and open the `Ubuntu` app (it should open an Ubuntu terminal)
4. Use the following command in the Ubuntu terminal:
```
sudo apt-get install build-essential cmake zlib1g-dev
```
5. Use the opened Ubuntu terminal for compiling, installing and running Ratatosk (see next section). See [Troubleshooting](#troubleshooting) if you have any problem during the installation.

1. Clone the Git repository
  ```
  git clone --recursive https://github.com/DecodeGenetics/Ratatosk.git
  cd Ratatosk
  ```
2. Install Ratatosk
  ```
  mkdir build && cd build
  cmake ..
  make
  make install
  ```

By default, the installation creates:
* a binary (*Ratatosk*)

**Notes**
`make install` might require `sudo` (`sudo make install`) to proceed. If you want to install Ratatosk in a non-default path, add the option `-DCMAKE_INSTALL_PREFIX=/some/path/ ..` to the `cmake` command where `/some/path/` is where you want to see the Ratatosk files installed. Do not forget to add this path to your environment variables (see [Troubleshooting](#troubleshooting)). If you encounter any problem during the installation, see the [Troubleshooting](#troubleshooting) section.

## Usage

**Before starting**
- Ratatosk works best with paired-end short reads in input (`-s`): **reads from the same pair must have the same FASTA/FASTQ name** (if the reads are extracted from a BAM file, use `samtools bam2fq -n`).
- Ratatosk was designed primarily for correcting ONT R9.4 reads for which the maximum base quality is 40 (default value). When correcting a different type of long reads, adjust the maximum base quality accordingly with `-Q`, e.g `-Q 90` must be used with ONT R10.
- Several temporary files are written in the same repository as the output file (`-o`) so make sure the output folder has plenty of free space.

### Single compute node - one step

```
Ratatosk correct -v -G -c 16 -s short_reads.fastq -l in_long_reads.fastq -o out_long_reads
```
Ratatosk corrects (`Ratatosk correct`) the long read file (`-l in_long_reads.fastq`) with 16 threads (`-c 16`) using an index built from the short read file (`-s short_reads.fastq`). Information messages are printed during the execution (`-v`) and the corrected long reads are written to the compressed file *out_long_reads.fastq.gz* (`-G -o out_long_reads`).

### Single compute node - two steps

The correction can be split in two steps which can be run on different compute nodes in the order given below. It can be beneficial if there is a time limit on the used compute nodes.
```
Ratatosk correct -1 -v -c 16 -s short_reads.fastq -l in_long_reads.fastq -o out_long_reads
Ratatosk correct -2 -v -G -c 16 -s short_reads.fastq -l out_long_reads.2.fastq -L in_long_reads.fastq -o out_long_reads
```
These commands split the correction in the two different correction passes of Ratatosk (`-1` and `-2`). The first command is likely to be the most memory and time consuming of the two.

### Single compute node - four steps

The correction can be split in four steps which can be run on different compute nodes in the order given below. It is sometimes beneficial if there is a time limit on the used compute nodes.
```
Ratatosk index -1 -v -c 16 -s short_reads.fastq -l in_long_reads.fastq -o out_long_reads
Ratatosk correct -1 -v -c 16 -g out_long_reads.index.k31.fasta -d out_long_reads.index.k31.rtsk -l in_long_reads.fastq -o out_long_reads
Ratatosk index -2 -v -c 16 -g out_long_reads.index.k63.fasta -l out_long_reads.2.fastq -o out_long_reads
Ratatosk correct -2 -v -G -c 16 -g out_long_reads.index.k63.fasta -d out_long_reads.index.k63.rtsk -l out_long_reads.2.fastq -L in_long_reads.fastq -o out_long_reads
```
These commands split the correction in the two different correction passes of Ratatosk (`-1` and `-2`) and each correction pass is split into its indexing part (`index`) and correction part (`correct`).

### Multiple compute nodes - cluster/cloud

See our [Nextflow pipeline](Ratatosk_nf) to run the correction in parallel over multiple compute nodes.

### Options

- **Base quality** (`-Q`)

  By default, Ratatosk considers that the maximum base quality value of long reads is 40 but this can be changed to a different value using `-Q`. For example, newer ONT dataset can use a maximum quality vaue of 90 in which case `-Q 90` **must** be used.

- **Insert size** (`-i`)

  By default, Ratatosk uses an estimated paired-end short reads insert size (read length * 2 + fragment size) of about 500bp. You can set the exact insert size with `-i`. If your input short reads are single end reads, set the insert size to the read length.

- **Trimming/Splitting** (`-t`)

  By default, Ratatosk outputs all bases (corrected and uncorrected). By using `-t`, bases with a low correction quality score are trimmed and split. Specifically, given a minimum quality score *Q* (`-t Q`), only subsequences of the corrected long reads for which the bases have a correction quality score >=*Q* are output. Each output subsequence will have `@name/i` as name for which `name` is the input name of the long read and `i` is an integer subsequence ID for read `name`. Note that only subsequences larger than the *k2*-mer size in Ratatosk (63) are output.

### Advanced options

The default *k1*/*k2*-mer lengths (1st/2nd correction passes) are 31/63. To work with larger *k*-mers (using `-k/-K`), you must compile Ratatosk with a larger `MAX_KMER_SIZE` parameter where `MAX_KMER_SIZE=round(k2 + 1, 32)`, i.e. (*k2* + 1) rounded to the larger multiple of 32. Specifying `MAX_KMER_SIZE` at compilation is done as follows when entering the `cmake` command:
```
cmake -DMAX_KMER_SIZE=96 ..
```
In this example, the maximum *k1*/*k2*-mer length allowed is 95.

## Variant calling

See [Variant calling](variant_calling.md) to call SNP and indels from Ratatosk-corrected long reads.

## Interface

```
Ratatosk --help
```
displays the command line interface:
```
Ratatosk x.y.z

Hybrid error correction of long reads using colored de Bruijn graphs

Usage: Ratatosk [COMMAND] [PARAMETERS]
Usage: Ratatosk --help
Usage: Ratatosk --version
Usage: Ratatosk --cite

[COMMAND]:

   correct                         Correct long reads with short reads
   index                           Prepare a Ratatosk index (advanced)

Use "Ratatosk [COMMAND] --help" to get a specific command help
```

Two commands are available: `correct` and `index`. Command `index` is only useful when correcting a data set in multiple steps or over multiple compute nodes (see [multiple machines de novo correction](scripts/multi_nodes_denovo_correction)).

```
Ratatosk correct --help
```
displays the command line interface for the `correct` command:
```
Ratatosk x.y.z

Hybrid error correction of long reads using colored de Bruijn graphs

Usage: Ratatosk [COMMAND] [PARAMETERS]
Usage: Ratatosk --help
Usage: Ratatosk --version
Usage: Ratatosk --cite

[COMMAND]: correct

[PARAMETERS]:

   > Mandatory with required argument:

   -s, --in-short                  Input short read file in fasta/fastq(.gz)
                                   List of input short read files (one file per line)
   -l, --in-long                   Input long read file to correct in fasta/fastq(.gz)
                                   List of input long read files to correct (one file per line)
   -o, --out-long                  Output corrected long read file

   > Optional with required argument:

   -c, --cores                     Number of cores (default: 1)
   -S, --subsampling               Rate of short reads subsampling (default: Auto)
   -t, --trim-split                Trim and split bases with quality score < t (default: no trim/split)
                                   Only sub-read with length >= 63 are output if used
   -u, --in-unmapped-short         Input read file of the unmapped short reads (FASTA/FASTQ possibly gzipped)
                                   List of input read files of the unmapped short reads (one file per line)
   -a, --in-accurate-long          Input high quality long read file (FASTA/FASTQ possibly gzipped)
                                   List of input high quality long read files (one file per line)
                                   (Those reads are NOT corrected but assist the correction of reads in input)
   -g, --in-graph                  Load graph file prepared with the index command
   -d, --in-unitig-data            Load unitig data file prepared with the index command
   -Q, --max-base-qual             Maximum base quality of input long reads (default: 40)

   > Optional with no argument:

   -G, --gzip-out                  Output file is compressed with gzip
   -O, --force-io-order            Force same long read input/output order
   -v, --verbose                   Print information

[ADVANCED PARAMETERS]:

   > Optional with required argument:

   -m, --min-conf-snp-corr         Minimum confidence threshold to correct a SNP (default: 0.9)
   -M, --min-conf-color2           Minimum confidence threshold to color vertices for 2nd pass (default: 0)
   -C, --min-len-color2            Minimum length of a long read to color vertices for 2nd pass (default: 3000)
   -i, --insert-sz                 Insert size of the input paired-end short reads (default: 500)
   -k, --k1                        Length of short k-mers for 1st pass (default: 31)
   -K, --k2                        Length of long k-mers for 2nd pass (default: 63)
   -w, --max-len-weak1             Do not correct non-solid regions >= w bases during 1st pass (default: 1000)
   -W, --max-len-weak2             Do not correct non-solid regions >= w bases during 2nd pass (default: 5000)

   > Optional with no argument:

   -1, --1st-pass-only             Perform *only* the 1st correction pass (default: false)
   -2, --2nd-pass-only             Perform *only* the 2nd correction pass (default: false)
   -F, --no-snp-correction         Disable SNP detection and correction
   -I, --no-graph-index            Disable graph index output

[EXPERIMENTAL PARAMETERS]:

   > Optional with required argument:

   -L, --in-long_raw               Input long read file from 1st pass (FASTA/FASTQ possibly gzipped)
                                   List of input long read files to correct (one file per line)
   -p, --in-short-phase            Input short read phasing file (diploid only)
                                   List of input short read phasing files (one file per line)
   -P, --in-long-phase             Input long read phasing file (diploid only)
                                   List of input long read phasing files (one file per line)
```

```
Ratatosk index --help
```
displays the command line interface for the `index` command:

```
Ratatosk x.y.z

Hybrid error correction of long reads using colored de Bruijn graphs

Usage: Ratatosk [COMMAND] [PARAMETERS]
Usage: Ratatosk --help
Usage: Ratatosk --version
Usage: Ratatosk --cite

[COMMAND]: index

[PARAMETERS]:

   > Mandatory with required argument:

   -s, --in-short                  Input short read file to correct (FASTA/FASTQ possibly gzipped)
                                   List of input short read files to correct (one file per line)
   -l, --in-long                   Input long read file to correct (FASTA/FASTQ possibly gzipped)
                                   List of input long read files to correct (one file per line)
   -o, --out-long                  Prefix of the output index file

   > Mandatory with no argument:

   -1, --1st-pass-only             Prepare index for the 1st correction pass (default: false)
   -2, --2nd-pass-only             Prepare index for the the 2nd correction pass (default: false)

   > Optional with required argument:

   -c, --cores                     Number of cores (default: 1)
   -S, --subsampling               Rate of short reads subsampling (default: No subsampling)
   -u, --in-unmapped-short         Input read file of the unmapped short reads (FASTA/FASTQ possibly gzipped)
                                   List of input read files of the unmapped short reads (one file per line)
   -a, --in-accurate-long          Input high quality long read file (FASTA/FASTQ possibly gzipped)
                                   List of input high quality long read files (one file per line)
                                   (Those reads are NOT corrected but assist the correction of reads in input)
   -g, --in-graph                  Load graph file prepared with the index command
   -Q, --max-base-qual             Maximum base quality of input long reads (default: 40)

   > Optional with no argument:

   -v, --verbose                   Print information

[ADVANCED PARAMETERS]:

   > Optional with required argument:

   -M, --min-conf-color2           Minimum confidence threshold to color vertices for 2nd pass (default: 0)
   -C, --min-len-color2            Minimum length of a long read to color vertices for 2nd pass (default: 3000)
   -i, --insert-sz                 Insert size of the input paired-end short reads (default: 500)
   -k, --k1                        Length of short k-mers for 1st pass (default: 31)
   -K, --k2                        Length of long k-mers for 2nd pass (default: 63)

   > Mandatory with no argument:

   -F, --no-snp-correction         Disable SNP detection and correction
   -I, --no-graph-index            Disable graph index output
```

## FAQ

**Can I provide multiple read files in input?**

Yes, use mutiple times parameters `-s`/`-l`/`-u`/`-a`, once for each input file.

**Can I provide a file which is a list of read files in input?**

Yes, a text file containing one input filename per line with no empty lines can be given in input.

**Does Ratatosk work with input short reads which are not paired-end reads?**

Yes, although Ratatosk works best with paired-end short reads in input. You can also mix paired and non-paired short reads in input.

**Are the output corrected long reads in the same order as the input uncorrected long reads?**

Not necessarily.

**Is it fine if my input *long* reads contain non-{A,C,G,T} characters?**

Yes.

**Is it fine if my input *short* reads contain non-{A,C,G,T} characters?**

Yes (the *k*-mers overlapping these characters will be discarded).

**Why does Ratatosk outputs non-{A,C,G,T} characters?**

Ratatosk automatically detects heterozygous SNP candidates from the input short and long reads. If Ratatosk corrects a SNP candidate site in a long read but cannot establish with enough confidence which of the possible bases is the right one, a IUPAC character is output to represent the ambiguity. For example, character `N` represents all possible bases while character `R` only represents `A` or `G`.

**My downstream pipeline does not handle the non-{A,C,G,T} characters output be Ratatosk, what should I do?**

Replace these characters by random {A,C,G,T} characters. The python script `scripts/replaceIUPAC.py` is provided to that effect.

**What are the quality scores output by Ratatosk?**

Ratatosk outputs a quality score for each corrected base indicating how confident is Ratatosk in the correction of that base. A score of 0 means the base is left uncorrected, a score of 1 means the base was corrected with a very low confidence while a score of 40 means that Ratatosk is very confident of the correction. The scores are output in the Phred33 format but the scoring scale is linear (rather than logarithmic for Phred33).

## Troubleshooting

The following might happen when environment variables are not set correctly on your system:

* compilation (`make`) fails because some header files (*.h*) are not found

Assuming the header files (*.h*) are located at the path */usr/local/include/*, the following command set the environment variables *C_INCLUDE_PATH* and *CPLUS_INCLUDE_PATH* correctly for the time of the session:
```
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/local/include/
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/include/
```

## Citation

```
@article{holley2021ratatosk,
  title={Ratatosk: hybrid error correction of long reads enables accurate variant calling and assembly},
  author={Holley, Guillaume and Beyter, Doruk and Ingimundardottir, Helga and M{\o}ller, Peter L and Kristmundsdottir, Sn{\ae}dis and Eggertsson, Hannes P and Halldorsson, Bjarni V},
  journal={Genome Biology},
  volume={22},
  number={1},
  pages={1--22},
  year={2021}
}
```

## Contact

For any question, feedback or problem, please feel free to file an issue on this GitHub repository and we will get back to you as soon as possible.

## License

* The wyhash library is Unlicense licensed (https://github.com/wangyi-fudan/wyhash)
* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)
* The libdivide library is zlib licensed (https://github.com/ridiculousfish/libdivide)
* The kseq library is copyrighted by Heng Li and released under the MIT license (http://lh3lh3.users.sourceforge.net/kseq.shtml)
* The CRoaring library is Apache 2.0 licensed (https://github.com/RoaringBitmap/CRoaring)
* The GetRSS library is Creative Commons Attribution 3.0 licensed
* The edlib library is MIT licensed (https://github.com/Martinsos/edlib)
* The Bifrost library is BSD2 licensed (https://github.com/pmelsted/bifrost)
