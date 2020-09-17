# Ratatosk

### Phased hybrid error correction of long reads using colored de Bruijn graphs

Ratatosk is a phased error correction tool for erroneous long reads based on compacted and colored de Bruijn graphs built from accurate short reads. Short and long reads color paths in the graph while vertices are annotated with candidate *de novo* Single Nucleotide Polymorphisms. Long reads are subsequently anchored on the graph using exact and inexact *k*-mer matches to find paths corresponding to corrected sequences.
We demonstrate that Ratatosk can reduce the raw error rate of Oxford Nanopore reads 6-fold on average with a median error rate as low as 0.28%. Ratatosk corrected data maintain nearly 99% accurate SNP calls and increase indel call accuracy by up to about 40% compared to the raw data. An assembly of the Ashkenazi individual HG002 created from Ratatosk corrected ONT reads yields a contig N50 of 43.22 Mbp and less misassemblies than an assembly created from PacBio HiFi reads.

## Table of Contents

* [Requirements](#requirements)
* [Installation](#installation)
* [Interface](#interface)
* [Usage](#usage)
* [FAQ](#faq)
* [Troubleshooting](#troubleshooting)
* [Citation](#citation)
* [Contact](#contact)
* [License](#license)

## Requirements

* C++11 compiler:
    * [GCC](https://gcc.gnu.org/) >= 4.8.5
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

## Installation

1. Clone the Git repository
  ```
  git clone --recursive https://github.com/GuillaumeHolley/Ratatosk.git
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
- Several temporary files are written to disk. These files have the same prefix name as the output file (`-o`) but are deleted at the end of Ratatosk execution. Given an input long read file (`-l`) of size *X* GB, ensure that the output folder has at least about *2.5X* GB of free space.


### Command line
```
Ratatosk --help
```

displays the command line interface:
```
Ratatosk x.y.z

Hybrid error correction of long reads using colored de Bruijn graphs

Usage: Ratatosk [PARAMETERS]

[PARAMETERS]:

   > Mandatory with required argument:

   -s, --in-short                  Input short read file to correct (FASTA/FASTQ possibly gzipped)
                                   List of input short read files to correct (one file per line)
   -l, --in-long                   Input long read file to correct (FASTA/FASTQ possibly gzipped)
                                   List of input long read files to correct (one file per line)
   -o, --out-long                  Output corrected long read file

   > Optional with required argument:

   -c, --cores                     Number of cores (default: 1)
   -q, --quality                   Output Quality Scores: corrected bases get QS >= t (default: no QS)
   -t, --trim-split                Trim and split bases with quality score < t (default: no trim/split)
                                   Only sub-read with length >= 63 are output if used
   -u, --in-unmapped-short         Input read file of the unmapped short reads (FASTA/FASTQ possibly gzipped)
                                   List of input read files of the unmapped short reads (one file per line)
   -a, --in-accurate-long          Input high quality long read file (FASTA/FASTQ possibly gzipped)
                                   List of input high quality long read files (one file per line)
                                   (Those reads are NOT corrected but assist the correction of reads in input).

   > Optional with no argument:

   -v, --verbose                   Print information during execution

[ADVANCED PARAMETERS]:

   > Optional with required argument:

   -i, --insert_sz                 Insert size of the input short reads (default: 500)
   -k, --k1                        Length of short k-mers for 1st pass correction (default: 31)
   -K, --k2                        Length of long k-mers for 2nd pass correction (default: 63)
   -w, --max_len_weak1             Do not correct weak regions >= w bases during 1st pass correction (default: 1000)
   -W, --max_len_weak2             Do not correct weak regions >= w bases during 2nd pass correction (default: 10000)
```

### ***de novo*** correction

```
Ratatosk -v -c 16 -s short_reads.fastq -l long_reads.fastq -o out_long_reads
```
Ratatosk corrects (`Ratatosk`) the long read file (`-l long_reads.fastq`) with 16 threads (`-c 16`) using an index built from the short read file (`-s short_reads.fastq`). Information messages are printed during the execution (`-v`) and the corrected long reads are written to file *out_long_reads_corrected* (`-o out_long_reads`).

### Reference-guided correction

See [reference-guided preprocessing](https://lsource2.decode.is/stat/ratatosk/tree/master/scripts/reference_guiding).

## Options

- **Insert size** (`-i`)

  By default, Ratatosk considers that the input paired-end short reads insert size is roughly 500 (read length * 2 + fragment size). You can set the exact insert size with `-i`. If your input short reads are single end reads, set the insert size to the read length.

- **Quality scores** (`-q`)

  By default, Ratatosk outputs the corrected long reads without quality scores in FASTA format. By using `-q`, corrected reads are output with quality scores in FASTQ format. The output quality score of a base reflects how confident is Ratatosk in the correction of that base. Given a minimum quality score *X* (`-q X`), all corrected bases have a quality score ranging from *X* to 40 while uncorrected bases have a quality score of 0. Note that Ratatosk uses Phred33 with scores ranging from 0 to 40.

- **Trimming/Splitting** (`-t`)

  By default, Ratatosk outputs all bases (corrected and uncorrected). By using `-t`, bases with a low correction quality score are trimmed and split. Specifically, given a minimum quality score *X* (`-t X`), only subsequences of the corrected long reads for which the bases have a correction quality score equal to or larger than *X* are output. Each output subsequence will have as name `>name/i` (FASTA) or `@name/i` (FASTQ) where `name` is the input name of the long read and `i` is an integer subsequence ID for read `name`. Note that only subsequences larger than the *k*-mer size in Ratatosk (63) are output.

## Advanced usage

The default *k1*/*k2*-mer lengths (1st/2nd correction passes) are 31/63. To work with larger *k*-mers (using `-k/-K`), you must compile Ratatosk with a larger `MAX_KMER_SIZE` parameter where `MAX_KMER_SIZE=round(k2 + 1, 32)`, i.e. (*k2* + 1) rounded to the larger multiple of 32. Specifying `MAX_KMER_SIZE` at compilation is done as follows when entering the `cmake` command:
```
cmake -DMAX_KMER_SIZE=96 ..
```
In this example, the maximum *k1*/*k2*-mer length allowed is 95.

## FAQ

**Can I provide multiple read files in input?**

Yes, use mutiple times parameters `-i`/`-l`/`-u`/`-p`, once for each input file.

**Can I provide a file which is a list of read files in input?**

Yes, a text file containing one input filename per line with no empty lines can be given in input.

**Does Ratatosk work with input short reads which are not paired-end reads?**

Yes, although Ratatosk works best with paired-end short reads in input. You can mix paired-end and non-paired-end reads in input as well.

**Are the output corrected long reads in the same order as the input uncorrected long reads?**

Yes if Ratatosk was run with a single thread, otherwise the output order is random.

**Is it fine if my input reads contain non-{A,C,G,T} characters?**

Yes but they won't be corrected nor used in the graph index.

## Troubleshooting

The following might happen when environment variables are not set correctly on your system:

* compilation (`make`) fails because some header files (*.h*) are not found

Assuming the header files (*.h*) are located at the path */usr/local/include/*, the following command set the environment variables *C_INCLUDE_PATH* and *CPLUS_INCLUDE_PATH* correctly for the time of the session:
```
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/local/include/
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/include/
```

## Citation

Coming soon

## Contact

For any question, feedback or problem, please feel free to file an issue on this GitHub repository and we will get back to you as soon as possible.

## License

* The xxHash library is BSD licensed (https://github.com/Cyan4973/xxHash)
* The popcount library is BSD licensed (https://github.com/kimwalisch/libpopcnt)
* The libdivide library is zlib licensed (https://github.com/ridiculousfish/libdivide)
* The kseq library is copyrighted by Heng Li and released under the MIT license (http://lh3lh3.users.sourceforge.net/kseq.shtml)
* The CRoaring library is Apache 2.0 licensed (https://github.com/RoaringBitmap/CRoaring)
* The GetRSS library is Creative Commons Attribution 3.0 licensed
* The edlib library is MIT licensed (https://github.com/Martinsos/edlib)
* The Bifrost library is BSD2 licensed (https://github.com/pmelsted/bifrost)
