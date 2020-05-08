# Ratatosk

### Phased hybrid error correction of long reads using colored de Bruijn graphs

Ratatosk is a phased error correction tool for erroneous long reads based on compacted and colored de Bruijn graphs built from accurate short reads. Short and long reads color paths in the graph while vertices are annotated with candidate *de novo* Single Nucleotide Polymorphisms. Long reads are subsequently anchored on the graph using exact and inexact *k*-mer matches to find paths corresponding to corrected sequences.
We demonstrate that Ratatosk can reduce the raw error rate of long reads down to X\,\% while preserving read phasing. Ratatosk corrected data allow for the same high quality SNP calls as for the raw data but indel calling precision and recall substantially increase by up to Y\,\% in the corrected data. Phased assemblies created from Ratatosk corrected data are also significantly more contiguous compared to the raw long reads.

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
5. Use the opened Ubuntu terminal for compiling, installing and running Bifrost (see next section). See [Troubleshooting](#troubleshooting) if you have any problem during the installation.

## Installation

  ```
  git clone --recursive https://github.com/GuillaumeHolley/Ratatosk.git
  cd Ratatosk && mkdir build && cd build
  cmake ..
  make
  make install
  ```

  `make install` might requires `sudo` (`sudo make install`) to proceed. If you want to install Ratatosk in a non-default path, replace `cmake ..` with `cmake -DCMAKE_INSTALL_PREFIX=/my/path/ ..` where `/my/path/` is where you want to see the Ratatosk binary installed. Do not forget to had this path to your environment variables (see [Troubleshooting](#troubleshooting)). If you encounter any problem during the installation, see the [Troubleshooting](#troubleshooting) section.

  By default, the installation creates a binary (*Ratatosk*)

## Usage:

```
Ratatosk
```

displays the command line interface:
```
Ratatosk x.y

Phased hybrid error correction of long reads using colored de Bruijn graphs

Usage: Ratatosk [PARAMETERS]

[PARAMETERS]:

   > Mandatory with required argument:

   -i, --in-short-files        Input short read files (FASTA/FASTQ possibly gzipped)
                               Input short read files can be provided as a list in a TXT file (one file per line)
   -l, --in-long-files         Input long read files (FASTA/FASTQ possibly gzipped)
                               Input long read files can be provided as a list in a TXT file (one file per line)
   -o, --out-file              Output corrected long read file

   > Optional with required argument:

   -c, --cores                 Number of cores (default is 1)
   -t, --trimming              Trim bases with Q-score < t with 0 <= t <= 40 (default is t=0, no trimming)
                               Note that only sub-read with length >= k are output if t > 0
   -q, --quality               Output quality scores such that corrected bases have Q-score >= t (default is t=0, no output)
   -u, --in-unmap-short-files  Input unmapped short read files (FASTA/FASTQ possibly gzipped)
                               Input unmapped short read files can be provided as a list in a TXT file (one file per line)
   -p, --in-helper-long-files  Input high quality long read files (FASTA/FASTQ possibly gzipped)
                               Input high quality long read files can be provided as a list in a TXT file (one file per line)
                               Those reads are *not* corrected but help the correction.
   -m, --in-unmap-graph-file   Input graph file of unmapped reads (default is no input graph)
   -g, --in-graph-file         Input graph file (default is no input graph)
   -w, --out-graph-file        Output graph file (default is no output graph)

   > Optional with no argument:

   -v, --verbose            Print information messages during execution
```

## ***de novo*** correction

```
Ratatosk -v -c 16 -i short_reads.fastq -l long_reads.fastq -o corrected_long_reads
```
Ratatosk corrects the long read file (`-l long_reads.fastq`) with 16 threads (`-c 16`) using an index built from the short read file (`-i short_reads.fastq`). Information messages are printed during the execution (`-v`) and the corrected long reads are written to file *corrected_long_reads.*

## Reference-guided correction

While Ratatosk is a reference-free method, a reference-guided preprocessing of the input data is proposed to improve the correction, redcue the running time and distribute the workload over many nodes of an HPC. First, short and long reads with good mapping quality are distributed over bins corresponding to segments of the reference genome they map to. Bins are corrected independently. Then, unmapped long reads and long reads with low mapping quality are corrected *de novo* using all short reads and the previously corrected long read bins. More details are provided in.

1. Map the long and short reads. We advise [bwa-mem](https://github.com/lh3/bwa) to map the short reads and [minimap2](https://github.com/lh3/minimap2) for the long reads. BAM files must be sorted and indexed (`samtools sort` and `samtools index`).

2. Bin the reads:
```
```


## Notes

- Ratatosk works best with paired-end short reads in input (`-i`): reads from the same pair **must** have the same FASTA/FASTQ name (if the reads are extracted from a BAM file, use `samtools bam2fq -n`). The order of the reads in the input file(s) does not matter.

- Several temporary files are written to disk. Those files have the same prefix name as the output file (`-o`) but are deleted at the end of Ratatosk execution. Given an input long read file (`-l`) of size *L* GB, ensure that the output folder has at least about *2.5L* GB of free space.

## FAQ

**Can I provide multiple read files in input?**

Yes, files must be separated by a space character for parameter `-i` and `-l`.

**Can I provide a file which is a list of read files in input?**

Yes, a text file containing one input filename per line with no empty lines can be given in input.

**Does Ratatosk work with input short reads which are not paired-end reads?**

Yes, although Ratatosk works best with input paired-end short reads. You can mix paired-end and non-paired-end reads in input as well.

**Are the output corrected long reads in the same order as the uncorrected input long reads?**

Yes if Ratatosk was run with a single thread, no otherwise.

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