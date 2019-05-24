
README for BART(1.0.3)

Introduction
============

BART (Binding Analysis for Regulation of Transcription) is a bioinformatics tool for predicting functional transcription factors (TFs) that bind at genomic cis-regulatory regions to regulate gene expression in the human or mouse genomes, given a query gene set or a ChIP-seq dataset as input. BART leverages over 7,000 human TF binding profiles and over 5,000 mouse TF binding profiles from the public domain (collected in Cistrome Data Browser) to make the prediction.

BART is implemented in Python and distributed as an open-source package along with necessary data libraries.

**BART is developed and maintained by the Chongzhi Zang Lab at the University of Virginia.**



# Installation
#### Prerequisites

BART uses Python's distutils tools for source installation. Before installing BART, please make sure Python3 (Python 3.3 or higher is recommended) is installed in the system, and the following python packages are installed:

- setuptools
- numpy
- pandas
- scipy

#### Install the full package (All data included, requires at least 15GB hard drive storage in the installation directory)

To install a source distribution of BART, unpack the distribution tarball and open up a command terminal. Go to the directory where you unpacked BART, and simply run the install script an install BART globally or locally. For example, if you want to install the package BART-v1.0.2-py3-full.tar.gz:

```shell
$ tar zxf BART-v1.0.2-py3-full.tar.gz
$ cd BART-v1.0.2-py3-full
```

Install with root/administrator permission (by default, the script will install python library and executable codes globally):

```shell
$ python setup.py install
```

If you want to install everything under your own directory, for example, a directory as /path/to/bart/, use these commands:

```shell
$ mkdir -p /path/to/bart/lib/pythonX.Y/site-packages 
$ export PYTHONPATH=/path/to/bart/lib/pythonX.Y/site-packages/:$PYTHONPATH 
$ python setup.py install --prefix /path/to/bart 
$ export PATH=/path/to/bart/bin/:$PATH
```

In this value, X.Y stands for the major–minor version of Python you are using (such as 3.5 ; you can find this with sys.version[:3] from a Python command line).

Configure environment variables

You’ll need to add those two lines in your bash file (varies on each platform, usually is `~/.bashrc` or `~/.bash_profile`) so that you can use the BART command line directly:

```shell
$ export PYTHONPATH= "/path/to/bart/lib/pythonX.Y/site-packages/:$PYTHONPATH"
$ export PATH="/path/to/bart/bin/:$PATH"
```



#### Install from source package without data libraries (recommended)

You can download the Human or Mouse Data Library separately under your own directory. In this case, you have to edit the config file (e.g. BART1.0.1/BART/bart.conf) after you unpack the source package to provide the directory for the data. For example, if you download the hg38_library.tar.gz (or mm10_library.tar.gz) and unpack it under /path/to/library, then you can modify the bart.conf file as:

`hg38_library_dir = /path/to/library/`

Then you can run the install script and install BART source package globally or locally same as the full package described above.



# Tutorial
Positional arguments {geneset,profile}


#### bart geneset

Given a query gene set (at least 100 genes recommended), predict functional transcription factors that regulate these genes.

**Usage**:	`bart geneset [-h] -i <file> [--refseq] -s <species> [-t <target>] [-p <processes>] [--nonorm] [--outdir <outdir>] [-o <ofilename>]`

**Example**:	bart geneset -i name_enhancer_prediction.txt -s hg38 -t target.txt -p 4 --outdir bart_output

**Input arguments**:

`-i <file>, --infile <file>`

Input file contains gene symbol in each row.

`--refseq` 

Input file contains 2 columns, one is geneID(refseqID) and the other is 1/0. 1 for traget and 0 for un-target. 

`-s <species>, --species <species>`

Species, please choose from "hg38" or "mm10".

`-t <target>, --target <target>`

Target transcription factors of interests, please put each TF in one line. BART will generate extra plots showing prediction results for each TF.

`-p <processes>, --processes <processes>`

Number of CPUs BART can use.

`--nonorm`

Whether or not do the standardization for each TF by all of its Wilcoxon statistic scores in our compendium. If set, BART will not do the normalization. Default: FALSE.

**Output arguments**:

`--outdir <outdir>`

If specified, all output files will be written to that directory. Default: the current working directory

`-o <ofilename>, --ofilename <ofilename>`

Name string of output files. Default: the base name of the input file.

**Notes**:

The input file for \<BART geneset\> :

a. Gene symbol version:

|  AR   |
| :---: |
| BAGE4 |
| CST1  |
| DLX3  |

b. Refseq version:

| NM_000397 |  1   |
| :-------: | :--: |
| NR_045675 |  0   |
| NM_033518 |  1   |
| NM_052939 |  1   |
| NM_002290 |  0   |
| NM_000495 |  0   |



#### bart profile

Given a ChIP-seq data file (bed or bam format mapped reads), predict transcription factors whose binding pattern associates with the input ChIP-seq profile.

**Usage**: 	`bart profile [-h] -i <file> -f <format> [-n <int>] -s <species> [-t <target>] [-p <processes>] [--nonorm][--outdir <outdir>] [-o <ofilename>]`

**Example**:	bart profile -i ChIP.bed -f bed -s hg38 -t target.txt -p 4
				--outdir bart_output

**Input files arguments**:

`-i <file>, --infile <file>`

Input ChIP-seq bed or bam file.

`-f <format>, --format <format>`

Specify "bed" or "bam" format.

`-n <int>, --fragmentsize <int>`

Fragment size of ChIP-seq reads, in bps. Default: 150.

`-s <species>, --species <species>`

Species, please choose from "hg38" or "mm10".

`-t <target>, --target <target>`

Target transcription factors of interests, please put each TF in one line. BART will generate extra plots showing prediction results for each TF.

`-p <processes>, --processes <processes>`

Number of CPUs BART can use.

`--nonorm`

Whether or not do the standardization for each TF by all of its Wilcoxon statistic scores in our compendium. If set, BART will not do the normalization. Default: FALSE.

**Output arguments**:

`--outdir <outdir>`

If specified, all output files will be written to that directory. Default: the current working directory

`-o <ofilename>, --ofilename <ofilename>`

Name string of output files. Default: the base name of input file.

**Notes**:

The input file for \<BART profile\> should be [BED](https://genome.ucsc.edu/FAQ/FAQformat#format1) or [BAM](http://samtools.github.io/hts-specs/SAMv1.pdf) format in either hg38 or mm10. 

Bed is a tab-delimited text file that defines the data lines, and the BED file format is described on [UCSC genome browser website](https://genome.ucsc.edu/FAQ/FAQformat). For BED format input, the first three columns should be chrom, chromStart, chromEnd, and the 6th column of strand information is required by BART. 

BAM is a binary version of [Sequence Alignment/Map(SAM)](http://samtools.sourceforge.net) format, and for more information about BAM custom tracks, please click [here](https://genome.ucsc.edu/goldenPath/help/bam.html). 

#### bart region

Given a non-overlapping bed file with score on 5th column, predict transcription factors whose binding pattern associates with the input profile.

**Usage**: 	`bart region [-h] -i <file> -s <species> [-t <target>] [-p <processes>] [--nonorm][--outdir <outdir>] [-o <ofilename>]`

#### Output files

1. **name_auc.txt** contains the ROC-AUC scores for all TF datasets in human/mouse, we use this score to measure the similarity of TF dataset to cis-regulatory profile, and all TFs are ranked decreasingly by scores. The file should be like this:

|   AR_56206   | AUC = 0.854 |
| :----------: | :---------: |
|   AR_56205   | AUC = 0.846 |
|   AR_69287   | AUC = 0.844 |
|   AR_68090   | AUC = 0.835 |
|  EZH2_41813  | AUC = 0.820 |
| SUZ12_74199  | AUC = 0.819 |
| JARID2_41807 | AUC = 0.819 |
|   AR_51837   | AUC = 0.818 |

2. **name_bart_results.txt** is a ranking list of all TFs, which includes the Wilcoxon statistic score, Wilcoxon p value, standard Wilcoxon statistic score (zscore), maximum ROC-AUC score, rank score (relative rank of z score, p value and max auc) and Irwin-Hall P-value for each TF. The most functional TFs of input data are ranked first. The file should be like this:

   |   TF   | Wilcoxon statistics | Wilcoxon P-value | Z-score | max_AUC | relative_rank | Irwin-Hall P-value |
   | :----: | :-----------------: | :--------------: | :-----: | :-----: | :-----------: | ------------------ |
   |   AR   |       22.917        |    1.572e-116    |  3.228  |  0.854  |     0.004     | 3.126e-07          |
   | POLR3D |        5.770        |    3.972e-09     |  3.038  |  0.582  |     0.021     | 4.125e-05          |
   |  CTCF  |       25.105        |    2.169e-139    |  3.026  |  0.560  |     0.023     | 5.332e-05          |
   |  BRDU  |        4.206        |    1.302e-05     |  2.952  |  0.592  |     0.025     | 7.065e-05          |
   |   PR   |        6.690        |    1.118e-11     |  2.647  |  0.640  |     0.030     | 1.158e-04          |
   | RAD21  |        8.222        |    1.003e-16     |  3.008  |  0.537  |     0.034     | 1.712e-04          |

3. **name_plot** is a folder which contains all the extra plots for the TFs listed in target files (*target.txt* file in test data). For each TF, we have rank dot plot, which shows the rank position and rank significance of this TF in all TFs (derived from the Irwin-Hall P-value score in *name_bart_results.txt*), and the cumulative distribution plot, which compares the distribution of ROC-AUC scores from datasets of this TF and the scores of all datasets (derived from the AUC scores in *name_auc.txt*).
