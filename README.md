

NOTES
===========

This is a standalone repository for the local command line version of BART (v2.0). 

For the source code of BART web interface, please refer to:
- <a href="https://github.com/zanglab/BARTweb_frontend">**BARTweb Frontend**</a> for receiving users’ job submission requests and displaying job execution information and results. 
- <a href="https://github.com/zanglab/BARTweb_backend">**BARTweb Backend**</a> for performing all of the computation.

For the custom code for processing the data and generating the figures in the manuscript of **BARTweb: a web server for transcriptional regulator association analysis**, please refer to:
- <a href="https://github.com/zanglab/BARTweb_custom_code">**BARTweb custom code**</a>




README for BART(v2.0)
===========


Introduction
============

BART (Binding Analysis for Regulation of Transcription) is a bioinformatics tool for predicting functional transcriptional regulators (TRs) that bind at genomic cis-regulatory regions to regulate gene expression in the human or mouse genomes, taking a query gene set, a ChIP-seq dataset or a scored genomic region set as input. BART leverages over 7,000 human TR binding profiles and over 5,000 mouse TR binding profiles from the public domain (collected in Cistrome Data Browser) to make the prediction.

BART is implemented in Python and distributed as an open-source package along with necessary data libraries.

BART web interface (Beta version) can be accessed <a href="http://bartweb.org/">here</a>.

**BART is developed and maintained by the <a href="https://faculty.virginia.edu/zanglab/">Chongzhi Zang Lab</a> at the University of Virginia.**



# Installation
#### Prerequisites

BART uses Python's distutils tools for source installation. Before installing BART, please make sure Python3 and the following python packages are installed. We highly recommend the <a href="https://docs.anaconda.com/anaconda/install/">Anaconda environment</a>, which include all the required python packages.

- setuptools
- numpy
- pandas
- scipy
- tables
- scikit-learn
- matplotlib

#### Download the source package and setup the configuration file

You have to download the Human or Mouse Data Library under your own directory before install BART. The unpacked libraries occupy 14GB hard drive storage in the download directory. 

```shell
wget https://virginia.box.com/shared/static/2kqczz9gixetcr9p4bl650uyrio5zd33.gz -O hg38_library.tar.gz
tar zxf hg38_library.tar.gz 
wget https://virginia.box.com/shared/static/bxdggnhp4bjz2l5h2zjlisnzp0ac7axf.gz -O mm10_library.tar.gz
tar zxf mm10_library.tar.gz
```

To install a source distribution of BART, unpack the distribution tarball and go to the directory where you unpacked BART.

```shell
wget https://virginia.box.com/shared/static/jvwc097d7ca0oart6ka4ivonqi4qxkox.zip -O bart21.zip
gunzip bart21.zip
cd bart2-master
```

Modify the configure file (bart2/bart.conf). For example, if you have the hg38_library (or mm10_library) downloaded in this directory: /abc/def/hg38_library (or /abc/def/mm10_library), then the bart.conf file should read:

```shell
[path]
hg38_library_dir = /abc/def/
mm10_library_dir = /abc/def/
```

#### Global installation 
Install with root/administrator permission, or you have the <a href="https://docs.anaconda.com/anaconda/install/">Anaconda environment</a> prepared. By default, the script will install python library and executable codes globally.

```shell
python setup.py install
```

#### Local installation 
If you want to install everything under a specific directory, for example, a directory as /path/to/bart2/, use the following commands.

```shell
mkdir -p /path/to/bart/lib/pythonX.Y/site-packages 
export PYTHONPATH=/path/to/bart/lib/pythonX.Y/site-packages/:$PYTHONPATH 
python setup.py install --prefix /path/to/bart 
export PATH=/path/to/bart/bin/:$PATH
```

In this value, X.Y stands for the major–minor version of Python you are using (such as 3.5 ; you can find this with sys.version[:3] from a Python command line).

You’ll need to modify the environment variables and add those lines in your bash file (varies on each platform, usually is ~/.bashrc or ~/.bash_profile).

```shell
export PYTHONPATH= "/path/to/bart/lib/pythonX.Y/site-packages/:$PYTHONPATH"
export PATH="/path/to/bart/bin/:$PATH"
```

# Tutorial
#### Positional arguments 
`{geneset, profile, region}`

#### bart geneset

Given a query gene set in official gene symbols (HGNC for human or MGI for mouse) in text format (each gene in a row, at least 100 genes recommended), predict functional TRs that regulate these genes.

**Usage**:	`bart2 geneset -i genelist.txt -s hg38 --outdir bart2_output`


#### bart profile

Given a ChIP-seq data file (mapped reads in 
<a href="http://samtools.github.io/hts-specs/SAMv1.pdf" target="_blank">BAM</a> 
or 
<a href="https://genome.ucsc.edu/FAQ/FAQformat#format1" target="_blank">BED</a> 
format in either hg38 or mm10), predict TRs whose binding pattern associates with the input ChIP-seq profile.

**Usage**: 	`bart2 profile -i ChIP.bam -f bam -s hg38 --outdir bart2_output`


#### bart region

Given a scored genomic region set (<a href="https://genome.ucsc.edu/FAQ/FAQformat#format1" target="_blank">BED</a> format
in either hg38 or mm10), predict TRs enriched in this genomic region set.

**Usage**: 	`bart2 region -i ChIPpeak.bed -c 4 -s hg38 --outdir bart2_output`

#### additional argument

`--binsize`:	Size of bin used to estimate ROC AUC. Default is 1000 for geneset and profile mode, 50 for region mode. Setting it lower will guarantees a more precise calculation, with greater time cost.

#### Output files

1. **\*_adaptive_lass_Info.txt** 
    provides regression information tells which representative H3K27ac samples are selected along with coefficients through adaptive lasso regression and sample annotations including cell line, cell type or tissue type. 
    This is the output only generated in geneset mode.
3. **\*_auc.txt** 
provides the association score of each of the TR ChIP-seq dataset with the genome cis-regulatory profile.
4. **\*_bart_results.txt** 
is a rank of all TRs with multiple quantification scores.  

Please refer 
<a href="http://bartweb.org/result?user_key=sample_15881335954485407" target="_blank">here</a>
for an example of BART2 results.


# Citation

If you use BART in your data analysis, please cite: 

<a href="https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty194/4956015" target="_blank">BART: a transcription factor prediction tool with query gene sets or epigenomic profiles</a> <br>
Zhenjia Wang, Mete Civelek, Clint Miller, Nathan Sheffield, Michael J. Guertin, Chongzhi Zang. <i><b>Bioinformatics</b></i> 34, 2867–2869 (2018)


<!---
If you use "geneset" mode, please also cite:
<a href="http://genome.cshlp.org/content/26/10/1417" target="_blank">Modeling cis-regulation with a compendium of genome-wide histone H3K27ac profiles</a> <br>
Su Wang, Chongzhi Zang, Tengfei Xiao, Jingyu Fan, Shenglin Mei, Qian Qin, Qiu Wu, Xujuan Li, Kexin Xu, Housheng Hansen He, Myles Brown, Clifford A. Meyer, X. Shirley Liu. <i><b>Genome Research</b></i> 26, 1417–1429 (2016)
-->

