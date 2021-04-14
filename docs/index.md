# AQUARIUM - Accurate QUAntification of circulaR Isoforms Using Model-based strategy

## Introduction

AQUARIUM (Accurate QUAntification of circulaR Isoforms Using Model-based strategy) is a bioinformatics analysis pipe line for circular RNA sequencing data.

This pipeline can determine the relative abundance of both linear and circRNA in RNA-seq data under the guidance of a configuration file.
Since the relative expression abundunce of circular and linear RNAs are unified under the same standard, it avoids the need to introduce other parameters to describe the expression of cyclic RNA in other count-based circular RNA quantification tools.

It accepts output of circRNA identification tools (CIRI, CIRI-full) or a BED-format file to specify the circular RNA transcripts. Then, it transforms all circular transcripts to pseudo-linear transcripts. Finally, it estimates the expression of both linear and circular transcripts using salmon framework.

Code written in python and R, and published under mit license.

## Installation

### install step 1: download code

This pipeline does not need compilation. You can download the package directly from the [project homepage](https://github.com/wanjun-group-seu/AQUARIUM) as a zip file and unzip it manually.

If [git](https://git-scm.com) is installed, you can also download this pipeline using following command.

``` bash
git clone https://github.com/wanjun-group-seu/AQUARIUM
```


### install step 2: prerequisites

[Python 3](https://www.python.org) and [R](https://www.r-project.org) need to be installed before you can proceed with the following steps.

#### Python packages

This pipeline requires [biopython](https://biopython.org) and [gffutils](https://pythonhosted.org/gffutils/).

But if you are not sure whether you have installed the relevant package or not.
you can use a versatile helper script `info.py`:

```bash
python path/to/your/AQUARIUM/info.py
```

When this script is invoked in command line, it first checks whether the software-related python packages are installed properly.
It will report an error once any python package is not installed properly. You can check the error message to find the missing package.

#### R packages

These R package may be required for analysis: `tidyverse`, `Biostrings` , `rtracklayer`, `GenomicFeatures`, the latter 3 packages are part of  [bioconductor](https://www.bioconductor.org).

Install them in the R console:

```R

# To install `tidyverse`:
install.packages("tidyverse")

# to install other 3 bioconductor packages

install.packages("BiocManager")
BiocManager::install("Biostrings")
BiocManager::install("rtracklayer")
BiocManager::install("GenomicFeatures")

```

If you need to install  `Biostrings` `rtracklayer` `GenomicFeatures` in an older version of R that does not work with `BiocManager`, please check `bioconductor` page of these packages.

#### gffread

gffread is part of Cufflinks (http://cole-trapnell-lab.github.io/cufflinks/). This GFF/GTF parsing utility is used to extract sequence.

### install step 3: adjust configuration template

In order to use this pipeline to analyze the data, a configuration file is essential.
But instead of writing it from scratch, you can call `new_config.py` to create a new one, and then modify it.

The default configure template lies in the `pysrc/body/default.cpysrc/body/default.cfg`.
When you invoke [new_config.py](#get-config-file), you're actually making a copy of this file.

This template is almost empty after a fresh installation.
You can modify this file to save your time in future use.

More details can be found in [instructions on config file](#configure-file)

## CLI interface and Usage

The main command line interface for circRNA profiling is `wf_profile_circRNA.py`.

You can luanch it as follows:

```bash
python wf_profile_circRNA.py config_file
```

The only argument is the path of config file. 
It contains the information needed to analyze a single sample, 
where is the input files , where is the output files, the parameters of certain steps, etc.

Before calling this interface, you can [get a config copy](#get-config-file) and then [modify it](#modify-modules)

## Example

Here we will show how to perform a simple but complete analysis after the installation of pipeline.

### data preparation

You can download a minimal dataset for testing from [MEGA](https://mega.nz/file/ht4mhACS#bSmu5jH_dUu92rK8qoVmpptObydjm89deVsf9cYfK_4).

It contains a tiny genome( mt.fa, mt.gtf),pair end sequence data (r1.fq, r2.fq), and a circRNA detection report(circ_mt.report), 
also a sample config file(your.cfg).

Run this pipeline using that sample config file will cause error, because file path in that cfg file is not correct. And creating a new config file is a more preferable way. But you use it as a guide by searching "TODO" to find what needs to be edited.

### get config file

By invoking new_config.py, you will get a copy of the template. then you can modify this config file to guide your analysis.

The format of invoking new_config.py is as following:

```shell
new_config.py target
```

`target` is the place where your configuration file will be found.

If this `target` is an existing folder, then a `default.cfg` file will appear in this folder.

### modify this configure file

Before you modify the configuration file, it is recommended that you read [instructions on configure file](#configure-file)

If you have modified the template after installation, you may notice that the changes you have made are already in effect on this copy.

With the sample config as a guide, you can modify it by referring to [instructions provided by info.py](#modify-modules).

Here, our goal is to quantify the RNA in the sequencing data. so you can query info.py like this:

```bash
python path/to/AQUARIUM/info.py profile_circRNA
```

From its output, we can see that the relevant module is [CIRC_PROFILE].

### Run the pipeline

After modifying the configuration file, it is time to try running the process.

```bash
python path/to/AQUARIUM/wf_profile_circRNA.py some.cfg
```

It will output a lot of log messages to the screen. So you can use the `nohup` command to make it run in the background, like:

```bash
nohup python path/to/AQUARIUM/wf_profile_circRNA.py some.cfg &
```

You can find a `quant.sf` with RNA expression levels in the folder specified by the `-o` parameter in the `CIRC_PROFILE` section of configration. 

## Troubleshooting

### Examine error messages

When an error occurs, first check the log message.  The error message will be at the end of the screen output. Search the words ERROR and WARING in the LOG message would also help.

If you come across a strange looking error message or find a bug, please do let us know. You submit new issues here: https://github.com/wanjun-group-seu/AQUARIUM/issues

### build the annotation database manually (optional)

During the process, we will create a binary database file for GTF in order to query the genome information faster. This '.db' database file will share base name with the GTF file and will locate in the same folder.

However, since there are different versions of GTF files for different species.
So creating the database may encounter some unexpected problems.

For this reason, a prudent solution is to create the relevant database manually with following steps:

1. Enter the python interactive environment
   
   ```bash
   python
   ```

2. call `gffutils.create_db` to build the database.
   ```python
    import gffutils
    gffutils.create_db(gtf_file, db_file_path)
   ```
    here, gtf_file is path to your gtf file, db_file_path is where the .db file will locate.

In general it takes about 15 minutes to process the human genome GTF, But if it takes too long, or if it shows a warning, some parameters needs to be tweaked. Get more information by

```python
help(gffutils.create_db)
```

## Configure File

We use an `ini` format configuration file to set all parameters in the process.

The user can call `new_config.py` to get a copy of the default configuration file, and modify it to meet specific requirements.

Typically, a single sample corresponds to a single configuration file.

For large-scale analysis of multiple samples, users can start by writing a configuration file template, and then use character substitution to obtain multiple configuration files.

### structure of config

The configuration file is a plain text file that follows the syntax of the `ini` file.

The file is divided into sections, with the name of each section stored in square brackets, and pairs of key-value pairs under the square brackets. It suppports syntax like `${A:b}` to quote the value of b in another section A.

``` ini
[section_name]
keyname1=value1
keyname2=value2
keyname3=${A:b}
```
You might see a lot of sections in this file, but don't worry, only some of them will be used in each  operation. 

These sections can be divided into two categories, first three sections are required for all analysis modules, and the other sections correspond to different module of analysis.

The first three sections meta/global/custom store the paths to the necessary files for the analysis process.

Section `META` contains parameters that are constant in the environment, such as the location of the executable and genome annotations.

Section `GLOBAL` includes information that is invariant in the analysis of the same batch of samples, such as location of sequencing data.

Section `CUSTOM` is used to store the information corresponding to a single sample. and some user-defined variables.

The subsequent sections are related to specific analysis modules, with more information available in `info.py`.

In each section you can see a lot of keys reserved as placeholder. You need to refer to the results of info.py to determine which ones are useful.

### modify modules

To facilitate your modification of the parameters of each module, you need to use `info.py`.

This is a multi-purpose script. During the installation process, you can invoke info.py to check whether the relevant package has been installed or not. After the installation is complete, info.py is also used to query the information for each analysis step.

If the software-related packages are installed correctly, invoking info.py without parameters will give you the following information:

``` bash
python path/to/your/AQUARIUM/info.py

this pipeline now has following wrappers, choose one to see more information

bwa     -        BWA and BWA MEM
ciri    -        CIRI : a circular RNA detection tool
ciri2   -        CIRI2 : a circular RNA detection tool, multiple core empowered
ciri_as         -        CIRI-AS: circular RNA Alternative Splicing Event detection tool
ciri_full       -        CIRI-FULL: a powerful circular RNA detection and rebuilding tool, which combines CIRI and CIRI-AS
detect_circRNA  -        home made pipeline for detecting circRNA
fastqc  -        fastqc , a QC tool
knife   -        KNIFE: a circular RNA detection tool
profile_circRNA         -        home made pipeline for profiling the circRNA
rsem    -        RSEM : a RNA seq quantification tool
sailfish        -        Sailfish: a RNA-seq quantification tool based on k-mer
salmon  -        Salmon: a RNA-seq quantification tool based on fragment
star    -        STAR : a junction sensitive aligner


 you could use 'grep -v "^#\ " | grep -v "^$"' to filter out the options and use them in config file

```

Here, it prompts you to select an item in the list for more information. Each of these items corresponding to a section in the configuration file.  `profile_circRNA` is the main step of the process. It correspond to `wf_profile_linear.py`

Here we look up the parameters of `profile_circRNA`:

``` bash

-> % python path/to/your/AQUARIUM/info.py profile_circRNA

[CIRC_PROFILE]
# this section [CIRC_PROFILE] contains following options:
#######  essential arguments

quantifier = salmon
# the back end quantifier: sailfish or salmon

-o = 
# output folder that contains the index built by sailfish and quantification results

#######  only one of these option pairs are needed

-g/--genomic_seqs_fasta
# path to genomic sequence fasta file

-a/--annotation
# path to gene annotation file, ie, .gtf or .gff files

-c/--ciri_bsj/--bed
# path to  circRNA detection (file for ciri, folder for ciri-full) to specify circular RNA

#######  optional arguments

ciri_as_prefix
# path prefix to CIRI-AS report of circular exons

--decoys
# path to decoy sequences for salmon

preserved_id_list
# path to preserved id list

--mll
# mean library length, this option is to fix up the effective length.

-k = 31
# k-mer size used by sailfish to built index. default is 21

-1
# path to raw pair-end reads, mate 1

-2
# path to raw pair-end reads, mate 2

-r
# path to single-end raw sequencing reads file.

additional_circ_ref
# path to additional circular RNA reference file (.fa),

additional_annotation
# path to additional circular RNA annotation file in gtf format

additional_linear_ref
# path to additional linear RNA reference file(.fa)

flag_use_linc_explicitly
# flag to specify whether linc RNA should be include in quantification result, commenting it out to
# set False

flag_reject_linear
# flag to specify whether to reject linear RNA during quantification, for example for a RNase R
# treated sample, comment it out to set False

#######  only one of these following options are optional




#######  these options are not allowed

# # -h

# # --help
```

As you can see, the output of the command line lists all parameters of the [CIRC_PROFILE] section of the configuration file, along with a brief description of each parameter.

In fact, you can just overwrite the relevant section of the profile with the texts above, and then fill in the key-value pairs afterwards. 

Also, You may notice some commented out key-value pairs. With the help of info.py, you can change them to you need.

***All path should be absolute path***

### FAQ

#### Q: How to complete the parameters in the META section?

A: Not all of these parameters are essential. 

If you just want to use the existing ciri-full result to quantify the RNA-seq data. These paramters are required:
`num_thread`, `genome_fa`,`genome_annotation`,`salmon_bin`. 

#### Q: how to define the library type of salmon ?

A: This depends on the situation you have encountered, You can read the following pages for complete information: [Salmon Library Type](https://salmon.readthedocs.io/en/latest/library_type.html)

#### Q: why file path use absolute path? isn't it to long to use?

A: This make the config file "portable". If you think the absolute paths are too long, and there are many duplicate parts. You can define some prefixes for the paths in custom, and refer to them later