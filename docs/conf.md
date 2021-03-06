
<!-- 这个文档用来解释config文件的结构与内容 -->

# configure file 

We use an `ini` format configuration file to set all parameters in the process.

The user can call `new_config.py` to get a copy of the default configuration file, and modify it to meet specific requirements.

Typically, a single sample corresponds to a single configuration file.

For large-scale analysis of multiple samples, users can start by writing a configuration file template, and then use character substitution to obtain multiple configuration files.

## file structure

The configuration file is a plain text file that follows the syntax of the `ini` file.

The file is divided into sections, with the name of each section stored in square brackets, and pairs of key-value pairs under the square brackets. It suppports syntax like `${A:b}` to quote the value of b in another section A.

``` ini
[section_name]
keyname1=value1
keyname2=value2
keyname3=${A:b}
```
You might see a lot of sections in this file, but don't worry, only some of them will be used in each  operation. These sections can be divided into two categories, the first three sections are required for all types of analysis, and the other sections correspond to different types of analysis.

The first three sections meta/global/custom store the paths to the necessary files for the analysis process.

Section `META` contains parameters that are constant in the environment, such as the location of the executable and genome annotations.

Section `GLOBAL` includes information that is invariant in the analysis of the same batch of samples, such as location of sequencing data.

Section `CUSTOM` is used to store the information corresponding to a single sample. and some user-defined variables.


The subsequent sections are related to specific analysis steps, with more information available in `info.py`.

## modify the configure file with the help of info.py

During the installation process, you can invoke info.py to check whether the relevant package has been installed or not. 

After the installation is complete, info.py is also used to query the information for each analysis step.

If the software-related packages are installed correctly, invoking info.py without parameters will give you the following information:

``` bash
python ~/project/working_pipe_line/info.py

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

Here, it prompts you to select an item in the list for more information. Each of these items corresponding to a section in the configuration file. `detect_circRNA` `profile_circRNA` are the two main steps of the process. They correspond to the scripts wf_detect_circRNA.py
and wf_profile_linear.py

Here we look up the parameters of profile_circRNA:

``` bash

-> % python ~/project/working_pipe_line/info.py detect_circRNA

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

## FAQ of configuration

### Q: How to complete the parameters in the META section?

A: Not all of these parameters are essential. 

If you just want to use the existing ciri-full result to quantify the RNA-seq data. These paramters are required:
`num_thread`, `genome_fa`,`genome_annotation`,`salmon_bin`. 

### Q: how to define the library type of salmon ?

A: This depends on the situation you have encountered, You can read the following pages for complete information: [Salmon Library Type](https://salmon.readthedocs.io/en/latest/library_type.html)

### Q: why file path use absolute path? isn't it to long to use?

A: This make the config file "portable". If you think the absolute paths are too long, and there are many duplicate parts. You can define some prefixes for the paths in custom, and refer to them later