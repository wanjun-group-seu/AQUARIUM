# design convention #

## top-level interface ##
the name convention of top level pipeline/workflow should begin with "wf_*"

### what should top level interface do ? ###
it should

    1. load the config files
    2. read the additional arguments
    3. sub-level of these module should be tool sets divide by usage. such as align tools, quantification tools, DE-detectors. etc. the final tool should be decide by the sub-level tool sets.
    4. top-level module should provide output check and input check for the sub-level module, that is pre-module function and post-module function.

### input arguments or input file  ###


#### bunch mode ####
hmm.. , a often seen job is so-called "bunch", that is , a bunch of reads file with the same index .

so we need *a input file* , with contains the locations of these reads, and *a setting file*, which stands for the setting for all reads in this bunch .

#### option checker ####
the config file remains mostly all the original cli options. 
a controller/checker should be established, and should provide following feature: 

    1.check the options before invoking external tools, show description message when error occurs
    2.give some query/info interface

every module which contains a option_checker should have this constant *OPTION_CHECKERS* , which is a list of opt-checker objects 


### output / logging part ###
yes , the output file is the logging file, and should use the standard logging module. ... but still undone yet . 
it is mostly the python logging manual's fault .....

### option checker and helper ###

a mechanism for checking opts and print help is needed . 

## QC tools ##
here, we assume the QC tool is fastqc.

and those qc tool should have following parts:

1. run_qc

> a function/method to invoke the QC pipeline


2. parse_qc_result

> a function to parse QC result

3. SECTION_QC_SETTING

> the name of QC setting section in configure file

4. OPTION_CHECKERS

> a list of option check , this will generate the default description string

## aligner ##

most aligner need some kind of index to run the mapping.

so all aligner wrapper should have following feature:

1. index

> function interface for establishing index

	input:
	output: file_path of indice, which can be understand by the aligner itself
if the index already exists, just return the path.

2. align

> function interface for mapping sequence reads to reference genome.

	input: sequence reads in fasta/fastq
	output: return a dict contains serval entries ,

3. is_path_contain_index

> function for checking whether the given path have a index
this is done by checking existence of related files.
those related files should be provide by aligner module itself.

	input : file path to possible index
	output : bool value indiacting whether there is a index .


4. SECTION_INDEX:

> a string shows index parameter is under configure[SECTION_INDEX]

5. SECTION_ALIGN:

> a string shows align parameters are under configure[SECTION_ALIGN]

6. interpret_index_path:

> function : translate a given index path to a dict which is compatible with configure file section.

	input: a string containing possible index path

	output: a dict , which can be accept by ChainMap to make a dict as input of index/align

	when input is None, return {}

7. interpret_seq_files

> function: translate string contains a path or paths of sequencing read file to a dict which is compatible with configure file section

	input: a string contains one or more sequencing reads file path
	output: a dict , which can be accept by ChainMap to make a dict as input of align function

	if input is None, return {}


8. is_map_result_already_exists
> function for checking whether the alignment result already exists
    
    input : file path, or path/prefix specify the path
    output: bool value shows whether there already have mapped result . 

9. get_align_result_path
> function to get the path of mapped result file 
    
    input : the same input as the align function 
    output: path or path/prefix showing where the output file should be 


### suggested paradigm ###

    def index/quant(config, ...)
        get_dict(...)
        check_opts
        get_cmd
        run_cmd

## quantification tools ##

similar with aligners, quantification tools also need index . so the main framework of quantification tools is the same as the aligners , only the spell differs .

1. index

> function interface for establishing index

	input: ... some massive dict
	output: return the file_path of indice

2. quantify

> function interface for calculating the expression

	input: sequence reads in fasta/fastq
	output: return a dict including the path of expression file

3. is_path_contain_index :

> function for checking whether the given path have index

	input : file path to possible index
	output : bool value indiacting whether there is a index .


4. SECTION_INDEX:

> a string shows index parameter is under configure[SECTION_INDEX]

5. SECTION_QUANTIFY:

> a string shows align parameters are under configure[SECTION_QUSANTIFY]

6. interpret_index_path:

> function : translate a given index path to a dict which is compatible with configure file section.

	input: a string containing possible index path

	output: a dict , which can be accept by ChainMap to make a dict as input of index/quant
	return {} if None is given



7. interpret_seq_files

> function: translate string contains a path or paths of sequencing read file to a dict which is compatible with configure file section

	input: a string contains one or more sequencing reads file path
	output: a dict , which can be accept by ChainMap to make a dict as input of quant function
	return {} if None is given






## Circular RNA detecting tools ##

may need some sequence reads as input , this will be aligner's output .
seems this tool do not need a index phrase . or it can be done manually.

1. detect

> main function interface ,

    input: sequence_file, path to sequence reads fils in fastq/fa format
    output: a dict(or other container data type), this should contain a path will the report of circular RNA detection

2. export_as_bed

> function to transform detection report to .bed file type

    input: a dict contains settings for this dectection job. 
    output: 1. a bed file contains the information of circRNA detection.
            2. optional gene-mapping file

3. SECTION_DETECT:

> a string shows where is the setting of this detection tool in config file.

4. interpret_seq_files

> function designed for work with files assigned by arguments,

5. which_external_aligner

> function to tell whether a external mapping should be performed before the detecting phase . 
> return "" if no aligner is needed . 
> else return the config section name of aligner such as  "BWA", "STAR", etc

6. bed_out

> [optional] file path to bed type result, summarized by module's to_bed function

### bed file format convention ###
bed file of circRNA are comparable with the default bed file format described in http://www.ensembl.org/info/website/upload/bed.html

meaning of each column are :
first 3 column are required field.
1. chrom  : name of chromosome
2. chromStart : start position of the circRNA
3. chromEnd : end position of the circRNA

4. name : should follow pattern : rna_name@gene_name_or_id, string after the first @ symbol will be taken as description,
5. score : number of junction reads detected by detection tools


## DE tool ##


# TODO:  need give lazy option for aligners 
# TODO: detector module should have a flag show whether it need an external aligner

