# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
# ChangeLog:
#
# 2020-01-09: the main work-horse moved from sailfish to salmon, need to tweak the interface of it .
# and current salmon version is : 1.1.0 (https://github.com/COMBINE-lab/salmon/releases)
#
#
import copy
import math
import os

import pysrc.body.cli_opts
import pysrc.body.logger
import pysrc.body.option_check
import pysrc.body.utilities
import pysrc.body.worker

__doc__ = ''' wrapper of SALMON , contains two phase: 1. index 2. quantify,
 quantify has two sub mode: align mode and reads mode .
'''
__author__ = 'zerodel'

SECTION_INDEX = "SALMON_INDEX"
SECTION_QUANTIFY = "SALMON_QUANTIFY"

_FILES_IN_INDEX_FOLDER = ['hash.bin',
                          'header.json',
                          'indexing.log',
                          'quasi_index.log',
                          # 'refInfo.json',
                          'rsd.bin',
                          'sa.bin',
                          'txpInfo.bin',
                          'versionInfo.json'
                          ]

_GENE_MAPPING_FILE = "--geneMap"

_QUANTIFY_OUTPUT = "--output"

_LIB_TYPE = "--libType"

_TARGET_TRANSCRIPTS = "--targets"

_ALIGNMENTS = "--alignments"

_GEN_CODE_FORMAT_FLAG = "--gencode"

_PERFECT_HASH = "--perfectHash"

_INDEX_TYPE = "--type"

_THREADS = "--threads"

_INDEX = "--index"

_TRANSCRIPTS = "--transcripts"

_SALMON_BIN = "salmon_bin"

_OPT_UN_MATED_READS = "--unmatedReads"

_OPT_MATE2 = "--mate2"

_OPT_MATE1 = "--mate1"

_OPT_KMER_LEN = "--kmerLen"

_logger = pysrc.body.logger.default_logger("SALMON")


def get_index_path(para_config=None, *args, **kwargs):
    updated_para = pysrc.body.cli_opts.merge_parameters(
        kwargs, para_config, "SALMON_INDEX")
    return updated_para["--out"]


def interpret_seq_files(input_files):
    if input_files:
        if isinstance(input_files, list):
            return {_OPT_MATE1: input_files[0],
                    _OPT_MATE2: input_files[-1]}
        elif isinstance(input_files, str):
            return {_OPT_UN_MATED_READS: input_files}
        else:
            raise TypeError(
                "Error@SALMON_interpret_seq_file: only list and string are allowed in assign SALMON input_files")
    else:
        raise FileNotFoundError(
            "Error@SALMON_interpret_seq_file: seems that no input file is assigned")


def is_path_contain_index(possible_index_path):
    for part in os.listdir(possible_index_path):
        if part not in _FILES_IN_INDEX_FOLDER:
            return False
    else:
        return True


def _get_cmd_make_index(para_dict):
    # warning : here lies a hard-coded name
    str_salmon_index_command = "{salmon_bin} index".format(
        salmon_bin=para_dict.pop(_SALMON_BIN))
    str_salmon_index_command += " " + \
        pysrc.body.cli_opts.enum_all_opts(para_dict)
    return str_salmon_index_command


def _is_legal_suffix_array_interval(x):
    s = x.strip()
    is_num = s.isdecimal()

    if not is_num:
        return False

    n = int(s)
    root = math.log2(n)
    int_root = round(root)
    tor = 0.00001
    if tor < math.fabs(root - int_root):
        return False
    else:
        if 2 ** int_root == n:
            return True


def get_cmd_quantify(opts):
    para_dict = copy.copy(opts)

    priority_order_alignment_mode = [_ALIGNMENTS, _TARGET_TRANSCRIPTS,
                                     _LIB_TYPE, _QUANTIFY_OUTPUT, _GENE_MAPPING_FILE]

    priority_order_reads_mode = [_INDEX, _OPT_MATE1, _OPT_MATE2, _OPT_UN_MATED_READS,
                                 _LIB_TYPE, _QUANTIFY_OUTPUT, _GENE_MAPPING_FILE]

    if "-a" in para_dict or _ALIGNMENTS in para_dict:
        return _build_salmon_cmd_with_order_given(para_dict, priority_order_alignment_mode)
    else:
        return _build_salmon_cmd_with_order_given(para_dict, priority_order_reads_mode)


def _build_salmon_cmd_with_order_given(para_dict, priority_order, phase="quant"):
    # -a in para_dict , input is bam file

    cmd_string = "{salmon_bin} {phase} ".format(salmon_bin=para_dict.pop(_SALMON_BIN),
                                                phase=phase)

    cmd_following_order = " ".join([pysrc.body.cli_opts.drop_key(key, para_dict)
                                    for key in priority_order if key in para_dict])

    cmd_remaining_part = pysrc.body.cli_opts.enum_all_opts(para_dict)

    cmd_string = cmd_string + cmd_following_order + " " + cmd_remaining_part

    return cmd_string

# todo need a accurate way to check salmon index folder .


def _check_valid_index(path):
    #    if os.path.exists(path) and os.path.isdir(path):
    #        files = os.listdir(path)
    #        return all([x in files for x in _FILES_IN_INDEX_FOLDER])
    return os.path.exists(path) and os.path.isdir(path)


def _is_legal_lib_type(str_lib_type):
    return str_lib_type in ["I", "IU", "U"]


def _check_index_options(para_dict=None):
    opt_checker = pysrc.body.option_check.OptionChecker(
        para_dict, name=SECTION_INDEX)
    opt_checker.one_and_only_one([_TRANSCRIPTS, "-t"], os.path.exists,
                                 FileNotFoundError(
                                     "Error@SALMON_INDEX: unable to find the transcript fa"),
                                 "Transcript fasta file(s)")

    opt_checker.one_and_only_one(["-i", _INDEX], pysrc.body.cli_opts.is_suitable_path_with_prefix,
                                 FileNotFoundError(
                                     "Error@SALMON_INDEX: wrong path for index output"),
                                 "Output stem [all files needed by SALMON will be of the form stem.*].")

    opt_checker.at_most_one(["-d", "--decoys"], lambda x: all([os.path.exists(y) for y in x.split()]),
                            FileNotFoundError("ERROR@SALMON_INDEX: unable to find decoy sequences"), "path to decoy sequences")

    opt_checker.at_most_one([_THREADS, "-p"], lambda x: x.isdecimal() and 0 < int(x) <= 2,
                            ValueError(
                                "ERROR@SALMON_INDEX: given incorrect threads number"),
                            "The number of threads to use concurrently. default is the number of your cpu cores")

    opt_checker.may_need(_INDEX_TYPE, lambda x: x.strip() == "quasi",
                         ValueError(
                             "ERROR@SALMON_INDEX: type of index is quasi , fmd will be removed"),
                         "The type of index to build")

    opt_checker.at_most_one(["-s", "--sasamp"], _is_legal_suffix_array_interval,
                            ValueError(
                                "ERROR@SALMON_INDEX: the interval of suffix array sampled should be a power of 2"),
                            """The interval at which the suffix array should be
                         sampled. Smaller values are faster, but produce a
                         larger index. The default should be OK, unless
                         your transcriptome is huge. This value should be a
                         power of 2.""")

    opt_checker.may_need(_PERFECT_HASH, lambda x: not x,
                         ValueError(
                             "ERROR@SALMON_INDEX: perfect-hash is a  flag only in quasi mode"),
                         """[quasi index only] Build the index using a perfect
                         hash rather than a dense hash.  This will require
                         less memory (especially during quantification),
                         but will take longer to construct""")

    opt_checker.may_need(_GEN_CODE_FORMAT_FLAG, lambda x: not x,
                         ValueError(
                             "ERROR@SALMON_INDEX: --gencode is only a flag , no value needed . "),
                         """This flag will expect the input transcript fasta
                         to be in GENCODE format, and will split the
                         transcript name at the first '|' character.  These
                         reduced names will be used in the output and when
                         looking for these transcripts in a gene to
                         transcript GTF.""")
    opt_checker.may_need(_OPT_KMER_LEN, lambda x: int(x) > 0,
                         ValueError(
                             "ERROR@SALMON_INDEX: invalid k-mer size given.... "),
                         """This int confirms the length of k-mer,this would be"""),
    opt_checker.forbid_these_args("-h", "--help", "-v", "--version")

    return opt_checker


def _check_quantify_options_alignment_mode(para_dict=None):
    opt_checker = pysrc.body.option_check.OptionChecker(
        para_dict, name=SECTION_QUANTIFY + "_ALIGNMENT_MODE")
    opt_checker.one_and_only_one(["-a", _ALIGNMENTS], os.path.exists,
                                 FileNotFoundError(
                                     "ERROR@SALMON_QUANTIFY: incorrect alignment file given "),
                                 "input alignment (BAM) file(s).")

    opt_checker.one_and_only_one(["-t", _TARGET_TRANSCRIPTS], os.path.exists,
                                 FileNotFoundError(
                                     "ERROR@SALMON_QUANTIFY: incorrect FASTA format file given"),
                                 " FASTA format file containing target transcripts.")

    opt_checker = _add_common_options_both_salmon_quantify_modes(
        opt_checker, para_dict)

    opt_checker.forbid_these_args("-i", _INDEX)

    return opt_checker


def _add_common_options_both_salmon_quantify_modes(opt_checker, para_dict):
    opt_checker.one_and_only_one(["-l", _LIB_TYPE], _is_legal_lib_type,
                                 KeyError(
                                     "Error@SALMON_QUANTIFY: no suitable libtype for SALMON"),
                                 "Format string describing the library type")

    opt_checker.one_and_only_one(["-o", _QUANTIFY_OUTPUT], os.path.exists,
                                 FileNotFoundError(
                                     "ERROR@SALMON_QUANTIFY: incorrect path to SALMON output "),
                                 "Output quantification file.")

    opt_checker.at_most_one([_GENE_MAPPING_FILE, "-g"], os.path.exists,
                            FileNotFoundError(
                                "ERROR@SALMON_QUANTIFY: incorrect gene map file (gtf ,csv) given"),
                            """File containing a mapping of transcripts to genes.""")

    opt_checker.forbid_these_args("-h", "--help", "-v", "--version")
    return opt_checker


def _check_quantify_options_reads_mode(para_dict=None):
    opt_checker = pysrc.body.option_check.OptionChecker(
        para_dict, name=SECTION_QUANTIFY + "_READS_MODE")
    opt_checker = _add_common_options_both_salmon_quantify_modes(
        opt_checker, para_dict)

    opt_checker.one_and_only_one(["-i", _INDEX], _check_valid_index,
                                 FileNotFoundError(
                                     "ERROR@SALMON_QUANTIFY_reads_mode: incorrect index path given "),
                                 "Salmon index")

    opt_checker.may_need(_OPT_UN_MATED_READS, os.path.exists,
                         FileNotFoundError(
                             "ERROR@SALMON_QUANTIFY_reads_mode: incorrect single-end file"),
                         "un-mated, ie, single-end sequence reads file")

    opt_checker.may_need(_OPT_MATE1, os.path.exists,
                         FileNotFoundError(
                             "ERROR@SALMON_QUANTIFY_reads_mode: incorrect paired-end mate1 file"),
                         "paired-end sequence reads file: mate1")
    opt_checker.may_need(_OPT_MATE2, os.path.exists,
                         FileNotFoundError(
                             "ERROR@SALMON_QUANTIFY_reads_mode: incorrect paired-end mate2 file"),
                         "paired-end sequence reads file: mate2")
    opt_checker.may_need("--validateMappings", lambda x: True,
                         ValueError(
                             "ERROR@SALMON_QUANTIFY_reads_mode: validation mapping"),
                         "flag to increase accuracy. became a default after salmon 1.0 version")

    if para_dict:
        if _OPT_MATE1 in para_dict and _OPT_MATE2 in para_dict:
            if _OPT_UN_MATED_READS in para_dict:
                raise KeyError(
                    "ERROR@SALMON_QUANTIFY_reads_mode: can not using paired-end and single-end at the same time")

    opt_checker.forbid_these_args("-a", _ALIGNMENTS, "-t", _TARGET_TRANSCRIPTS)

    return opt_checker


def _check_quantify_options(para_dict):
    if "-a" in para_dict or _ALIGNMENTS in para_dict:
        _check_quantify_options_alignment_mode(para_dict)
    else:
        _check_quantify_options_reads_mode(para_dict)


#

# %% start of main parts
opt_checker_index = _check_index_options()
opt_checker_quantify_align = _check_quantify_options_alignment_mode()
opt_checker_quantify_reads = _check_quantify_options_reads_mode()

OPTION_CHECKERS = [opt_checker_index,
                   opt_checker_quantify_reads, opt_checker_quantify_align]


def index(para_config=None, *args, **kwargs):
    opts_raw = pysrc.body.cli_opts.merge_parameters(
        kwargs, para_config, SECTION_INDEX)
    opts = copy.copy(opts_raw)

    _check_index_options(opts)

    cmd = _get_cmd_make_index(opts)

    pysrc.body.worker.run(cmd)

    return opts_raw


def quantify(para_config=None, *args, **kwargs):
    opts_raw = pysrc.body.cli_opts.merge_parameters(
        kwargs, para_config, SECTION_QUANTIFY)
    opts = copy.copy(opts_raw)

    para_dict = copy.copy(opts_raw)
    if "-a" in para_dict or _ALIGNMENTS in para_dict:
        opt_checker_quantify_align.check(para_dict)
    else:
        opt_checker_quantify_reads.check(para_dict)

    cmd = get_cmd_quantify(opts)

    pysrc.body.worker.run(cmd)

    return opts_raw


if __name__ == "__main__":
    print(__doc__)
    print(opt_checker_index)
    print(opt_checker_quantify_align)
    print(opt_checker_quantify_reads)
