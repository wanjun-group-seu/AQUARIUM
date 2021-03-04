# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:

import copy
import os
import os.path

import pysrc.body.cli_opts
import pysrc.body.logger
import pysrc.body.option_check
import pysrc.body.utilities
import pysrc.body.worker

_OPT_KMER_SIZE = "--kmerSize"

__doc__ = ''' wrapper for Sailfish , has two phase: 1. index , 2. quantify
'''
__author__ = 'zerodel'

SECTION_INDEX = "SAILFISH_INDEX"
SECTION_QUANTIFY = "SAILFISH_QUANTIFY"

# sailfish version now is 0.9.0
_FILES_IN_INDEX_FOLDER = ['hash.bin',
                          'header.json',
                          'logs',
                          'rsd.bin',
                          'sa.bin',
                          'txpInfo.bin',
                          'versionInfo.json']

_GENE_MAPPING_FILE = "--geneMap"

_QUANTIFY_OUTPUT = "--output"

_LIB_TYPE = "--libType"

_THREADS = "--threads"

_INDEX = "--index"

_TRANSCRIPTS = "--transcripts"

_SAILFISH_BIN = "sailfish_bin"

_OPT_UN_MATED_READS = "--unmatedReads"

_OPT_MATE2 = "--mate2"

_OPT_MATE1 = "--mate1"

_logger = pysrc.body.logger.default_logger("SAILFISH")


def get_index_path(para_config=None, *args, **kwargs):
    updated_para = pysrc.body.cli_opts.merge_parameters(kwargs, para_config, "SAILFISH_INDEX")
    return updated_para["--out"]


def interpret_seq_files(input_files):
    if input_files:
        if isinstance(input_files, list):
            return {_OPT_MATE1: input_files[0],
                    _OPT_MATE2: input_files[-1]}
        elif isinstance(input_files, str):
            return {_OPT_UN_MATED_READS: input_files}
        else:
            raise TypeError("Error@sailfish: only list and string are allowed in assign sailfish input_files")
    else:
        _logger.warning("seems that no input file is assigned")
        return {}


def interpret_index_path(index_path):
    if index_path and os.path.exists(index_path):
        return {_INDEX: index_path}
    else:
        return {}


def is_path_contain_index(supposed_path_to_index):
    if os.path.isdir(supposed_path_to_index):
        files_under_index_dir = os.listdir(supposed_path_to_index)
        if files_under_index_dir:
            return all([os.path.exists(os.path.join(supposed_path_to_index, p)) for p in files_under_index_dir])

    return False  # false for all other cases


def get_cmd_make_index(para_dict):
    # warning: hard coded here
    cmd_sailfish_index = "{sailfish_bin} index".format(sailfish_bin=para_dict.pop(_SAILFISH_BIN))
    cmd_sailfish_index += " " + pysrc.body.cli_opts.enum_all_opts(para_dict)
    return cmd_sailfish_index


def _get_cmd_quantify(para_dict):
    priority_order = [_INDEX, _OPT_MATE1, _OPT_MATE2, _OPT_UN_MATED_READS,
                      _LIB_TYPE, _QUANTIFY_OUTPUT, _GENE_MAPPING_FILE]

    cmd_string = "{salmon_bin} {phase} ".format(salmon_bin=para_dict.pop(_SAILFISH_BIN),
                                                phase="quant")

    cmd_following_order = " ".join([pysrc.body.cli_opts.drop_key(key, para_dict)
                                    for key in priority_order if key in para_dict])

    cmd_remaining_part = pysrc.body.cli_opts.enum_all_opts(para_dict)
    cmd_string = cmd_string + cmd_following_order + " " + cmd_remaining_part
    return cmd_string


def _check_index_options(para_dict=None):
    opt_checker = pysrc.body.option_check.OptionChecker(para_dict, name=SECTION_INDEX)
    opt_checker.must_have("--transcripts", os.path.exists,
                          FileNotFoundError("Error@SAILFISH_INDEX: unable to find the transcript fa"),
                          "Transcript fasta file(s)")

    opt_checker.must_have("--out", pysrc.body.cli_opts.is_suitable_path_with_prefix,
                          FileNotFoundError("Error@SAILFISH_INDEX: wrong path for index output"),
                          "Output stem [all files needed by Sailfish will be of the form stem.*].")

    opt_checker.may_need(_OPT_KMER_SIZE, lambda x: x.isdecimal() and 0 < int(x) < 32,
                         ValueError("ERROR@SAILFISH_INDEX:  incorrect Kmer length"),
                         "Kmer size. default is 31")

    opt_checker.may_need(_THREADS, lambda x: x.isdecimal() and 0 < int(x) <= pysrc.body.utilities.core_numbers_of_cpu(),
                         ValueError("ERROR@SAILFISH_INDEX:  incorrect threads number "),
                         "The number of threads to use concurrently. default is the number of your cpu cores")

    opt_checker.forbid_these_args("-h", "--help", "-v", "--version")
    return opt_checker


def _check_quantify_options(para_dict=None):
    check_args = pysrc.body.option_check.OptionChecker(para_dict, name=SECTION_QUANTIFY)
    check_args.must_have(_INDEX, os.path.exists,
                         FileNotFoundError("Error@SAILFISH_QUANTIFY: unable to find sailfish index"),
                         "path to Sailfish index")

    check_args.must_have(_LIB_TYPE, _is_legal_lib_type,
                         KeyError("Error@SAILFISH_QUANTIFY: no suitable libtype for sailfish"),
                         "Format string describing the library type")

    check_args.must_have(_QUANTIFY_OUTPUT, os.path.exists,
                         FileNotFoundError("ERROR@SAILFISH_QUANTIFY: incorrect path to sailfish output "),
                         "Output quantification file.")

    check_args.may_need(_GENE_MAPPING_FILE, os.path.exists,
                        FileNotFoundError("ERROR@SAILFISH_QUANTIFY: incorrect gene map file (gtf ,csv) given"),
                        """File containing a mapping of transcripts to genes.""")

    check_args.may_need(_OPT_UN_MATED_READS, os.path.exists,
                        FileNotFoundError("ERROR@SAILFISH_QUANTIFY_reads_mode: incorrect single-end file"),
                        "un-mated, ie, single-end sequence reads file")

    check_args.may_need(_OPT_MATE1, os.path.exists,
                        FileNotFoundError("ERROR@SAILFISH_QUANTIFY_reads_mode: incorrect paired-end mate1 file"),
                        "paired-end sequence reads file: mate1")
    check_args.may_need(_OPT_MATE2, os.path.exists,
                        FileNotFoundError("ERROR@SAILFISH_QUANTIFY_reads_mode: incorrect paired-end mate2 file"),
                        "paired-end sequence reads file: mate2")
    if para_dict:
        if _OPT_MATE1 in para_dict and _OPT_MATE2 in para_dict:
            if _OPT_UN_MATED_READS in para_dict:
                raise KeyError(
                    "ERROR@SAILFISH_QUANTIFY_reads_mode: can not using paired-end and single-end at the same time")

    check_args.forbid_these_args("-h", "--help", "-v", "--version")
    return check_args


def _is_legal_lib_type(str_lib_type):
    return str_lib_type.strip('\"\'') in ["I", "IU", "U"]


opt_checker_index = _check_index_options()
opt_checker_quantify = _check_quantify_options()

OPTION_CHECKERS = [opt_checker_index, opt_checker_quantify]


def index(para_config=None, *args, **kwargs):
    opts_raw = pysrc.body.cli_opts.merge_parameters(kwargs, para_config, SECTION_INDEX)

    opts_raw = __infer_proper_threads_num(opts_raw)

    opts = copy.copy(opts_raw)

    opt_checker_index.check(copy.copy(opts_raw))

    cmd = get_cmd_make_index(opts)
    pysrc.body.worker.run(cmd)

    return opts_raw


def quantify(para_config=None, *args, **kwargs):
    opts_raw = pysrc.body.cli_opts.merge_parameters(kwargs, para_config, SECTION_QUANTIFY)

    opts_raw = __infer_proper_threads_num(opts_raw)  # if some thing is wrong , use sailfish default
    opts = copy.copy(opts_raw)
    opt_checker_quantify.check(copy.copy(opts_raw))

    cmd = _get_cmd_quantify(opts)
    pysrc.body.worker.run(cmd)

    return opts_raw


def __infer_proper_threads_num(dict_par):
    raw_x = dict_par.get(_THREADS, None)
    if raw_x:
        if raw_x.isdecimal() and 0 < int(raw_x) <= pysrc.body.utilities.core_numbers_of_cpu():
            _logger.debug("seems user defined thread number is legit : {num_threads}".format(num_threads=raw_x))
        else:
            _logger.warning("user defined thread number is illegal : {num}, and been removed, use sailfish "
                            "default".format(num=raw_x))
            dict_par.pop(_THREADS)
    else:
        _logger.debug("thread number not specified, use default of sailfish itself")
    return dict_par


if __name__ == "__main__":
    print(__doc__)
    print(opt_checker_index)
    print(opt_checker_quantify)
