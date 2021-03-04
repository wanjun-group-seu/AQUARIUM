# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import copy
import os
import os.path

import pysrc.body.cli_opts
import pysrc.body.logger
import pysrc.body.option_check
import pysrc.body.utilities
import pysrc.body.worker

_OPT_IS_PAIRED_END = "--paired-end"

_DESC_SAMPLE_NAME = """The name of the sample analyzed.
All output files are prefixed by this name (e.g., sample_name.genes.results)"""

_DESC_EXPRESSION_REFERENCE_NAME = """The name of the reference used.
The user must have run 'rsem-prepare-reference' with this reference_name before running this program."""

_DESC_REFERENCE_NAME = """The name of the reference used.
RSEM will generate several reference-related files that are prefixed by this name.
This name can contain path information (e.g. '/ref/mm9').
"""

_DESC_REFERENCE_FASTA_FILES = """Either a comma-separated list of Multi-FASTA formatted files OR a directory name.
If a directory name is specified, RSEM will read all files with suffix ".fa" or ".fasta" in this directory.
The files should contain either the sequences of transcripts or an entire genome, depending on whether the '--gtf' option is used.
"""

__doc__ = ''' wrapper of RSEM cmd line options, contains two phase : 1 index,2 quantify
'''
__author__ = 'zerodel'

SECTION_INDEX = "RSEM_INDEX"
SECTION_QUANTIFY = "RSEM_QUANTIFY"

_ESSENTIAL_ARGS_REFERENCE = "{reference_fasta_files} {reference_name}"
_ESSENTIAL_ARGS_QUANTIFICATION = "{reference_name} {sample_name}"
_RSEM_REFERENCE_SUFFIX = [".ti", ".grp", ".seq"]


_OPT_BAM = "--bam"

_FILE_SUFFIX_BAM = "Aligned.toTranscriptome.out.bam"


def confirm_rsem_style_name_folder(rsem_style_folder_prefix):
    folder, prefix = os.path.split(rsem_style_folder_prefix)
    if not os.path.isdir(folder):
        raise FileNotFoundError("ERROR@RSEM: RSEM style reference/output folder not exist")


def is_path_contain_index(prefix_reference):
    files = [prefix_reference + suffix for suffix in _RSEM_REFERENCE_SUFFIX]
    return all([os.path.exists(refer_file) for refer_file in files])


def get_index_path(updated_para):
    return updated_para["reference_name"]


def interpret_index_path(index_path):
    if index_path:
        return {"reference_name": index_path}
    else:
        return {}


def interpret_seq_files(input_files):
    if isinstance(input_files, list):
        return {"upstream_read_file": input_files[0],
                "downstream_read_file": input_files[-1]}
    elif isinstance(input_files, str):
        if input_files.endswith(".bam"):
            return {"input": input_files}
        else:
            return {"upstream_read_file": input_files}
    else:
        return {}


def _get_cmd_make_index(updated_para):
    cmd_prepare_reference = "{rsem_bin_prepare_reference}".format(**updated_para)

    pysrc.body.utilities.check_binary_executable(cmd_prepare_reference)

    cmd_prepare_reference += pysrc.body.cli_opts.cat_options(updated_para, updated_para)
    cmd_prepare_reference += _ESSENTIAL_ARGS_REFERENCE.format(**updated_para)

    return cmd_prepare_reference


def _get_cmd_calculate_expression(options):
    cmd_string_quantify = "{}".format(options.pop("rsem_bin_calculate_expression"))
    pysrc.body.utilities.check_binary_executable(cmd_string_quantify.strip())

    cmd_suffix_end = "{reference_name} {sample_name}".format(reference_name=options.pop("reference_name"),
                                                             sample_name=options.pop("sample_name"))

    # decide input mode .
    if "upstream_read_file" in options and options["upstream_read_file"]:
        if "downstream_read_file" in options and options["downstream_read_file"]:
            cmd_specify_input = "--paired-end {upstream_read_file} {downstream_read_file} ".format(
                upstream_read_file=options.pop("upstream_read_file"),
                downstream_read_file=options.pop("downstream_read_file"))

        else:
            cmd_specify_input = "{upstream_read_file} ".format(
                upstream_read_file=options.pop("upstream_read_file"))

    elif "alignment" in options and options["alignment"]:
        mapping_result_path = options.pop("alignment")
        is_paired_end = _OPT_IS_PAIRED_END if pysrc.body.utilities.is_paired(mapping_result_path) else ""
        cmd_specify_input = " ".join(["--alignments ", is_paired_end, mapping_result_path, cmd_suffix_end])

    elif _OPT_BAM in options and options[_OPT_BAM]:
        bam_format_alignment = options.pop(_OPT_BAM)
        is_paired_end = _OPT_IS_PAIRED_END  # todo: change this dumb implement into a real detection result
        cmd_specify_input = " ".join([_OPT_BAM, bam_format_alignment, is_paired_end])

    elif "seqs" in options and options["seqs"]:
        seqs = [x for x in options.pop("seqs").strip().split() if x]
        is_paired_seqs = _OPT_IS_PAIRED_END if len(seqs) == 2 else ""
        cmd_specify_input = is_paired_seqs + " " + " ".join(seqs)
    else:
        raise KeyError("Error : no argument specifying input file for rsem quantification")

    cmd_string_quantify = " ".join([cmd_string_quantify,
                                    pysrc.body.cli_opts.all_options(options),
                                    cmd_specify_input,
                                    cmd_suffix_end])

    return cmd_string_quantify


def _option_check_index_rsem(updated_para=None):
    opt_checker = pysrc.body.option_check.OptionChecker(updated_para, name=SECTION_INDEX)
    opt_checker.must_have("reference_fasta_files", os.path.exists,
                          FileNotFoundError("Error: fasta reference file not found as rsem reference"),
                          _DESC_REFERENCE_FASTA_FILES)

    opt_checker.must_have("reference_name", pysrc.body.cli_opts.is_suitable_path_with_prefix,
                          FileNotFoundError("Error: reference_name not specified"),
                          _DESC_REFERENCE_NAME)
    return opt_checker


def _option_check_quantify_rsem(updated_para=None):
    quantify_opt_check = pysrc.body.option_check.OptionChecker(updated_para, name=SECTION_QUANTIFY)
    quantify_opt_check.must_have("reference_name", is_path_contain_index,
                                 FileNotFoundError("ERROR: rsem reference file not specified"),
                                 _DESC_EXPRESSION_REFERENCE_NAME)
    quantify_opt_check.must_have("sample_name", pysrc.body.cli_opts.is_suitable_path_with_prefix,
                                 FileNotFoundError("ERROR: rsem sample file not specified"), _DESC_SAMPLE_NAME)
    return quantify_opt_check


_logger = pysrc.body.logger.default_logger("RSEM")
opt_checker_index = _option_check_index_rsem()
opt_checker_quantify = _option_check_quantify_rsem()

OPTION_CHECKERS = [opt_checker_index, opt_checker_quantify]


def index(para_config=None, *args, **kwargs):
    opts_index_rsem_raw = pysrc.body.cli_opts.merge_parameters(kwargs, para_config, SECTION_INDEX)

    opts_index_rsem = copy.copy(opts_index_rsem_raw)

    opt_checker_index.check(opts_index_rsem)

    cmd = _get_cmd_make_index(opts_index_rsem)

    if not is_path_contain_index(get_index_path(opts_index_rsem)):
        pysrc.body.worker.run(cmd)
    else:
        print("Report: already have a RSEM index in {}".format(get_index_path(opts_index_rsem)))

    return opts_index_rsem_raw


def quantify(para_config=None, *args, **kwargs):
    opts_quantify_rsem = pysrc.body.cli_opts.merge_parameters(kwargs, para_config, SECTION_QUANTIFY)

    opt_checker_quantify.check(copy.copy(opts_quantify_rsem))

    cmd = _get_cmd_calculate_expression(copy.copy(opts_quantify_rsem))

    pysrc.body.worker.run(cmd)

    return opts_quantify_rsem


if __name__ == "__main__":
    print(__doc__)
    print(opt_checker_index)
    print(opt_checker_quantify)
