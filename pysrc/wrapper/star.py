# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import copy
import os

import pysrc.body.cli_opts
import pysrc.body.logger
import pysrc.body.option_check
import pysrc.body.utilities
import pysrc.body.worker

__doc__ = ''' STAR commander line wrapper, contains two phase: 1.index, 2. align
'''
__author__ = 'zerodel'

_OPT_QUANT_MODE = "--quantMode"

_OPT_CHIM_SEGMENT_MIN = "--chimSegmentMin"

_OPT_OUT_FILE_NAME_PREFIX = "--outFileNamePrefix"

_OPT_SJDB_OVERHANG = "--sjdbOverhang"

_OPT_SJDB_GTF_FILE = "--sjdbGTFfile"

_OPT_GENOME_FASTA_FILES = "--genomeFastaFiles"

_OPT_RUN_MODE = "--runMode"

_OPT_THREAD_N = "--runThreadN"

_OPT_GENOME_DIR = "--genomeDir"

_OPT_READS_FILE = "--readFilesIn"

_OPT_STAR_BIN = "star_bin"

SECTION_INDEX = "STAR_INDEX"
SECTION_ALIGN = "STAR_ALIGN"

_run_mode_for_star = ["alignReads",
                      "genomeGenerate",
                      "inputAlignmentsFromBAM"]

_star_index_files = ['chrLength.txt',
                     'chrStart.txt',
                     'genomeParameters.txt',
                     'sjdbInfo.txt',
                     'transcriptInfo.tab',
                     'chrNameLength.txt',
                     'exonInfo.tab',
                     'SA',
                     'sjdbList.fromGTF.out.tab',
                     'chrName.txt',
                     'Genome',
                     'SAindex',
                     'sjdbList.out.tab']

_logger = pysrc.body.logger.default_logger("STAR")


def interpret_seq_files(input_files):
    if input_files:
        if isinstance(input_files, list):
            paths = " ".join(input_files)
        elif isinstance(input_files, str):
            paths = input_files
        else:
            raise TypeError("Error@interpret_seq_file_@STAR: only list and string are allowed in assign input_files")
    else:
        return {}

    return {_OPT_READS_FILE: paths}


def interpret_index_path(path_index):
    if not path_index:
        return {}
    if isinstance(path_index, list):
        paths = " ".join(path_index)
    elif isinstance(path_index, str):
        paths = path_index
    else:
        paths = ""
    return {_OPT_GENOME_DIR: paths}


def get_index_path(updated_para):
    return updated_para[_OPT_GENOME_DIR]


def is_path_contain_index(path):
    if os.path.isdir(path) and os.path.exists(path):
        dir_contents = os.listdir(path)
        if dir_contents:
            return all([os.path.exists(os.path.join(path, p)) for p in dir_contents])
        else:
            return False
    else:
        raise FileNotFoundError("Error: this is not a folder: {}".format(path))


def check_star_run_mode(x):
    return x in _run_mode_for_star


def _get_cmd_index(opt_index_phrase):
    cmd_star_str = "{}".format(opt_index_phrase.pop(_OPT_STAR_BIN))
    cmd_star_str += " " + pysrc.body.cli_opts.enum_all_opts(opt_index_phrase)
    return cmd_star_str.strip()


_DESC_OUT_FILE_NAME_PREFIX = """string: output files name prefix (including full or relative path).
Can only be defined on the command line"""

_DESC_READ_FILES_IN = "string(s): paths to files that contain input read1 (and, if needed, read2)"

_DESC_SJDB_OVER_HANG = """int>=0: length of the donor/acceptor sequence on each side of the junctions, ideally = (mate length - 1)
                                if =0, splice junction database is not used"""

_DESC_SJDB_GTF_FILE = "string: path to the GTF file with annotations"

_DESC_GENOME_FASTA_FILES = "string(s): path(s) to genomic fasta file for genome generation, separated by spaces."

_DESC_GENOME_DIR = "string: path to the directory where genome files are stored, or will be generated"

_DESC_STAR_BIN = "binary file path of STAR"

_DESC_RUN_THREAD = "number of cpu cores"

_DESC_RUN_MODE_STAR = """
alignReads: map reads
genomeGenerate: generate genome files
inputAlignmentsFromBAM: input alignments from BAM. Presently only works with –outWigType and –bamRemoveDuplicates."""


def _get_cmd_align(updated_para):
    cmd_star_align = "{}".format(updated_para.pop(_OPT_STAR_BIN))
    cmd_star_align += " " + pysrc.body.cli_opts.enum_all_opts(updated_para)
    return cmd_star_align.strip()


def is_map_result_already_exists(para_config=None, **kwargs):
    align_phrase_options = pysrc.body.cli_opts.merge_parameters(kwargs, para_config, SECTION_ALIGN)

    align_file_path = get_align_result_path(para_config, **kwargs)

    sub_dir, prefix = os.path.split(align_file_path)

    if not os.path.isdir(sub_dir) or not os.path.exists(sub_dir):
        return False

    if not os.path.exists(align_file_path + "SJ.out.tab"):
        return False

    # the STAR gene fusion is switched on .  should have chimeric.out.sam/ and chimeric.out.junction
    if _OPT_CHIM_SEGMENT_MIN in align_phrase_options:
        has_chim_junction = os.path.exists(align_file_path + "Chimeric.out.junction")
        has_chim_alignments = os.path.exists(align_file_path + "Chimeric.out.sam") or os.path.exists(
            align_file_path + "Chimeric.out.sorted.bam")

        if has_chim_alignments and has_chim_junction:
            pass
        else:
            return False

    if _OPT_QUANT_MODE in align_phrase_options and not os.path.exists(
                    align_file_path + "Aligned.toTranscriptome.out.bam"):
        return False

    return True


def get_align_result_path(para_config=None, **kwargs):
    align_phrase_options_raw = pysrc.body.cli_opts.merge_parameters(kwargs, para_config, SECTION_ALIGN)
    try:
        target = align_phrase_options_raw[_OPT_OUT_FILE_NAME_PREFIX]
    except KeyError:
        raise KeyError("Error@STAR: alignment output file name prefix not specified at option : '--outFileNamePrefix'")
    return target


def _option_check_index_phrase(opt_index_phrase=None):
    template_err = "Error@STAR_INDEX: %s"
    opt_checker = pysrc.body.option_check.OptionChecker(opt_index_phrase, name=SECTION_INDEX)

    opt_checker.must_have(_OPT_STAR_BIN, pysrc.body.utilities.which,
                          FileNotFoundError(template_err % "no such binary file "),
                          _DESC_STAR_BIN)

    opt_checker.must_have(_OPT_THREAD_N, pysrc.body.utilities.is_thread_num_less_than_core_number,
                          OSError(template_err % "too may threads"),
                          _DESC_RUN_THREAD)

    opt_checker.must_have(_OPT_RUN_MODE, lambda x: x.strip() == "genomeGenerate",
                          ValueError(template_err % "wrong run mode for STAR"),
                          _DESC_RUN_MODE_STAR)

    opt_checker.must_have(_OPT_GENOME_DIR, lambda x: os.path.exists(x) and os.path.isdir(x),
                          FileNotFoundError(template_err % " genomeDir not exist"),
                          _DESC_GENOME_DIR)

    opt_checker.must_have(_OPT_GENOME_FASTA_FILES, pysrc.body.cli_opts.check_if_these_files_exist,
                          FileNotFoundError(template_err % "Unable to find genome fasta file for STAR"),
                          _DESC_GENOME_FASTA_FILES)

    opt_checker.must_have(_OPT_SJDB_GTF_FILE, os.path.exists,
                          FileNotFoundError(template_err % "Unable to find genome Annotation file"),
                          _DESC_SJDB_GTF_FILE)

    opt_checker.must_have(_OPT_SJDB_OVERHANG, lambda x: int(x) >= 0,
                          ValueError(template_err % "STAR overhang should be an non-negative integer"),
                          _DESC_SJDB_OVER_HANG)

    return opt_checker


def _option_check_align_phrase(updated_para=None):
    err_temp = "Error@STAR_ALIGN: %s"
    opt_checker = pysrc.body.option_check.OptionChecker(updated_para, name=SECTION_ALIGN)
    opt_checker.must_have(_OPT_STAR_BIN, pysrc.body.utilities.which,
                          FileNotFoundError(err_temp % "no such STAR binary file"),
                          _DESC_STAR_BIN)

    opt_checker.must_have(_OPT_THREAD_N, pysrc.body.utilities.is_thread_num_less_than_core_number,
                          OSError(err_temp % "too may threads"),
                          _DESC_RUN_THREAD)

    opt_checker.must_have(_OPT_GENOME_DIR, is_path_contain_index,
                          FileExistsError(err_temp % "STAR index path 'genomeDir' not exist "),
                          _DESC_GENOME_DIR)

    opt_checker.must_have(_OPT_READS_FILE, pysrc.body.cli_opts.check_if_these_files_exist,
                          FileNotFoundError(err_temp % "Unable to find read files for STAR "),
                          _DESC_READ_FILES_IN)

    opt_checker.must_have(_OPT_OUT_FILE_NAME_PREFIX, pysrc.body.cli_opts.is_suitable_path_with_prefix,
                          FileNotFoundError(err_temp % "unable to set the alignment output path for STAR "),
                          _DESC_OUT_FILE_NAME_PREFIX)

    opt_checker.may_need(_OPT_QUANT_MODE, lambda x: x in ["-", "TranscriptomeSAM"],
                         KeyError(err_temp % "incorrect quant mode given "),
                         "specify how the output sam file should be ")

    return opt_checker


opt_checker_index = _option_check_index_phrase()
opt_checker_align = _option_check_align_phrase()
OPTION_CHECKERS = [opt_checker_index, opt_checker_align]


def index(para_config=None, *args, **kwargs):
    opts_of_index_phase_raw = pysrc.body.cli_opts.merge_parameters(kwargs, para_config, SECTION_INDEX)

    opts_of_index_phase = copy.copy(opts_of_index_phase_raw)
    opt_checker_index.check(copy.copy(opts_of_index_phase_raw))

    dir_index = get_index_path(opts_of_index_phase)

    if not is_path_contain_index(dir_index):
        cmd_index = _get_cmd_index(opts_of_index_phase)
        pysrc.body.worker.run(cmd_index)
    else:
        print("Report: already have a STAR index in {}".format(dir_index))

    return opts_of_index_phase_raw


def align(para_config=None, *args, **kwargs):
    align_phrase_options_raw = pysrc.body.cli_opts.merge_parameters(kwargs, para_config, SECTION_ALIGN)

    align_phrase_options = copy.copy(align_phrase_options_raw)

    opt_checker_align.check(copy.copy(align_phrase_options_raw))

    if not is_map_result_already_exists(para_config, **kwargs):
        cmd = _get_cmd_align(align_phrase_options)
        pysrc.body.worker.run(cmd)

    return align_phrase_options_raw


if __name__ == "__main__":
    print(__doc__)
    print(opt_checker_index)
    print(opt_checker_align)
