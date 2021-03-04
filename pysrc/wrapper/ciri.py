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
import pysrc.file_format.ciri_entry
import pysrc.file_format.fa
import pysrc.file_format.fq
import pysrc.file_format.sam
import pysrc.wrapper.bwa



BWA_T_LONGER_60BP = "19"
BWA_T_SHORT_READS = "15"

__doc__ = ''' this is the wrapper of CIRI version 1, it contains one phase: detection
'''

__author__ = 'zerodel'

_OPT_BED_OUTPUT = "bed_out"
_OPT_FILTER_READS_NUMBER = "filter_junction_reads"

_OPT_ANNOTATION = "--anno"

_OPT_INPUT = "--in"

_OPT_OUTPUT = "--out"

_OPT_REF_DIR = "--ref-dir"
_OPT_REF_FILE = "--ref-file"
_OPT_SHOW_ALL = "--no-strigency"


_OPT_REF_DIR_CIRI2 = "--ref_dir"
_OPT_REF_FILE_CIRI2 = "--ref_file"
_OPT_SHOW_ALL_CIRI2 = "--no_strigency"

_ESSENTIAL_ARGUMENTS = [_OPT_INPUT, _OPT_OUTPUT, _OPT_REF_FILE, _OPT_ANNOTATION]


_ESSENTIAL_ARGUMENTS_CIRI2 = [_OPT_INPUT, _OPT_OUTPUT, _OPT_REF_FILE_CIRI2, _OPT_ANNOTATION]


# def _is_this_path_contains_valid_folder(path):
#     folder, out_file = os.path.split(path)
#     if out_file:
#         return os.path.exists(folder) and os.path.isdir(folder)
#     else:
#         raise FileNotFoundError("Error: CIRI output should be a file")

#
# def _export_mapping_of_circular_isoform(some_ciri_entry):
#     if some_ciri_entry.id and some_ciri_entry.gene_id:
#         return "\t".join([some_ciri_entry.id, some_ciri_entry.gene_id])
#     else:
#         return ""

# def get_alignment(opts):
#     if _OPT_INPUT in opts:
#         return opts[_OPT_INPUT]
#     else:
#         raise KeyError("Error: no input file specified for CIRI")


def check_ref(ref_path):
    if not os.path.exists(ref_path):
        raise FileNotFoundError("Error: Path not exist: {}".format(ref_path))
    else:
        if os.path.isdir(ref_path):
            return any([pysrc.file_format.fa.is_fasta(part) for part in os.listdir(ref_path)])

        elif os.path.isfile(ref_path):
            return pysrc.file_format.fa.is_fasta(ref_path)
        else:
            raise FileNotFoundError("Error: Path is not file or folder: {} ".format(ref_path))


# def _is_there_alignment_already(opts):
#     is_there_alignments = _OPT_INPUT in opts and os.path.exists(opts[_OPT_INPUT])
#     return is_there_alignments


# # end of helper functions . ###################################


is_general_aligner_needed = False

SECTION_DETECT = "CIRI"

"""ciri heavily depends on BWA aligner....but RSEM prefers bowtie/STAR,
and here we adopt the long-format arguments
"""

_logger = pysrc.body.logger.default_logger(SECTION_DETECT)


def _check_opts(args_dict=None):
    oc = pysrc.body.option_check.OptionChecker(args_dict, name=SECTION_DETECT)
    oc.may_need("bwa_bin", pysrc.body.utilities.which,
                FileNotFoundError("ERROR@CIRI: incorrect bwa binary path for bwa"),
                "binary file path of BWA aligner")

    oc.may_need("bwa_index", os.path.exists,
                FileNotFoundError("ERROR@CIRI: incorrect bwa index path"),
                "index path for BWA aligner")

    oc.must_have("ciri_path", os.path.exists,
                 FileNotFoundError("Error@CIRI: no ciri script "),
                 "file path to ciri script")

    oc.may_need("--seqs", pysrc.body.cli_opts.check_if_these_files_exist,
                FileNotFoundError("ERROR@CIRI: incorrect reads files provided "),
                "sequence reads files need analysis")

    oc.may_need(_OPT_ANNOTATION, os.path.exists,
                FileNotFoundError("ERROR@CIRI: incorrect annotation file "),
                "genomic annotation file")

    oc.must_have(_OPT_INPUT, pysrc.file_format.sam.is_valid_sam,
                 FileNotFoundError("Error : unable to find CIRI input file"),
                 "path to alignments in SAM file type")

    oc.must_have(_OPT_OUTPUT, pysrc.body.utilities.is_path_creatable,
                 FileNotFoundError("Error : incorrect CIRI output file"),
                 "CIRI detection report file")

    oc.one_and_only_one([_OPT_REF_DIR, _OPT_REF_FILE, _OPT_REF_FILE_CIRI2, _OPT_REF_DIR_CIRI2], check_ref,
                        FileNotFoundError("Error: incorrect CIRI ref file"),
                        "reference file for CIRI/CIRI2")

    oc.may_need(_OPT_BED_OUTPUT, pysrc.body.utilities.is_path_creatable,
                FileNotFoundError("Error: unable to create a bed file"),
                des_str="""summarized report in .bed format. indicating the region of putative circRNA"""
                )

    def is_suitable_reads_limit(x):
        try:
            ix = int(x)
        except:
            return False
        else:
            return ix >= 0

    oc.may_need(_OPT_FILTER_READS_NUMBER, is_suitable_reads_limit,
                ValueError("Error: lower limit of junction reads should be a non-negative integer"),
                des_str="""lower limit of junction reads detected by CIRI , 
                only bsj have more reads can be taken as a true positive """
                )

    oc.forbid_these_args("--help", "-H")

    def _ciri_input_check(opts):
        no_sam_for_ciri = not pysrc.file_format.sam.is_sam_from_bwa(opts[_OPT_INPUT])
        if no_sam_for_ciri:
            if "bwa_bin" not in opts or not pysrc.body.utilities.which(opts["bwa_bin"]):
                raise KeyError("ERROR@CIRI@NO_SAM: incorrect BWA binary file provided")

            if "bwa_index" not in opts or not os.path.exists(opts["bwa_index"]):
                raise KeyError("ERROR@CIRI@NO_SAM: incorrect BWA index file provided")

            if "--seqs" not in opts or not pysrc.body.cli_opts.check_if_these_files_exist(opts["--seqs"]):
                raise KeyError("ERROR@CIRI@NO_SAM: incorrect sequence reads for BWA alignment ")

    oc.custom_condition(_ciri_input_check, "check the options when no sam file for CIRI")

    return oc


opt_checker = _check_opts()

OPTION_CHECKERS = [opt_checker]


def _get_detect_cmd(opts_raw, pe):
    opts = copy.copy(opts_raw)
    cmd_corp = "perl {ciri_path}".format(ciri_path=opts.pop("ciri_path"))
    cmd_main = " ".join([pysrc.body.cli_opts.drop_key(key, opts) for key in _ESSENTIAL_ARGUMENTS])
    cmd_latter = pysrc.body.cli_opts.enum_all_opts(opts)
    cmd_pe_se = "--PE" if pe else "--SE"
    cmd_raw = " ".join([x for x in [cmd_corp, cmd_main, cmd_latter, cmd_pe_se] if x])
    return cmd_raw


def detect(par_dict=None, **kwargs):
    opts_of_index_phase_raw = pysrc.body.cli_opts.merge_parameters(kwargs, par_dict, SECTION_DETECT)

    opts = copy.copy(opts_of_index_phase_raw)

    bwa_bin_path = pysrc.body.cli_opts.extract_one(opts, key="bwa_bin")
    bwa_index_path = pysrc.body.cli_opts.extract_one(opts, key="bwa_index")
    reads = pysrc.body.cli_opts.extract_one(opts, "--seqs")
    is_pair_end = False if len(pysrc.body.cli_opts.extract_entries_from_value_str(reads)) < 2 else True

    bwa_score_cutoff = pysrc.body.cli_opts.extract_one(opts, "bwa_score", None)

    if _OPT_INPUT in opts and os.path.exists(opts[_OPT_INPUT]):
        _logger.debug("already have alignment file at {path_to_align_file}".format(path_to_align_file=opts.get(
            _OPT_INPUT, "")))
    else:
        if reads:
            try:
                map_job_setting = {
                    "bwa_bin": bwa_bin_path,
                    "bwa_index": bwa_index_path,
                    "read_file": reads,
                    "sam": opts[_OPT_INPUT],
                }
                if bwa_score_cutoff:
                    map_job_setting["bwa_score"] = bwa_score_cutoff

                _logger.info("bwa mapping with parameter: %s" % str(map_job_setting))

                bwa_mapping_ciri_only(map_job_setting)  # perform the mapping using bwa
            except Exception as e:
                _logger.error(" ERROR occurs during preparing the SAM file for CIRI")
                raise e  # no idea how to handle it

        else:
            _logger.error("no SAM file and no Reads")
            raise FileNotFoundError("Error@CIRI: no reads for BWA mapping")

    opt_checker.check(copy.copy(opts_of_index_phase_raw))

    ciri_report = opts.get(_OPT_OUTPUT)

    path_bed_out = pysrc.body.cli_opts.extract_one(opts, key=_OPT_BED_OUTPUT,
                                                   default=".".join([os.path.splitext(opts[_OPT_OUTPUT])[0], "bed"]))

    reads_lower_limit = pysrc.body.cli_opts.extract_one(opts,
                                                        _OPT_FILTER_READS_NUMBER,
                                                        default=pysrc.file_format.ciri_entry.JUNCTION_READS_LIMIT)

    # start invoking CIRI
    cmd_detect = _get_detect_cmd(opts, is_pair_end)
    _logger.info("raw command for CIRI is : %s" % cmd_detect)
    pysrc.body.worker.run(cmd_detect)  # perform the CIRI job here

    _logger.info("CIRI is done, now export the result into 'standard' .bed file as %s" % path_bed_out)
    pysrc.file_format.ciri_entry.transform_ciri_to_bed(ciri_output_file=ciri_report,
                                                       num_read_lower_limit=reads_lower_limit,
                                                       output_bed_file=path_bed_out)
    return opts_of_index_phase_raw


def bwa_mapping_ciri_only(meta_setting, **kwargs):
    try:
        bwa_bin = meta_setting["bwa_bin"] if "bwa_bin" in meta_setting else kwargs.pop("bwa_bin")
    except KeyError:
        raise KeyError("Error@CIRI: no binary file for BWA aligner")

    try:
        index_bwa = meta_setting["bwa_index"] if "bwa_index" in meta_setting else kwargs.pop("bwa_index")
    except KeyError:
        raise KeyError("Error@CIRI: no index file for BWA ")

    try:
        read_file = meta_setting["read_file"] if "read_file" in meta_setting else kwargs.pop("read_file")
    except KeyError:
        raise KeyError("Error@CIRI: no sequence reads file for BWA")

    try:
        sam = meta_setting["sam"] if "sam" in meta_setting else kwargs.pop("sam")
    except KeyError:
        raise KeyError("Error@CIRI: no sam output for BWA")

    try:
        core_used = meta_setting["thread_num"] if "thread_num" in meta_setting else kwargs.pop("thread_num")
    except KeyError:
        cores_num_total = int(pysrc.body.utilities.core_numbers_of_cpu())
        core_used = cores_num_total - 4 if cores_num_total > 5 else 1

    if not pysrc.wrapper.bwa.is_path_contain_index(index_bwa):
            pysrc.wrapper.bwa.index(para_config={"bwa_bin": bwa_bin, "in_fasta": index_bwa, "-a": "bwtsw"})

    fq_test = pysrc.body.cli_opts.extract_entries_from_value_str(read_file)[0]
    read_length = pysrc.file_format.fq.get_read_length(fq_test)

    # determine bwa score, 19 for normal data, 15 for short data.
    bwa_score_cutoff_given = meta_setting.get("bwa_score", None)
    bwa_score_cutoff_default = BWA_T_SHORT_READS if read_length < 60 else BWA_T_LONGER_60BP
    bwa_score_cutoff = bwa_score_cutoff_given if bwa_score_cutoff_given else bwa_score_cutoff_default

    pysrc.wrapper.bwa.align(para_config={
        "bwa_bin": bwa_bin,
        "read_file": read_file,
        "bwa_index": index_bwa,
        "-T": bwa_score_cutoff,
        "-t": str(core_used),  # todo: here we use all the cpu cores, greedy...
    }, output=sam)


def export_as_bed(par_dict=None, **kwargs):
    opts_of_index_phase_raw = pysrc.body.cli_opts.merge_parameters(kwargs, par_dict, SECTION_DETECT)

    opts = copy.copy(opts_of_index_phase_raw)

    ciri_report = opts[_OPT_OUTPUT]

    path_bed_out = opts.get(_OPT_BED_OUTPUT, os.path.join(os.path.splitext(opts[_OPT_OUTPUT])[0], ".bed"))

    reads_lower_limit = opts.get(_OPT_FILTER_READS_NUMBER, pysrc.file_format.ciri_entry.JUNCTION_READS_LIMIT)

    pysrc.file_format.ciri_entry.transform_ciri_to_bed(ciri_output_file=ciri_report,
                                                       num_read_lower_limit=reads_lower_limit,
                                                       output_bed_file=path_bed_out)


