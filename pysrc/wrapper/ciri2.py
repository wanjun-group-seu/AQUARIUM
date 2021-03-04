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
import pysrc.wrapper.ciri


BWA_T_LONGER_60BP = "19"
BWA_T_SHORT_READS = "15"

__doc__ = ''' this is the wrapper of CIRI version 2, it contains one phase: detection
but it is abandoned since 19-03-13
'''

__author__ = 'zerodel'

_OPT_BED_OUTPUT = "bed_out"
_OPT_FILTER_READS_NUMBER = "filter_junction_reads"

_OPT_ANNOTATION = "--anno"

_OPT_INPUT = "--in"

_OPT_OUTPUT = "--out"

_OPT_REF_DIR = "--ref_dir"
_OPT_REF_FILE = "--ref_file"
_OPT_SHOW_ALL = "--no_strigency"

# BWA_T_LONGER_60BP = BWA_T_LONGER_60BP
# BWA_T_SHORT_READS = BWA_T_SHORT_READS

_ESSENTIAL_ARGUMENTS = [_OPT_INPUT, _OPT_OUTPUT, _OPT_REF_FILE, _OPT_ANNOTATION]


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


is_general_aligner_needed = False

SECTION_DETECT = "CIRI"

_logger = pysrc.body.logger.default_logger(SECTION_DETECT)


def _check_opts(args_dict=None):
    oc = pysrc.body.option_check.OptionChecker(args_dict, name=SECTION_DETECT)
    oc.may_need("bwa_bin", pysrc.body.utilities.which,
                FileNotFoundError("ERROR@CIRI2: incorrect bwa binary path for bwa"),
                "binary file path of BWA aligner")

    oc.may_need("bwa_index", os.path.exists,
                FileNotFoundError("ERROR@CIRI2: incorrect bwa index path"),
                "index path for BWA aligner")

    oc.must_have("ciri_path", os.path.exists,
                 FileNotFoundError("Error@CIRI2: no ciri script "),
                 "file path to ciri script")

    oc.may_need("--seqs", pysrc.body.cli_opts.check_if_these_files_exist,
                FileNotFoundError("ERROR@CIRI2: incorrect reads files provided "),
                "sequence reads files need analysis")

    oc.may_need(_OPT_ANNOTATION, os.path.exists,
                FileNotFoundError("ERROR@CIRI2: incorrect annotation file "),
                "genomic annotation file")

    oc.at_most_one(["--thread_num", "-T"], lambda x: x.isdecimal(),
                   ValueError("Error@CIRI2: thread num should be a number"),
                   "thread number. CIRI2 only ")

    oc.must_have(_OPT_INPUT, pysrc.file_format.sam.is_valid_sam,
                 FileNotFoundError("Error : unable to find CIRI input file"),
                 "path to alignments in SAM file type")

    oc.must_have(_OPT_OUTPUT, pysrc.body.utilities.is_path_creatable,
                 FileNotFoundError("Error: incorrect CIRI output file"),
                 "CIRI detection report file")

    oc.one_and_only_one([_OPT_REF_DIR, _OPT_REF_FILE, ], check_ref,
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

    oc.may_need(_OPT_SHOW_ALL, lambda x: True,
                ValueError("ERROR@CIRI: whether to use no-stringency is a flag, No value needed"),
                des_str="""a flag to decide whether all the possible BSJ should be shown 
                in final CIRI detection report""")

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


def _get_detect_cmd(opts_raw):
    opts = copy.copy(opts_raw)
    cmd_corp = "perl {ciri_path}".format(ciri_path=opts.pop("ciri_path"))
    cmd_main = " ".join([pysrc.body.cli_opts.drop_key(key, opts) for key in _ESSENTIAL_ARGUMENTS])
    cmd_latter = pysrc.body.cli_opts.enum_all_opts(opts)
    cmd_raw = " ".join([x for x in [cmd_corp, cmd_main, cmd_latter] if x])
    return cmd_raw


def detect(par_dict=None, **kwargs):
    opts_of_index_phase_raw = pysrc.body.cli_opts.merge_parameters(kwargs, par_dict, SECTION_DETECT)

    opts = copy.copy(opts_of_index_phase_raw)

    bwa_bin_path = pysrc.body.cli_opts.extract_one(opts, key="bwa_bin")
    bwa_index_path = pysrc.body.cli_opts.extract_one(opts, key="bwa_index")
    reads = pysrc.body.cli_opts.extract_one(opts, "--seqs")

    thread_num = opts.get("--thread_num", "1")

    bwa_score_cutoff = pysrc.body.cli_opts.extract_one(opts, "bwa_score", None)

    if _OPT_INPUT in opts and os.path.exists(opts[_OPT_INPUT]):
        _logger.debug("already have alignment file at {path_to_align_file}".format(path_to_align_file=opts.get(
            _OPT_INPUT, "")))
    else:
        _logger.debug("no alignment file ,  mapping using bwa")

        if reads:
            try:
                map_job_setting = {
                    "bwa_bin": bwa_bin_path,
                    "bwa_index": bwa_index_path,
                    "read_file": reads,
                    "sam": opts[_OPT_INPUT],
                    "thread_num": thread_num,
                }
                if bwa_score_cutoff:
                    map_job_setting["bwa_score"] = bwa_score_cutoff

                _logger.info("bwa mapping with parameter: %s" % str(map_job_setting))

                pysrc.wrapper.bwa.bwa_mapping_ciri_only(map_job_setting)  # perform the mapping using bwa
                _logger.debug("mapping completed")

            except Exception as e:
                _logger.error(" ERROR occurs during preparing the SAM file for CIRI")
                raise e  # no idea how to handle it

        else:
            _logger.error("no SAM file and no Reads")
            raise FileNotFoundError("Error@CIRI: no reads for BWA mapping")

    opt_checker.check(copy.copy(opts_of_index_phase_raw))

    ciri_report = opts.get(_OPT_OUTPUT)

    path_bed_out = pysrc.body.cli_opts.extract_one(opts, key=_OPT_BED_OUTPUT,
                                                   default=os.path.splitext(opts[_OPT_OUTPUT])[0] + ".bed")

    reads_lower_limit = pysrc.body.cli_opts.extract_one(opts,
                                                        _OPT_FILTER_READS_NUMBER,
                                                        default=pysrc.file_format.ciri_entry.JUNCTION_READS_LIMIT)

    # start invoking CIRI
    cmd_detect = _get_detect_cmd(opts)
    _logger.info("raw command for CIRI2 is : %s" % cmd_detect)
    pysrc.body.worker.run(cmd_detect)  # perform the CIRI job here

    _logger.info("CIRI2 is done, now export the result into 'standard' .bed file")
    pysrc.file_format.ciri_entry.transform_ciri_to_bed(ciri_output_file=ciri_report,
                                                       num_read_lower_limit=reads_lower_limit,
                                                       output_bed_file=path_bed_out)
    return opts_of_index_phase_raw


def export_as_bed(par_dict=None, **kwargs):
    opts_of_index_phase_raw = pysrc.body.cli_opts.merge_parameters(kwargs, par_dict, SECTION_DETECT)

    opts = copy.copy(opts_of_index_phase_raw)

    ciri_report = opts[_OPT_OUTPUT]

    path_bed_out = opts.get(_OPT_BED_OUTPUT, os.path.join(os.path.splitext(opts[_OPT_OUTPUT])[0], ".bed"))

    reads_lower_limit = opts.get(_OPT_FILTER_READS_NUMBER, pysrc.file_format.ciri_entry.JUNCTION_READS_LIMIT)

    pysrc.file_format.ciri_entry.transform_ciri_to_bed(ciri_output_file=ciri_report,
                                                       num_read_lower_limit=reads_lower_limit,
                                                       output_bed_file=path_bed_out)
