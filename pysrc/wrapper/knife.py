# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import copy
import os
import os.path
import re

import pysrc.body.cli_opts
import pysrc.body.logger
import pysrc.body.option_check
import pysrc.body.utilities
import pysrc.body.worker
import pysrc.file_format.fq

_OPT_DECOY_RATIO = "filter_decoy_ratio"

_OPT_NAIVE_PV = "filter_naive_pv"

_OPT_GLM_PV = "filter_glm_pv"

PV_GLM_DEFAULT = 0.9

PV_NAIVE_DEFAULT = 0.9

DECOY_RATIO_DEFAULT = 0.1

_OPT_BASH_BIN = "bash_bin"

_OPT_KNIFE_SCRIPT = "knife_script"

_OPT_INDEX_PATH = "index_path"

_OPT_JUNCTION_OVERLAP = "junction_overlap"

_OPT_DATA_SET_NAME = "dataset_name"

_OPT_ALIGNMENT_PARENT_DIRECTORY = "alignment_parent_directory"

_OPT_READ_ID_STYLE = "read_id_style"

_OPT_READ_DIRECTORY = "read_directory"

_OPT_BED_OUTPUT = "bed_out"

_READ_ID_STYLE_SINGE_END_OR_IDENTICAL = "complete"

_READ_ID_STYLE_PAIRED_END_UN_EQUAL = "appended"

__doc__ = ''' this is a wrapper for KNIFE(Known and Novel IsoForm Explorer)
'''
__author__ = 'zerodel'

SECTION_DETECT = "KNIFE"

_logger = pysrc.body.logger.default_logger(SECTION_DETECT)


def _is_this_folder_existing(x):
    return os.path.exists(x) and os.path.isdir(x)


def _check_valid_read_directory(some_path):
    if os.path.exists(some_path) and os.path.isdir(some_path):
        def is_fq(x):
            f_name, f_suffix = os.path.splitext(x)
            if f_name and f_suffix and f_suffix in [".fastq", ".fq.gz", ".fq"]:
                return True
            else:
                return False

        files = os.listdir(some_path)
        fqs = [f for f in files if is_fq(f)]

        base_names = [os.path.splitext(f)[0] for f in fqs]
        samples = set([base_name[:-1] for base_name in base_names])
        if samples:
            at_least_a_pair = False
            for sample in samples:
                sample1 = "%s1" % sample
                sample2 = "%s2" % sample

                # todo: need a better check on fastq reads
                if sample1 in base_names:
                    return True

                    # if (sample1 in base_names) != (sample2 in base_names):
                    #     return False
                    # if (sample1 in base_names) and (sample2 in base_names):
                    #     at_least_a_pair = True
            else:
                return at_least_a_pair

    return False


def _get_detect_cmd(para_dict):
    opts = copy.copy(para_dict)

    priority_order = [_OPT_BASH_BIN, _OPT_KNIFE_SCRIPT, _OPT_READ_DIRECTORY, _OPT_READ_ID_STYLE,
                      _OPT_ALIGNMENT_PARENT_DIRECTORY,
                      _OPT_DATA_SET_NAME, _OPT_JUNCTION_OVERLAP, _OPT_INDEX_PATH]

    command_string_must_have = " ".join([opts.pop(key, "") for key in priority_order])
    command_latter_part = " ".join([opts[key] for key in opts])

    if command_latter_part:
        command_string_must_have = command_string_must_have + " " + command_latter_part

    return command_string_must_have


def _predict_id_style_and_overlap_length(opts):
    directory_fqs = opts[_OPT_READ_DIRECTORY]
    files_in_reads_folder = os.listdir(directory_fqs)
    mate1, mate2 = _get_mate_files(files_in_reads_folder)

    mate1 = os.path.join(directory_fqs, mate1) if mate1 else ""
    mate2 = os.path.join(directory_fqs, mate2) if mate2 else ""

    opts[_OPT_READ_ID_STYLE] = _predict_reads_id_type_for_fq_files(mate1, mate2)
    opts[_OPT_JUNCTION_OVERLAP] = _pick_suitable_overlap(mate1, mate2)

    return opts


def _pick_suitable_overlap(mate1, mate2):
    if not mate1:
        raise FileNotFoundError("Error@picking_suitable_overlap_parameter: seems not fastq file")
    is_single_end = not mate2
    read_length = pysrc.file_format.fq.get_read_length(mate1)
    threshold_length_read = 70
    if is_single_end:
        overlap_len = 10 if read_length < threshold_length_read else 15
    else:
        overlap_len = 8 if read_length < threshold_length_read else 13

    return str(overlap_len)


def _predict_reads_id_type_for_fq_files(mate1, mate2):
    if (not mate2) or pysrc.file_format.fq.is_pair_end_fastq_id_identical(mate1, mate2):
        return _READ_ID_STYLE_SINGE_END_OR_IDENTICAL
    else:
        return _READ_ID_STYLE_PAIRED_END_UN_EQUAL


def report_suitable_overlap(mate1, mate2):
    print(_pick_suitable_overlap(mate1, mate2))


def report_id_type(mate1, mate2):
    print(_predict_reads_id_type_for_fq_files(mate1, mate2))


def _get_mate_files(files_in_reads_folder):
    _PATTERN_SAMPLE = ".*(?=[-_]?[12]\.f.{0,5}q$)"
    sample_reo = re.compile(_PATTERN_SAMPLE)
    fqs_samples = [sample_reo.findall(x) for x in files_in_reads_folder if sample_reo.match(x)]
    samples_id = list(set([fq[0] for fq in fqs_samples if fq[0]]))
    # here we only support one-sample-one-folder
    if len(samples_id) < 1:
        raise FileNotFoundError("Error@predict_fq_read_id_type: seems you do not have fastq files")
    if len(samples_id) > 1:
        raise FileExistsError("Error@predict_fq_read_id_type: multiple sample in a single folder is not supported ")
    ss_id = samples_id[0]
    pattern_mate1_this_sample = "%s[-_]?1\.f.{1,5}$" % ss_id
    pattern_mate2_this_sample = "%s[-_]?2\.f.{1,5}$" % ss_id
    mate1 = [x for x in files_in_reads_folder if re.findall(pattern_mate1_this_sample, x)]
    mate2 = [x for x in files_in_reads_folder if re.findall(pattern_mate2_this_sample, x)]

    mate1_f = mate1[0] if mate1 else ""
    mate2_f = mate2[0] if mate2 else ""

    if not (mate1_f or mate2_f):
        raise FileNotFoundError("Error@finding mate files : seems no mate file ....")

    return mate1_f, mate2_f


def safe_split_knife_report_file_line(line):
    return line.strip("\n").split("\t")


class ReportChecker(object):
    def __init__(self, glm_pv=PV_GLM_DEFAULT, naive_pv=PV_NAIVE_DEFAULT, decoy_ratio=DECOY_RATIO_DEFAULT):
        self.glm_pv_lower_limit = glm_pv
        self.naive_pv_lower_limit = naive_pv
        self.decoy_ratio_up_limit = decoy_ratio

    def check_naive_report_line(self, line_in):
        try:
            parts = safe_split_knife_report_file_line(line_in)
            r_circ, r_decoy, r_pv = parts[5], parts[6], parts[7]
        except:
            _logger.error("encounter an error at : {}".format(line_in))
            return False

        try:
            circ, decoy, pv = int(r_circ), int(r_decoy), float(r_pv)
        except ValueError:
            _logger.log(0, "meet a - at %s" % "\t".join(parts))
            return False
        return pv >= self.naive_pv_lower_limit and decoy / (circ + 0.0) < self.decoy_ratio_up_limit

    def check_glm_report_line(self, line_in):
        try:
            parts = safe_split_knife_report_file_line(line_in)
            pv = float(parts[2])
        except:
            _logger.error("encounter an error at : {}".format(line_in))
            return False
        return pv >= self.glm_pv_lower_limit

    @staticmethod
    def is_this_bsj_str_circular(line_in):
        try:
            info_str = safe_split_knife_report_file_line(line_in)[0]
            seq_name, gene_splice_1, gene_splice_2, junction_type, strand = info_str.strip().split("|")
        except:
            return False
        if junction_type == "reg":
            return False
        return True


def _convert_naive_report(naive_report, line_checker):
    _logger.debug("converting naive report : %s" % naive_report)

    with open(naive_report) as nr:
        bsj_lines = nr.readlines()[2:]

        _logger.debug("there is {num_lines} raw bsj entries in {report}".format(num_lines=len(bsj_lines),
                                                                                report=naive_report))

        bed_lines = [bsj_line_to_bed(line) for line in bsj_lines
                     if line_checker.is_this_bsj_str_circular(line) and line_checker.check_naive_report_line(line)]

    return bed_lines


def _convert_glm_report(glm_report, line_checker):
    _logger.debug("converting glm report : %s" % glm_report)
    with open(glm_report) as nr:
        bsj_lines = nr.readlines()[1:]

        _logger.debug("there is {num_lines} raw bsj entries in {report}".format(num_lines=len(bsj_lines),
                                                                                report=glm_report))

        bed_lines = [bsj_line_to_bed(line) for line in bsj_lines
                     if line_checker.is_this_bsj_str_circular(line) and line_checker.check_glm_report_line(line)]

    return bed_lines


def bsj_line_to_bed(line):
    parts = safe_split_knife_report_file_line(line)
    return _bsj_junction_to_bed(parts[0])


def _bsj_junction_to_bed(info_str):
    """junction: chr|gene1_symbol:splice_position|gene2_symbol:splice_position|junction_type|strand
        junction types are reg (linear),
        rev (circle formed from 2 or more exons),
        or dup (circle formed from single exon)
    """

    seq_name, gene_splice_1, gene_splice_2, junction_type, strand = info_str.strip().split("|")
    if junction_type == "reg":
        return None
    else:
        gene1, splice_1 = gene_splice_1.strip().split(":")
        gene2, splice_2 = gene_splice_2.strip().split(":")

        start_point = splice_1 if int(splice_1) < int(splice_2) else splice_2
        end_point = splice_2 if int(splice_1) < int(splice_2) else splice_1

        # name_bsj = "{chr}_{start}_{end}_{gene_from}_{gene_to}".format(chr=seq_name,
        #                                                               start=start_point,
        #                                                               end=end_point,
        #                                                               gene_from=gene1,
        #                                                               gene_to=gene2)

        name_bsj = info_str.strip()
        return "\t".join([seq_name, start_point, end_point, name_bsj, "0", strand])


def _extract_bed_from_knife_report_path(output_bed_file_path, path_of_knife_result, line_checker):
    report_path = os.path.join(path_of_knife_result, "circReads")
    naive_report_folder = os.path.join(report_path, "reports")
    annotated_junction_report_folder = os.path.join(report_path, "glmReports")

    all_bed_lines = []
    is_naive_report_included = False
    is_glm_report_included = False
    for report in os.listdir(naive_report_folder):
        current_naive_report_file = os.path.join(naive_report_folder, report)
        all_bed_lines.extend(_convert_naive_report(current_naive_report_file, line_checker))
        is_naive_report_included = True

    if is_naive_report_included:
        _logger.debug("naive report path: %s" % naive_report_folder)
    else:
        _logger.warning("NO DENOVO REPORT : {}".format(naive_report_folder))

    for report in os.listdir(annotated_junction_report_folder):
        current_glm_report_file = os.path.join(annotated_junction_report_folder, report)
        all_bed_lines.extend(_convert_glm_report(current_glm_report_file, line_checker))
        is_glm_report_included = True

    if is_glm_report_included:
        _logger.debug("glm report path : %s" % annotated_junction_report_folder)
    else:
        _logger.warning("NO GLM REPORT: {}".format(annotated_junction_report_folder))

    _logger.info("after all , there is {} putative bsj junctions".format(len(all_bed_lines)))

    with open(output_bed_file_path, "w") as op:
        for line in all_bed_lines:
            op.write("{}\n".format(line.strip()))


def _check_opts(opts=None):
    temp_error = "Error@KNIFE: %s"
    oc = pysrc.body.option_check.OptionChecker(opts, name=SECTION_DETECT)
    oc.must_have(_OPT_BASH_BIN, pysrc.body.utilities.which,
                 FileNotFoundError(temp_error % "bash binary not found"),
                 "abs path to bash binary")

    oc.must_have(_OPT_KNIFE_SCRIPT, os.path.exists,
                 FileNotFoundError(temp_error % "incorrect KNIFE script path given "),
                 "knife_script: absolute path to KNIFE executive script")

    oc.must_have(_OPT_READ_DIRECTORY, _check_valid_read_directory,
                 FileNotFoundError(temp_error % "incorrect read_directory"),
                 """absolute path to directory containing fastq files for alignment.
                 Paired-end reads (PE) must have read1 and read2 in separate files.""")

    oc.may_need(_OPT_READ_ID_STYLE, lambda x: x.strip() in ["appended", "complete"],
                KeyError(temp_error % "incorrect read_id_type given"),
                """read_id_style, complete|appended (use complete for single end).""")

    oc.must_have(_OPT_ALIGNMENT_PARENT_DIRECTORY, _is_this_folder_existing,
                 NotADirectoryError(temp_error % "must have a directory for all result "),
                 """alignment_parent_directory:
                 absolute path to directory where the dataset analysis output and log files will be stored.
                 This directory must already exist,
                 and a directory named dataset_name (see below) will be created under this directory
                  for all output files.""")

    oc.must_have(_OPT_DATA_SET_NAME, lambda x: True,
                 NotADirectoryError(temp_error % "incorrect dataset_name given"),
                 """string identifier for this dataset.
                 A folder of this name will be created under alignment_parent_directory (see above)
                 and all output for this run will be stored in this directory.""")

    oc.must_have(_OPT_JUNCTION_OVERLAP, lambda x: x.isdecimal(),
                 KeyError(temp_error % "junction_overlap should be a integer"),
                 """minimum number of bases in the read
                 which must be on each side of the junction to consider
                 that the read is overlapping the junction.
                  Values that have empirically worked well are 8 for paired-end (PE) reads of length < 70,
                  13 for longer PE,
                  and 10 for single-end (SE) reads of length < 70,
                  15 for longer SE reads.""")

    oc.may_need(_OPT_INDEX_PATH, _is_this_folder_existing,
                NotADirectoryError(temp_error % "incorrect KNIFE pre-fabricated index path"),
                """index_path: a path to KNIFE index directory. """)

    oc.may_need(_OPT_BED_OUTPUT, pysrc.body.utilities.is_path_creatable,
                FileNotFoundError(temp_error % "can not create a bed file there"),
                des_str="""path to filtered .bed format report""")

    oc.may_need(_OPT_GLM_PV, pysrc.body.utilities.is_ratio_like,
                ValueError(temp_error % "glm pv criteria should be a float between 0 and 1"),
                des_str="""a numeric ratio filters the GLM model BSJ""")

    oc.may_need(_OPT_NAIVE_PV, pysrc.body.utilities.is_ratio_like,
                ValueError(temp_error % "naive pv criteria should be a float between 0 and 1"),
                des_str="""a numeric ratio filters the naive BSJ result""")

    oc.may_need(_OPT_DECOY_RATIO, pysrc.body.utilities.is_ratio_like,
                ValueError(temp_error % "decoy ratio criteria should be a float between 0 and 1"),
                des_str="""a numeric ratio filters the decoy ratio of naive BSJ result""")

    return oc


opt_checker = _check_opts()

OPTION_CHECKERS = [opt_checker]


def detect(par_dict=None, **kwargs):
    opts_of_index_phase_raw = pysrc.body.cli_opts.merge_parameters(kwargs, par_dict, SECTION_DETECT)
    opts = copy.copy(opts_of_index_phase_raw)

    opts = _predict_id_style_and_overlap_length(opts)

    opt_checker.check(copy.copy(opts_of_index_phase_raw))

    path_of_knife_result = os.path.join(opts[_OPT_ALIGNMENT_PARENT_DIRECTORY],
                                        opts[_OPT_DATA_SET_NAME])

    path_output_bed = pysrc.body.cli_opts.extract_one(opts, _OPT_BED_OUTPUT,
                                                      default=os.path.join(path_of_knife_result, "main.bed"))

    report_line_checker = ReportChecker(glm_pv=pysrc.body.cli_opts.extract_one(opts, _OPT_GLM_PV, PV_GLM_DEFAULT),
                                        naive_pv=pysrc.body.cli_opts.extract_one(opts, _OPT_NAIVE_PV, PV_NAIVE_DEFAULT),
                                        decoy_ratio=pysrc.body.cli_opts.extract_one(opts, _OPT_DECOY_RATIO,
                                                                                    DECOY_RATIO_DEFAULT))

    # start invoking KNIFE
    cmd_detect = _get_detect_cmd(opts)
    pysrc.body.worker.run(cmd_detect)  # KNIFE running here

    # post-detection transforming
    _extract_bed_from_knife_report_path(output_bed_file_path=path_output_bed,
                                        path_of_knife_result=path_of_knife_result,
                                        line_checker=report_line_checker)
    return opts_of_index_phase_raw


def export_as_bed(par_dict=None, **kwargs):
    # todo: should change into a 3rd party function
    opts_of_index_phase_raw = pysrc.body.cli_opts.merge_parameters(kwargs, par_dict, SECTION_DETECT)

    knife_opts_dict = copy.copy(opts_of_index_phase_raw)

    # use combined-report as primary source
    path_of_knife_result = os.path.join(knife_opts_dict[_OPT_ALIGNMENT_PARENT_DIRECTORY],
                                        knife_opts_dict[_OPT_DATA_SET_NAME])

    path_output_bed = knife_opts_dict.get(_OPT_BED_OUTPUT, os.path.join(path_of_knife_result, "main.bed"))

    report_line_checker = ReportChecker(glm_pv=knife_opts_dict[_OPT_GLM_PV],
                                        naive_pv=knife_opts_dict[_OPT_NAIVE_PV],
                                        decoy_ratio=knife_opts_dict[_OPT_DECOY_RATIO])

    _extract_bed_from_knife_report_path(output_bed_file_path=path_output_bed,
                                        path_of_knife_result=path_of_knife_result,
                                        line_checker=report_line_checker)
