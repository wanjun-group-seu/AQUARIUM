# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import os
import argparse

import pysrc.body.config
import pysrc.body.logger
import pysrc.wrapper.knife
import pysrc.file_format.gtf

_SUB_DIR_COMBINED_REPORT = "circReads/combinedReports"
_SUB_DIR_READ1_GLM_REPORT = "circReads/glmReports"
_SUB_DIR_READ1_DENOVO_REPORT = "circReads/reports"

BED_FORMAT_OUTPUT_FILE_NAME = "out.bed"

__doc__ = ''' this script will transform KNIFE output into a .bed file under the same path of KNIFE folder
'''
__author__ = 'zerodel'

WORKFLOW_NAME = "TRANSFORM_KNIFE_REPORT"

_logger = pysrc.body.logger.default_logger(WORKFLOW_NAME)


def load_mapping_info_from_gtf(gtf_in):
    mapping = {}
    _logger.debug("starting loading {gtf}".format(gtf=gtf_in))
    with open(gtf_in) as gtf:
        for line in gtf:
            some_entry = pysrc.file_format.gtf.GTFItem(line.strip())
            attr_this = some_entry.get_attribute()
            if "gene_name" in attr_this and "gene_id" in attr_this:
                mapping[attr_this.get("gene_name")] = attr_this.get("gene_id")

    _logger.debug("finish loading {gtf}".format(gtf=gtf_in))
    return mapping


def _safe_split_knife_report_file_line(line):
    return line.strip("\n").split("\t")


def _process_line(line_in_file, processed_result_list, func_check_positive, gene_id_mapping=None):
    parts = _safe_split_knife_report_file_line(line_in_file)
    if parts:
        if func_check_positive(parts):
            line_bed = _bsj_junction_to_bed(parts[0], gene_id_mapping)
            if line_bed:
                processed_result_list.append(line_bed)


def _is_this_r1_only_denovo_bsj_positive(parts):
    try:
        r_circ, r_decoy, r_pv = parts[5], parts[6], parts[7]
    except Exception as e:
        raise e

    try:
        circ, decoy, pv = int(r_circ), int(r_decoy), float(r_pv)
    except ValueError:
        return False

    # this is from KNIFE github page
    return pv >= 0.9 and decoy < circ * 0.1


def _is_this_r1_only_glm_bsj_positive(parts):
    if not parts:
        return False
    pv = float(parts[2])
    return pv >= 0.9  # this is also from KNIFE github page


def _is_this_combined_glm_bsj_positive(parts):
    if not parts:
        return False
    threshold_posterior = 0.9
    p_orig = float(parts[2])
    p_swap = float(parts[4])
    return p_orig >= threshold_posterior or p_swap >= threshold_posterior


def _is_this_combined_denovo_bsj_positive(parts):
    """ junction	orig_circOrLinear	orig_decoyOrAnom	orig_unmapped	orig_pval	swapped_circOrLinear	swapped_decoyOrAnom	swapped_unmapped	swapped_pval	total_reads
    """
    if not parts:
        return False

    pv_threshold = 0.9
    reads_rate_threshold = 0.1
    orig_circOrLinear, orig_decoyOrAnom, orig_pval = int(parts[1]), int(parts[2]), float(parts[4])
    swapped_circOrLinear, swapped_decoyOrAnom, swapped_pval = int(parts[5]), int(parts[6]), float(parts[8])

    by_r1 = orig_pval >= pv_threshold  # and orig_decoyOrAnom < orig_circOrLinear * reads_rate_threshold
    by_r2 = swapped_pval >= pv_threshold  # and swapped_decoyOrAnom < swapped_circOrLinear * reads_rate_threshold
    return by_r1 or by_r2


def _bsj_junction_to_bed(info_str, gene_id_mapping=None):
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

        gene_name = gene1 if gene1 == gene2 else "n/a"
        field_gene = gene_id_mapping.get(gene_name, "n/a") if gene_id_mapping else gene_name

        name_bsj = "{chr}:{start_nt}|{end_nt}@{gene}".format(chr=seq_name, start_nt=str(int(start_point) + 1),
                                                             end_nt=end_point,
                                                             gene=field_gene)
        return "\t".join([seq_name, start_point, end_point, name_bsj, "0", strand])


def __dump_bed_line_list_to(all_bed_lines, output_bed_file_path):
    with open(output_bed_file_path, "w") as op:
        for line in all_bed_lines:
            op.write("{}\n".format(line.strip()))


def pick_glm_report(file_name):
    return file_name.strip().endswith("circJuncProbs.txt")


def pick_denovo_report(file_name):
    return file_name.strip().endswith("report.txt") and "unaligned" not in file_name


def __get_file_under(tap_root, sub_dir, func_to_pick_report):
    report_sub_dir = os.path.join(tap_root, sub_dir)
    if os.path.exists(report_sub_dir) and os.path.isdir(report_sub_dir):
        contents_file = os.listdir(report_sub_dir)
        if contents_file:
            file_you_need = [x.strip() for x in contents_file if func_to_pick_report(x.strip())][0]
            return os.path.join(report_sub_dir, file_you_need)
        else:
            raise FileNotFoundError("no files under this folder: {}".format(report_sub_dir))
    else:
        raise NotADirectoryError("no such folder: {}.".format(report_sub_dir))


def transform_each_line(report_file, func_to_check_bsj, header_line_num, gene_id_mapping=None):
    bed_lines_this_file = []
    with open(report_file) as report_it:
        for x in range(header_line_num):
            report_it.readline()
        for line in report_it:
            _process_line(line, bed_lines_this_file, func_to_check_bsj, gene_id_mapping)
    return bed_lines_this_file


def __report_glm_combined(dir_knife, gene_id_mapping=None):
    # here we assume the combined glm report is that file with suffix: circJuncProbs.txt under the
    # circReads/combinedReports folder .
    report_file = __get_file_under(dir_knife, _SUB_DIR_COMBINED_REPORT, pick_glm_report)

    _logger.debug("now loading {}".format(report_file))

    bed_lines_this_file = transform_each_line(report_file, _is_this_combined_glm_bsj_positive, 1, gene_id_mapping)

    return bed_lines_this_file


def __report_glm_read1_only(dir_knife, gene_id_mapping=None):
    report_file = __get_file_under(dir_knife, _SUB_DIR_READ1_GLM_REPORT, pick_glm_report)
    _logger.debug("now loading {}".format(report_file))
    bed_outputs = transform_each_line(report_file, _is_this_r1_only_glm_bsj_positive, 1, gene_id_mapping)
    return bed_outputs


def __report_denovo_combined(dir_knife, gene_id_mapping=None):
    report_file = __get_file_under(dir_knife, _SUB_DIR_COMBINED_REPORT, pick_denovo_report)
    _logger.debug("now loading {}".format(report_file))
    bed_contents = transform_each_line(report_file, _is_this_combined_denovo_bsj_positive, 1, gene_id_mapping)
    return bed_contents


def __report_denovo_read1_only(dir_knife, gene_id_mapping=None):
    report_file = __get_file_under(dir_knife, _SUB_DIR_READ1_DENOVO_REPORT, pick_denovo_report)
    _logger.debug("now loading {}".format(report_file))
    bed_lines = transform_each_line(report_file, _is_this_r1_only_denovo_bsj_positive, 2, gene_id_mapping)
    return bed_lines


def __infer_knife_output_directory(cfg_file):  # todo: here we confirmed that this is a configuration file
    if os.path.exists(cfg_file) and os.path.isfile(cfg_file):
        user_config = pysrc.body.config.config(cfg_file)
        dir_par = user_config[pysrc.wrapper.knife.SECTION_DETECT][pysrc.wrapper.knife._OPT_ALIGNMENT_PARENT_DIRECTORY]
        dir_object = user_config[pysrc.wrapper.knife.SECTION_DETECT][pysrc.wrapper.knife._OPT_DATA_SET_NAME]
        putative_dir = os.path.join(dir_par, dir_object)
        if os.path.exists(putative_dir) and os.path.isdir(putative_dir):
            return putative_dir
        else:
            raise NotADirectoryError("Unable to get the right path from : {}.".format(cfg_file))
    else:
        raise FileNotFoundError("Unable to find the cfg_file in : {}.".format(cfg_file))


def main(cfg_path, is_conserve, is_use_swapped, anno_file):
    dir_knife_output = __infer_knife_output_directory(cfg_path)

    gene_name_id_mapping = load_mapping_info_from_gtf(anno_file) if anno_file and os.path.exists(anno_file) else None

    if is_conserve:  # use only glm report
        if is_use_swapped:  # glm_combined
            bed_lines = __report_glm_combined(dir_knife_output, gene_name_id_mapping)
        else:  # glm_r1
            bed_lines = __report_glm_read1_only(dir_knife_output, gene_name_id_mapping)
    else:
        if is_use_swapped:  # glm_combined, and denovo_combined
            bed_lines = __report_glm_combined(dir_knife_output, gene_name_id_mapping)
            bed_lines.extend(__report_denovo_combined(dir_knife_output, gene_name_id_mapping))
        else:  # glm_r1 and denovo_r1
            bed_lines = __report_glm_read1_only(dir_knife_output, gene_name_id_mapping)
            bed_lines.extend(__report_denovo_read1_only(dir_knife_output, gene_name_id_mapping))

    bed_path = os.path.join(dir_knife_output, BED_FORMAT_OUTPUT_FILE_NAME)
    __dump_bed_line_list_to(bed_lines, bed_path)
    return bed_path


def __cli_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("cfg_file", help="path to the config file, KNIFE setting is required")
    parser.add_argument("-c", "--conserve", help="use annotation powered glm report only", action="store_true")
    parser.add_argument("-s", "--use_swapped", help="use both reads instead only read 1", action="store_true")
    parser.add_argument("-a", "--anno", help="path for annotation file , which provide the name id mapping info",
                        default="")
    return parser


if __name__ == "__main__":
    arg_parser = __cli_arg_parser()
    args = arg_parser.parse_args()
    _logger.debug("config_file is {}".format(args.cfg_file))
    _logger.debug("use glm only : {}".format(args.conserve))
    _logger.debug("use combined report : {}".format(args.use_swapped))
    _logger.debug("user defined annotation file: {}".format(args.anno))
    print(main(args.cfg_file, args.conserve, args.use_swapped, args.anno))
