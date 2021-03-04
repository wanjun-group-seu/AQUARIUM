# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import os
import os.path
import argparse

import pysam
import Bio.SeqIO

__doc__ = '''
'''
__author__ = 'zerodel'


def __get_parser():
    parser_your_args = argparse.ArgumentParser()
    parser_your_args.add_argument("--ciri", "-C", default="", help="file path to CIRI output")
    parser_your_args.add_argument("--polyester", "-P", default="",
                                  help="file path to polyester simulation setting file")
    return parser_your_args


def strip_junction_str(str_junction):
    if "." in str_junction:
        return str_junction.strip().split(".")[0]
    else:
        return str_junction


def __load_ciri_output(ciri_output_file):
    bsj_and_their_reads = {}
    with open(ciri_output_file) as ciri:
        ciri.readline()  # here we assume the header of ciri should be skipped
        for line in ciri:
            bsj_str = strip_junction_str(line.strip().split()[0])
            circ_reads = line.strip().split()[-1].strip().strip(",").split(",")
            bsj_and_their_reads[bsj_str] = circ_reads

    return bsj_and_their_reads


def __load_polyester_setting(polyester_setting_file, is_all_transcript_include=False):
    bsj_and_their_number, bsj_and_their_length = {}, {}
    # we assume the header of the setting file is "Name    Length  NumRead"
    with open(polyester_setting_file) as psf:
        header = psf.readline().strip().split()
        for line in psf:
            name, length, num_read = line.strip().split()
            name = strip_junction_str(name)
            if is_all_transcript_include or name.startswith("chr"):
                bsj_and_their_length[name] = int(length)
                bsj_and_their_number[name] = int(num_read)

    class NumLen(object):
        def __init__(self, len_dict, num_dict):
            self.len = len_dict
            self.num = num_dict

    return NumLen(bsj_and_their_length, bsj_and_their_number)


def extract_junction_str(bsj_str):
    bsj_str = bsj_str.strip()
    chr_name, start_end = bsj_str.split(":")
    start, end = start_end.split("|")
    return chr_name, int(start), int(end)


def origin_of_simulated_read(read_id):
    num, origin_str = read_id.strip().split("/")
    return strip_junction_str(origin_str)


def pick_out_junction_region(bam_file_path, junction_list):
    with open(bam_file_path) as bam_file_op:
        bam = pysam.AlignmentFile(bam_file_op)

        for bsj in junction_list:
            chr_id, start, end = extract_junction_str(bsj)

            reads_5 = dict([(x.query_name, x) for x in bam.fetch(chr_id, start, start + 20)])
            reads_3 = dict([(x.query_name, x) for x in bam.fetch(chr_id, end - 20, end)])
            reads_id_common = sorted(list(set([x for x in reads_5]) & set([x for x in reads_3])))
            yield (bsj, reads_id_common)


def fetch_total_junction_list(ciri_out, setting_file):
    ciri_mapping = __load_ciri_output(ciri_out)
    poly_setting = __load_polyester_setting(setting_file).num
    total_junctions = sorted(list(set([x for x in ciri_mapping]) | set([x for x in poly_setting])))
    return total_junctions


def make_extract_junction_read_table(ciri_out, setting_file, bam_file, junction_info_table):
    total_junctions = fetch_total_junction_list(ciri_out, setting_file)

    with open(junction_info_table, "w") as summarize_it:
        for bsj, reads_over_bsj in pick_out_junction_region(bam_file, total_junctions):
            if reads_over_bsj and len(reads_over_bsj) > 0:
                for one_read in reads_over_bsj:
                    summarize_it.write("%s\t%s\n" % (bsj, one_read))


def __bsj_total_circ_design(bsj_id, reads_overlap_bsj):
    num_reads_from_circ = 0
    num_reads_designed = 0  # this means this read is from the same junction as the simulation setting
    num_reads_this_junction = 0
    for read_id in reads_overlap_bsj:
        num_reads_this_junction += 1
        read_origin = origin_of_simulated_read(read_id)
        if read_origin.startswith("chr"):
            num_reads_from_circ += 1
        if read_origin == strip_junction_str(bsj_id):
            num_reads_designed += 1
    report_line_for_this_bsj = "{bsj}\t{total}\t{circ}\t{design}\n".format(bsj=bsj_id,
                                                                           total=num_reads_this_junction,
                                                                           circ=num_reads_from_circ,
                                                                           design=num_reads_designed)
    return report_line_for_this_bsj


def make_a_summary_each_junction(ciri_out, setting_poly, bam_file, output_summary):
    total_junctions = fetch_total_junction_list(ciri_out, setting_poly)

    bsj_reads_summarize(bam_file, output_summary, total_junctions)


def bsj_reads_summarize(bam_file, output_summary, total_junctions):
    with open(output_summary, "w") as summarize_it:
        summarize_it.write("bsj\ttotal\tcirc\tdesign\n")

        for bsj, reads_overlap in pick_out_junction_region(bam_file, total_junctions):
            report_for_this_bsj = __bsj_total_circ_design(bsj, reads_overlap)

            summarize_it.write(report_for_this_bsj)


def summary_ciri_reads(ciri_file, output_summary):
    ciri = __load_ciri_output(ciri_file)

    with open(output_summary, "w") as summarize_it:
        summarize_it.write("bsj\ttotal\tcirc\tdesign\n")

        for bsj in ciri:
            summarize_it.write(__bsj_total_circ_design(bsj, ciri[bsj]))


def _piece_together_bsj_seq(genomic_seq_file, bsj_list, k=20):
    # load fa dict

    fas = dict([(fa.id, fa.seq) for fa in Bio.SeqIO.parse(genomic_seq_file, "fasta")])
    for bsj in bsj_list:
        chr_name, start, end = extract_junction_str(strip_junction_str(bsj))

        seq_5 = fas[chr_name][start - 1: start - 1 + k]
        seq_3 = fas[chr_name][end - k: end]
        seq_bsj = ">{seq_name}\n{seq3}{seq5}\n".format(seq_name=strip_junction_str(bsj),
                                                       seq3=seq_3,
                                                       seq5=seq_5)
        yield seq_bsj


def dump_bsj_fa(genomic_seq, bsj_list, k, output_fa):
    with open(output_fa, "w") as out:
        for bsj_junction_seq in _piece_together_bsj_seq(genomic_seq, bsj_list, k):
            out.write(bsj_junction_seq)


def bump_bsj_as_bed(bsj_list, output_bed):
    with open(output_bed) as below_bed:
        for bsj in bsj_list:
            chr_id, start, end = extract_junction_str(strip_junction_str(bsj))
            below_bed.write("%s\t%d\t%d\t%s\n" % (chr_id, start, end, bsj))
