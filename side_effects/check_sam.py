# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import sys
import os

upper_root = os.path.abspath(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

sys.path.append(upper_root)

__doc__ = ''' this is a pilot script for sam file summarization , part of detection prototype
'''

__author__ = 'zerodel'

import argparse
import collections

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pysam


def __get_args_parser():
    this_parser = argparse.ArgumentParser()
    this_parser.add_argument("--sam", default="", help="path to .sam file")
    this_parser.add_argument("--st", default="", help="path to a site table file")
    this_parser.add_argument("--out", default="", help="")
    return this_parser


def __get_sam_obj(path_sam):
    return pysam.AlignmentFile(path_sam)


def load_star_junction_info(path_site_table):
    obj_table = pd.read_table(path_site_table, header=None)
    obj_table.columns = ["chr.donor",
                         "start.donor",
                         "strand.donor",
                         "chr.acceptor",
                         "start.acceptor",
                         'strand.acceptor',
                         "junction.type",
                         "repeat.len.left",
                         "repeat.len.right",
                         "read.name",
                         "start.first",
                         "cigar.first",
                         "start.second",
                         "cigar.second"]

    return obj_table


# def parse_bsj(bam_obj, chr_name, start, end, width=20):
#     half_width = width / 2
#     r5s = [r for r in bam_obj.fetch(chr_name, start - half_width, start + half_width)]
#     r3s = [r for r in bam_obj.fetch(chr_name, end - half_width, end + half_width)]
#
#     pass

# here we use some junction:
# chr12:109046048|109048186


def get_region_str(chr_name, start, end):
    region_str = "{chr}:{start}-{end}".format(chr=chr_name,
                                              start=start,
                                              end=end)
    return region_str


read_stat = collections.namedtuple("read_stat", ["is_multiple",
                                                 "is_unmapped",
                                                 "is_secondary",
                                                 "is_duplicate",
                                                 "is_supplementary",
                                                 "is_primary"])


def is_multiple_mapping(read):
    return read.flag & 0x1


def check_flag(read):
    is_multiple, is_unmapped, is_secondary, is_duplicate, is_supplementary, is_primary = [False] * 6
    flag = read.flag
    if flag & 0x4 == 1:  # unmapped :
        is_unmapped = True
    else:
        if flag & 0x100 == 1:
            is_secondary = True
        if flag & 0x800 == 1:
            is_supplementary = True
    if flag & 0x1 == 1:
        is_multiple = True
    if flag & 0x400:
        is_duplicate = True
    if flag & 0x900 == 0:
        is_primary = True

    return read_stat(is_multiple, is_unmapped, is_secondary, is_duplicate, is_supplementary, is_primary)


def report_site(bam_obj, region_str, func):
    return map(func, bam_obj.fetch(region=region_str))


def get_region_info(bam_obj, chr_name, start, end):
    region_str = "{chr}:{start}-{end}".format(chr=chr_name,
                                              start=start,
                                              end=end)

    return bam_obj.count_coverage(region=region_str)

#    for r in bam_obj.fetch(region = region_str):









# here are when you need to use command line
# if __name__ == "__main__":
#     args = __get_args_parser().parse_args()
#
