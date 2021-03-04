# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import os
import sys
import pysrc.file_format.bsj_gtf
import argparse

__doc__ = ''' make a annotation file for bed file containing circular RNA identification
'''
__author__ = 'zerodel'


def __this_args():
    args_parser = argparse.ArgumentParser()
    args_parser.add_argument("-g", "--gtf", help="file path to a genomic annotation file in .gtf format", default="")
    args_parser.add_argument("-c", "--circ", help="file path to a circular RNA identification file, in .bed file "
                                                  "format", default="")
    args_parser.add_argument("-o", "--output", help="output file path", default="")
    return args_parser


def check_your_args(args):
    return os.path.exists(args.gtf) and os.path.exists(args.circ)


if __name__ == "__main__":
    import sys

    if 1 == len(sys.argv):
        print(__doc__)
    else:
        args_parser = __this_args()
        args = args_parser.parse_args()

        if check_your_args(args):
            pysrc.file_format.bsj_gtf.do_make_gtf_for_circular_prediction_greedy(gff_db=args.gtf,
                                                                                 circular_candidate_regions=args.circ,
                                                                                 output_gtf_path_name=args.output)
