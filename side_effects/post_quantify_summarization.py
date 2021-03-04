# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import argparse
import os
import sys

upper_root = os.path.abspath(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

sys.path.append(upper_root)

import pysrc.sub_module.summary_quant as sq


parser = argparse.ArgumentParser()
parser.add_argument("annotation", help="path to annotation file", default="")
parser.add_argument("-q", help="path to quantification file", default="")
parser.add_argument("-o", help="output file location", default="")

__doc__ = '''
'''

__author__ = 'zerodel'

if __name__ == "__main__":
    args = parser.parse_args()

    if args.annotation and args.q and args.o:
        sq.aggregate_isoform_quantify_result(quant_sf=args.q, summarized_output=args.o, gtf_annotation=args.annotation)

    else:
        if args.annotation and args.o:
            isoform_ownership = sq.TranscriptOwnership()
            isoform_ownership.parse_gtf(args.annotation)
            if args.o:
                with open(args.o, "w") as dump_it:
                    for line in isoform_ownership.to_text_table_lines():
                        dump_it.write("{}\n".format(line.strip()))
            else:
                # use stdout
                for line in isoform_ownership.to_text_table_lines():
                    print("{}\n".format(line.strip()))
