# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import argparse
import os
import os.path
import sys

upper_root = os.path.abspath(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

sys.path.append(upper_root)

import pysrc.sub_module.summary_quant as sq

parser = argparse.ArgumentParser()
parser.add_argument("anno", help="""path to annotation file or quantification result directory with each sample has a final.gtf in
                                 its sub-dir""",
                    default="",
                    nargs="+")
parser.add_argument("-o", help="output file location", default="")

__doc__ = ''' summarize the host gene information from 
'''

__author__ = 'zerodel'

if __name__ == "__main__":
    args = parser.parse_args()
    if args.o and args.anno:
        isoform_ownership = sq.TranscriptOwnership()

        for entry_annotation in args.anno:
            if os.path.isfile(entry_annotation):
                isoform_ownership.parse_gtf(entry_annotation)
            elif os.path.isdir(entry_annotation):
                # print(entry_annotation)

                gtf_under_this_folder = [os.path.join(entry_annotation, x, "final.gtf") for x in os.listdir(entry_annotation) if
                      os.path.isdir(os.path.join(entry_annotation, x)) and not x.startswith(".")]
                # for x in os.listdir(entry_annotation):
                #     print(x)
                #     if os.path.isdir(os.path.join(entry_annotation, x)) and not x.startswith("."):
                #         print(os.path.join(entry_annotation, x, "final.gtf"))

                # print(gtf_under_this_folder)
                for path_gtf in gtf_under_this_folder:
                    if os.path.exists(path_gtf) and os.path.isfile(path_gtf):
                        isoform_ownership.parse_gtf(path_gtf)

            else:
                print("a wrong path is given : {}".format(entry_annotation))
                continue
        with open(args.o, "w") as dump_it:
            for line in isoform_ownership.to_text_table_lines():
                dump_it.write("{}\n".format(line.strip()))
    else:
        print("wrong arguments given , please check usage by using '-h' option")
        exit(-1)
