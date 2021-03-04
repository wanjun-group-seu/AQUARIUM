# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import argparse
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(os.path.realpath(__file__)))))

import pysrc.sub_module.summary_quant
import pysrc.body.logger

__doc__ = ''' this script is to test the ownership class in pysrc.summary_quant module  , 
try use -h option to get more information 
'''

__author__ = 'zerodel'

_logger = pysrc.body.logger.default_logger("TEST_TRANSCRIPT_OWNERSHIP")


def get_args_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--anno", help="path to .gtf file", default="")
    parser.add_argument("-c", "--ciri", help="path to ciri report file , v1.2", default="")
    parser.add_argument("-o", "--output", help="path to output table file", default="")
    return parser


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print(__doc__)

    else:
        args = get_args_parser().parse_args()

        _logger.info("your args are {args_str}".format(args_str=str(args)))

        test_obj = pysrc.sub_module.summary_quant.TranscriptOwnership()
        if args.anno:
            test_obj.parse_gtf(args.anno)
            _logger.info("success in parsing annotation file : {gtf}".format(gtf=args.anno))
        else:
            _logger.warning("no annotation given !")

        if args.ciri:
            test_obj.parse_ciri(args.ciri)
            _logger.info("success in parsing ciri output: {ciri_report}".format(ciri_report=args.ciri))

        with open(args.output, "w") as dump_it:
            _logger.info("start dump info into {output_file}".format(output_file=args.output))
            for obj_str_line in test_obj.to_text_table_lines():
                dump_it.write("{}\n".format(obj_str_line.strip()))

        _logger.info("all finished.... ")
