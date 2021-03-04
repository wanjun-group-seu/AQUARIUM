# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import argparse
import Bio.SeqIO

__doc__ = '''
'''
__author__ = 'zerodel'


def __get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fa", "-i", default="", help="input fasta file path")
    parser.add_argument("--fq", "-o", default="", help="output fastq file path")
    return parser


TEMPLATE_FQ = "@{seq_id}\n{seq}\n+\n{quality}\n"


def fa2fq(input_fa, output_fq):
    with open(output_fq, "w") as out:
        for fa in Bio.SeqIO.parse(input_fa, "fasta"):
            fq_entry = TEMPLATE_FQ.format(seq_id=fa.id,
                                          seq=fa.seq.strip(),
                                          quality="g" * len(fa.seq.strip()))
            out.write(fq_entry)


if __name__ == "__main__":
    get_args = __get_parser()
    args = get_args.parse_args()
    fa2fq(args.fa, args.fq)
