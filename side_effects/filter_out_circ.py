# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import pysam
import Bio.SeqIO


__doc__ = '''pick out the circular reads and reference
'''
__author__ = 'zerodel'


def pick_out_circular_reference(genomic_in, circular_out):
    with open(circular_out, "w") as dump_circ:
        for fa in Bio.SeqIO.parse(genomic_in, "fasta"):
            if fa.id.strip().startswith("chr"):
                dump_circ.write(">%s\n%s\n" % (fa.id, fa.seq.strip()))


def pick_out_reads_circular(fa_in, circ_out):
    with open(circ_out, "w") as dump_circ:
        for fa in Bio.SeqIO.parse(fa_in, "fasta"):
            read_num, seq_id = fa.id.strip().split("/")
            if seq_id.strip().startswith("chr"):
                dump_circ.write(">%s\n%s\n" % (fa.id, fa.seq.strip()))
