# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import os
import sys

import pysrc.wrapper.gffread
import pysrc.file_format.gtf

__doc__ = '''this file contains utilities works for lincRNA functionality
'''

__author__ = 'zerodel'

__ANNOTATION_SECTION_BIOTYPE = "gene_biotype"
__BIOTYPE_STR_LINC = ["lincRNA", "lncRNA"]

__linc_sign = 'gene_biotype "lincRNA"'


def prepare_linc_annotation(original_gff, target_linc_annotation):
    with open(target_linc_annotation, "w") as dump_it:
        dump_it.writelines(pysrc.file_format.gtf.filter_gff_by_source(gff=original_gff, bio_type=__BIOTYPE_STR_LINC))


def prepare_linc_transcriptome_seq(linc_annotation, genomic_seq, target_fa):
    pysrc.wrapper.gffread.do_extract_non_coding_transcript(gff=linc_annotation,
                                                           path_ref_sequence_file=genomic_seq,
                                                           output=target_fa)


if __name__ == "__main__":
    pass
