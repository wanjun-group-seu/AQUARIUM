# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import os

import pysrc.body
import pysrc.body.logger
import pysrc.body.worker

__doc__ = ''' wrapper of gffread from cufflinks
'''
__author__ = 'zerodel'

_logger = pysrc.body.logger.default_logger("gffread_process")


def build_cmd_gffread(annotation_files, genomic_seqs, output_fasta, transcript_filter="CME", gffread_cmd="gffread"):
    annotation_files = os.path.abspath(annotation_files)
    genomic_seqs = os.path.abspath(genomic_seqs)
    output_fasta = os.path.abspath(output_fasta)

    cmd_parts = [gffread_cmd, annotation_files, "-g", genomic_seqs]

    if transcript_filter:
        cmd_parts.append("-{}".format(transcript_filter))
    else:
        pass

    cmd_parts.append("-w")
    cmd_parts.append(output_fasta)

    return cmd_parts


def do_extract_classic_message_transcript(gff, path_ref_sequence_file, output):
    cmd_parts = build_cmd_gffread(annotation_files=gff, genomic_seqs=path_ref_sequence_file, output_fasta=output)

    _logger.debug("raw cmd for extracting linear transcript is : %s" % str(cmd_parts))

    pysrc.body.worker.run(cmd_parts)


def do_extract_non_coding_transcript(gff, path_ref_sequence_file, output):
    cmd_parts = build_cmd_gffread(annotation_files=gff, genomic_seqs=path_ref_sequence_file, output_fasta=output,
                                  transcript_filter="")

    _logger.debug("raw cmd for extracting circular RNA transcripts : %s" % str(cmd_parts))

    pysrc.body.worker.run(cmd_parts)
