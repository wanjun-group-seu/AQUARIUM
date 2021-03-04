# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import os

import pysrc.body.logger

__doc__ = '''
'''
__author__ = 'zerodel'

_logger = pysrc.body.logger.default_logger("sam")

try:
    import pysam
except ImportError:
    _logger.error("ImportERROR: unable to load pysam module")
    raise ImportError("unable to load pysam")


class AlignEntry(object):
    def __init__(self, str_line_sam):
        pass


def is_valid_sam(sam_file):
    sub_dir, sam = os.path.split(sam_file)
    if os.path.isdir(sub_dir) and os.path.exists(sub_dir):
        if sam.endswith(".sam"):
            return True
    else:
        return False


def is_sam_from_bwa(path_sam):
    def _is_alignment_from_bwa(align_sam):
        with open(align_sam) as sr:
            for line in sr:
                if line.strip().startswith("@"):
                    if line.strip().startswith("@PG"):
                        parts = line.strip().split()
                        for part in parts:
                            if part.startswith("ID:"):
                                return part.split(":")[-1].lower() == "bwa"
                elif len(line.strip()) > 2:
                    return False

    return os.path.exists(path_sam) and path_sam.endswith(".sam") and _is_alignment_from_bwa(path_sam)
