# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import Bio.SeqIO as Bo

__doc__ = '''
'''
__author__ = 'zerodel'

FASTQ_FILE_EXTENSION = [".fastq", ".fq"]


def is_pair_end_fastq_id_identical(fq1_path, fq2_path):
    num_id = 3
    counter = 0
    fq_entries = zip(Bo.parse(fq1_path, "fastq"), Bo.parse(fq2_path, "fastq"))

    for r1, r2 in fq_entries:
        counter += 1
        if counter > num_id:
            return True
        # or (str(r1.description).strip() != str(r2.description).strip()):
        if r1.id != r2.id:
            return False

    if 0 == counter:
        raise KeyError(
            "Error@fq_checking_id_style: empty fastq files or wrong file format.")
    else:
        return True


def get_read_length(fq):
    # here we assume all reads in a single fastq file shares the same length .
    return len(next(Bo.parse(fq, "fastq")).seq)


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print(__doc__)
    else:
        pass
