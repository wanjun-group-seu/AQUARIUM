# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import functools
import os

import Bio.Seq
import Bio.SeqIO

import pysrc.body.logger
import pysrc.file_format

__doc__ = '''
'''
__author__ = 'zerodel'

FASTA_FILE_EXTENSION = [".fasta", ".fa", ".fna", ".csfasta", "csfa"]

_logger = pysrc.body.logger.default_logger("FASTA_FILE_OPERATION")


def convert_all_entries_in_fasta(fa_in, fa_out, convert_fun):
    res = [convert_fun(fa) for fa in Bio.SeqIO.parse(fa_in, "fasta")]
    with open(fa_out, "w") as fa_export:
        Bio.SeqIO.write(res, fa_export, "fasta")


def make_adapter(k):
    def fun_inner(fa, kmer_len):
        kmer_len_fixed = kmer_len - 1 if kmer_len >= 1 else 0
        fa_seq = str(fa.seq)
        fa_seq = "{}{}".format(fa_seq[-kmer_len_fixed:], fa_seq)
        fa.seq = Bio.Seq.Seq(fa_seq)
        return fa

    return functools.partial(fun_inner, kmer_len=k)


def pad_for_effective_length(length_needed):
    def fun_inner(fa, len_you_need):
        fa_seq = str(fa.seq)
        fa_seq = "{}{}".format("N" * len_you_need, fa_seq)
        fa.seq = Bio.Seq.Seq(fa_seq)
        return fa

    return functools.partial(fun_inner, len_you_need=length_needed)


def is_fasta(ref_path):
    basename, extension_with_dot = os.path.splitext(ref_path)
    return extension_with_dot in pysrc.file_format.fa.FASTA_FILE_EXTENSION


def incremental_updating(target_fa, list_of_fa_file):
    if list_of_fa_file:
        _logger.debug("we will combine the following files : {}".format('\n'.join(list_of_fa_file)))
        whole_dict_fa = {}
        for single_fa_file in list_of_fa_file:
            if not os.path.exists(single_fa_file):
                _logger.warning("WARNING: not valid fa file: {}".format(single_fa_file))

            else:
                for x in Bio.SeqIO.parse(single_fa_file, "fasta"):
                    whole_dict_fa[x.id.strip()] = str(x.seq).strip()

        # output to target_fa,
        with open(target_fa, "a") as dump_fa:
            for x in whole_dict_fa:
                dump_fa.write(">{name}\n{sequence}\n".format(name=x,
                                                             sequence=whole_dict_fa[x]))
    else:
        _logger.warning("list contains no fa file , so no action taken ")
        import pathlib
        pathlib.Path(target_fa).touch()


def pool_fa_by_seq_content(path_target, lst_fa):
    # seq using id in latter fa files of lst_fa
    if lst_fa:
        _logger.debug("now pool those fa together, {}".format(
            "\n".join(lst_fa)))
        dict_seq = {}
        for one_fa in lst_fa:
            if not os.path.exists(one_fa):
                _logger.warning(" no such file in {}".format(one_fa))
            else:
                for x in Bio.SeqIO.parse(one_fa, "fasta"):
                    dict_seq[str(x.seq).strip()] = x.id.strip()

        dict_seq_final = {v: k for k, v in dict_seq.items()}
        with open(path_target, "a") as dump_fa:
            for x in dict_seq_final:
                dump_fa.write(">{name}\n{sequence}\n".format(name=x,
                                                             sequence=dict_seq_final[x]))

    else:
        _logger.warning("list of fa has NONE content , no action taken")
        import pathlib
        pathlib.Path(path_target).touch()


def filter_fa_by_id(path_raw, path_target, id_lst):

    with open(path_target, "a") as dump_fa:
        for x in Bio.SeqIO.parse(path_raw, "fasta"):
            if x.id.strip() not in id_lst:
                dump_fa.write(">{name}\n{sequence}\n".format(name=x.id.strip(),
                                                             sequence=str(x.seq).strip()))
