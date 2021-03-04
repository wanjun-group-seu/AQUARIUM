# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import os
import sys
import argparse
import random

from collections import namedtuple

import pysrc.body.utilities
import pysrc.file_format.gtf

import Bio.SeqIO
import pysrc.wrapper.gffread
import pysrc.file_format.fa
import pysrc.body.logger


__doc__ = '''decorate circular RNA annotation file with alternative splice events
now supported: exon skipping
'''
__logger = pysrc.body.logger.default_logger("simulation for AS event")
__author__ = 'zerodel'

random.seed(666)

K = 21
READ_LENGTH = 250
LEN_MLL = 250
LEN_FIX_LENGTH = LEN_MLL - K + 1
SEQ_DEPTH = 50

context_of_isoform = namedtuple("context_of_isoform", ["gene_of", "anno_lines"])
info_simu_setting = namedtuple("info_simu_setting", ["Name", "Length", "NumRead"])


def dump_gtf_lines_to(kv_transcript_gtf_lines, out_gff):
    with open(out_gff, "w") as dump_it:
        for rna in kv_transcript_gtf_lines:
            for exon_string in kv_transcript_gtf_lines[rna]:
                dump_it.write("%s\n" % exon_string)


def load_circular_annotation(circ_gff):
    kv_transcript_gene = {}
    kv_transcript_gtf_lines = {}
    with open(circ_gff) as load_gtf:
        for line in load_gtf:
            try:
                exon_entry = pysrc.file_format.gtf.GTFItem(line.strip())

            except pysrc.file_format.gtf.AttributionIncomplete:
                continue

            kv_transcript_gtf_lines.setdefault(exon_entry.get_transcript_id(), set()).add(str(exon_entry))
            kv_transcript_gene[exon_entry.get_transcript_id()] = exon_entry.get_gene_id()

    return context_of_isoform(gene_of=kv_transcript_gene, anno_lines=kv_transcript_gtf_lines)


def __fabricate_exon_skipping_for_non_na_circular_rna(context_of_transcript_isoform):
    gene_of, kv_isoform_gtf_lines = context_of_transcript_isoform
    for rna in [x for x in kv_isoform_gtf_lines.keys()]:
        if not gene_of[rna] == "n/a":
            exons_string_this_rna = kv_isoform_gtf_lines.pop(rna)
            exons_this_rna = sorted([pysrc.file_format.gtf.GTFItem(x) for x in exons_string_this_rna],
                                    key=lambda x: x.starts())

            if len(exons_this_rna) < 3:
                kv_isoform_gtf_lines[rna] = [str(x) for x in exons_this_rna]
                continue

            exon_id_be_skipped = random.choice([1 + i for i in range(len(exons_this_rna) - 2)])

            iso_id_0 = rna + ".0"
            iso_id_1 = rna + ".1"

            exons_iso_0 = []
            exons_iso_1 = []
            for index_exon, exon in enumerate(exons_this_rna):
                neo_exon_iso_0 = pysrc.file_format.gtf.GTFItem(str(exon))
                neo_exon_iso_0.set_transcript_id(iso_id_0)

                neo_exon_iso_1 = pysrc.file_format.gtf.GTFItem(str(exon))
                neo_exon_iso_1.set_transcript_id(iso_id_1)

                exons_iso_0.append(str(neo_exon_iso_0))

                if index_exon == exon_id_be_skipped:
                    # print("pop %d" % index_exon)
                    pass
                else:
                    exons_iso_1.append(str(neo_exon_iso_1))

            kv_isoform_gtf_lines[iso_id_0] = exons_iso_0
            kv_isoform_gtf_lines[iso_id_1] = exons_iso_1
    return kv_isoform_gtf_lines


def prepare_annotation_as_event(path_circular_rna_annotation_file, out_gff):
    context_of_isoform_from = load_circular_annotation(path_circular_rna_annotation_file)
    as_event_added_context = __fabricate_exon_skipping_for_non_na_circular_rna(context_of_isoform_from)
    dump_gtf_lines_to(as_event_added_context, out_gff)


def prepare_seq(linear_gtf, circular_gtf, genomic_seqs, target_folder):
    extractor = pysrc.wrapper.gffread

    path_temp_linear_seq = os.path.join(target_folder, "linear.fa")
    path_temp_circular_seq = os.path.join(target_folder, "circular_raw.fa")

    extractor.do_extract_classic_message_transcript(gff=linear_gtf, path_ref_sequence_file=genomic_seqs,
                                                    output=path_temp_linear_seq)
    extractor.do_extract_non_coding_transcript(gff=circular_gtf, path_ref_sequence_file=genomic_seqs,
                                               output=path_temp_circular_seq)

    path_temp_circular_simu = os.path.join(target_folder, "simu_circular.fa")
    path_temp_circular_quant = os.path.join(target_folder, "quant_circular.fa")

    raw = Bio.SeqIO.parse(path_temp_circular_seq, "fasta")
    with open(path_temp_circular_simu, "w") as fa_simu:
        with open(path_temp_circular_quant, "w") as fa_quant:
            for r in raw:
                trans_seq = str(r.seq).strip()
                seq_simu = "%s%s" % (trans_seq[-READ_LENGTH:], trans_seq)
                seq_quant = "%s%s%s" % ("N" * LEN_FIX_LENGTH, trans_seq[-(K - 1):], trans_seq)
                fa_simu.write(">%s\n%s\n" % (r.id, seq_simu))
                fa_quant.write(">%s\n%s\n" % (r.id, seq_quant))

    combined_fa_simu = os.path.join(target_folder, "simu.fa")
    pysrc.body.utilities.do_merge_files(combined_fa_simu, (path_temp_linear_seq, path_temp_circular_simu))
    combined_fa_quant = os.path.join(target_folder, "quant.fa")
    pysrc.body.utilities.do_merge_files(combined_fa_quant, (path_temp_linear_seq, path_temp_circular_quant))

    os.remove(path_temp_circular_seq)
    os.remove(path_temp_linear_seq)
    os.remove(path_temp_circular_quant)
    os.remove(path_temp_circular_simu)
    return combined_fa_simu


def make_reads_num_assignment(seqs, depth):
    parse_seq = Bio.SeqIO.parse(seqs, format="fasta")
    seq_id, seq_length = zip(*([(fa.id.strip(), len(fa.seq)) for fa in parse_seq]))

    sum_all_this_simulation = round((depth * sum(seq_length)) / (READ_LENGTH + 0.0))
    num_rate_each_loci = [l * random.random() for l in seq_length]
    sum_all_num_rate = sum(num_rate_each_loci)
    fraction_each_loci = [l / (sum_all_num_rate + 0.0) for l in num_rate_each_loci]
    final_num_each_loci = [int(round(f * sum_all_this_simulation)) for f in fraction_each_loci]
    return info_simu_setting(Name=seq_id, Length=seq_length, NumRead=final_num_each_loci)


def dump_setting(info_simu, target_file):
    names = info_simu.Name
    lengths = info_simu.Length
    num_reads = info_simu.NumRead

    with open(target_file, "w") as output_f:
        output_f.write("Name\tLength\tNumRead\n")
        for line in zip(names, lengths, num_reads):
            name, length, num = line
            output_f.write("{name}\t{len}\t{num}\n".format(name=name,
                                                           len=length,
                                                           num=num))


def __get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--circular_gtf", help="annotation for circular annotation in GTF format ", default="")
    parser.add_argument("-l", "--linear_gtf", help="annotation of linear rna in .gtf format", default="")
    parser.add_argument("-g", "--genomic_seq", help="genomic sequence in FASTA format", default="")
    parser.add_argument("-t", "--target_dir", help="where to put output results", default="")
    parser.add_argument("-d", "--depth", help="sequencing depth , default is 30", default=30)

    return parser


if __name__ == "__main__":
    args_parser = __get_parser()
    if len(sys.argv) < 2:
        print(__doc__)
    else:
        args = args_parser.parse_args()
        if args.target_dir and not os.path.exists(args.target_dir):
            os.mkdir(args.target_dir)

        if os.path.exists(args.linear_gtf) and os.path.exists(args.circular_gtf) and os.path.exists(args.genomic_seq):
            circ_as_gtf = os.path.join(args.target_dir, "as.gtf")
            prepare_annotation_as_event(args.circular_gtf, circ_as_gtf)
            fa_simu = prepare_seq(linear_gtf=args.linear_gtf,
                                  circular_gtf=circ_as_gtf,
                                  genomic_seqs=args.genomic_seq,
                                  target_folder=args.target_dir)
            path_of_simu_setting = os.path.join(args.target_dir, "simu.setting")
            info_this_simu_setting = make_reads_num_assignment(fa_simu, int(args.depth))
            dump_setting(info_this_simu_setting, path_of_simu_setting)

        elif not args.circular_gtf and os.path.exists(args.genomic_seq) and os.path.exists(args.linear_gtf):
            __logger.debug("make a simulation with only linear transcripts")
            tmp_linear_transcriptome = os.path.join(args.target_dir, "linear_transcriptome.fa")
            pysrc.wrapper.gffread.do_extract_classic_message_transcript(args.linear_gtf,
                                                                        path_ref_sequence_file=args.genomic_seq,
                                                                        output=tmp_linear_transcriptome)
            info_linear_simu = make_reads_num_assignment(tmp_linear_transcriptome, int(args.depth))
            dump_setting(info_linear_simu, os.path.join(args.target_dir, "simu.setting"))

        else:
            raise FileExistsError("some path is not existing")
