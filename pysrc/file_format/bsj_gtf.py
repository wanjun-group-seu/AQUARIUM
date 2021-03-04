# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import itertools
import multiprocessing
import os

import gffutils
import pysrc.file_format.ciri_entry
import pysrc.file_format.gtf
import pysrc.body.logger
from pysrc.body.utilities import guess_num_core

__doc__ = '''
'''
__author__ = 'zerodel'

OVERLAP_WINDOW_WIDTH = 3
valid_feature_type_circ = ["processed_transcript", "protein_coding"]

_logger = pysrc.body.logger.default_logger("BSJ")


class PredictedCircularRegion(object):
    def __init__(self, args_tuple, **kwargs):
        if args_tuple:
            self.predict_id, self.seqid, self.start, self.end, self.strand = args_tuple
            self.start = int(self.start)
            self.end = int(self.end)
        elif kwargs:
            self.seqid = kwargs.get("chromosome")
            self.start = int(kwargs.get("start"))
            self.end = int(kwargs.get("end"))
            self.predict_id = kwargs.get("given_id")
            self.strand = kwargs.get("strand")
        else:
            raise NameError("Error@PredictedCircularRegion: wrong arguments given")

    def is_flanking(self, gff_feature):
        return self.start <= int(gff_feature.start) and self.end >= int(gff_feature.end)

    def extract_flanked_linear_entries(self, gffutils_database):
        # extract all isoform part from some database of gtf file
        # here we assume all exon has attribution of 'transcript_id'
        transcript_exon_dict = {}
        for linear_isoform in gffutils_database.region(seqid=self.seqid, start=self.start, end=self.end,
                                                       strand=self.strand,
                                                       featuretype="transcript"):
            corresponding_circular_exons = [exon for exon in
                                            gffutils_database.children(linear_isoform.id, featuretype="exon",
                                                                       order_by="start",
                                                                       limit=(self.seqid, self.start, self.end),
                                                                       completely_within=True)]
            if corresponding_circular_exons:
                transcript_exon_dict.setdefault(linear_isoform.id, corresponding_circular_exons)

        return transcript_exon_dict

    @staticmethod
    def generate_exon_for_circular_isoform(host_seqname, start, end, host_gene_id, host_tran_id, strand="+", frame="."):
        artificial_exon = pysrc.file_format.gtf.GTFItem()
        artificial_exon.set_start(int(start))
        artificial_exon.set_end(int(end))
        artificial_exon.set_gene_id(host_gene_id)
        artificial_exon.set_transcript_id(host_tran_id)
        artificial_exon.set_seqname(host_seqname)
        artificial_exon.set_source("ciri")
        artificial_exon.set_feature("exon")
        artificial_exon.set_strand(strand)
        artificial_exon.set_frame(frame)
        return artificial_exon

    @staticmethod
    def guess_feature_type(feature):
        if feature.source in valid_feature_type_circ:
            return feature.source
        else:
            return feature.attributes.get("gene_biotype", ["protein_coding"])[0]

    def arrange_exons_the_naive_way(self, db):
        exons_raw = list(set(
            [(exon.seqid, self.guess_feature_type(exon), exon.start, exon.end, exon.strand, exon.frame)
             for exon in db.region(seqid=self.seqid,
                                   start=int(self.start),
                                   end=int(self.end),
                                   strand=self.strand,
                                   featuretype="exon")]))

        exon_filtered = []  # start filter exon objects
        for exon in exons_raw:
            exon_seqid, exon_source, exon_start, exon_end, exon_strand, exon_frame = exon

            if exon_seqid == 'chrM':
                continue

            if exon_source not in valid_feature_type_circ:
                continue

            if exon_start < self.start - OVERLAP_WINDOW_WIDTH:
                exon_start = self.start

            if exon_end > self.end + OVERLAP_WINDOW_WIDTH:
                exon_end = self.end

            exon_filtered.append((exon_seqid, exon_source, exon_start, exon_end, exon_strand.strip(), exon_frame))

        exon_filtered = sorted(exon_filtered, key=lambda exon_str: exon_str[2])

        artificial_exons = []

        transcript_id = self.predict_id
        gene_id = "n/a"

        if len(self.predict_id.split("@")) == 2:
            transcript_id, gene_id = self.predict_id.split("@")

        for exon_locus in exon_filtered:
            exon_seqid, exon_source, exon_start, exon_end, exon_strand, exon_frame = exon_locus

            # [2019_10_15 17:34], add .r only when the region strand is unclear and exon are on reverse strand
            is_this_region_stand_clear = self.strand in {"+", "-"}
            transcript_on_reverse_strand = exon_strand == "-" and not is_this_region_stand_clear
            transcript_id_show_strand = "%s.r" % transcript_id.strip() if transcript_on_reverse_strand else transcript_id.strip()

            neo_isoform = self.generate_exon_for_circular_isoform(host_seqname=exon_seqid, start=exon_start,
                                                                  end=exon_end,
                                                                  host_gene_id=gene_id,
                                                                  host_tran_id=transcript_id_show_strand,
                                                                  strand=exon_strand, frame=exon_frame)
            artificial_exons.append(neo_isoform
                                    )
        return artificial_exons

    def mark_extracted_exons(self, dict_transcript_exon):
        # this function is after the extract_flanked.... function
        marked_exons = []
        for transcript_id in dict_transcript_exon.keys():
            for exon in dict_transcript_exon[transcript_id]:
                neo_exon = self.simplify_this_feature(exon, new_source="circRNA",
                                                      new_transcript_id="%s@%s" % (self.predict_id, transcript_id))
                marked_exons.append(neo_exon)

        return marked_exons

    @staticmethod
    def simplify_this_feature(feature_from_gffutils_db, new_source="", new_transcript_id=""):
        artificial_exon = pysrc.file_format.gtf.GTFItem(str(feature_from_gffutils_db))
        formal_gene_id = artificial_exon.get_gene_id()
        formal_trans_id = artificial_exon.get_transcript_id()
        artificial_exon.init_null_attribute()
        artificial_exon.set_gene_id(formal_gene_id)
        artificial_exon.set_transcript_id(formal_trans_id)

        if new_source:
            artificial_exon.set_source(new_source)
        if new_transcript_id:
            artificial_exon.set_transcript_id(new_transcript_id)

        return artificial_exon


def parse_bed_line(line):
    parts = line.strip().split()
    if len(parts) < 4:
        raise KeyError("Error: not right bed file type, not enough columns ")
    chr_name, start, end, isoform_id = parts[:4]

    strand = parts[5] if len(parts) > 5 else None
    return isoform_id, chr_name, start, end, strand


def parse_ciri_line(line):
    entry_ciri = pysrc.file_format.ciri_entry.CIRIEntry(line)
    return entry_ciri.id_show_host, entry_ciri.obj.chr, entry_ciri.obj.circRNA_start, entry_ciri.obj.circRNA_end, \
           entry_ciri.obj.strand


def parse_ciri_as_region(ciri_output, comment_char="#"):
    with open(ciri_output) as ciri_reader:
        ciri_reader.readline()
        for line in ciri_reader:
            if not line.strip().startswith(comment_char):
                yield PredictedCircularRegion(parse_ciri_line(line))


def parse_bed_as_region(bed_output_no_header, comment_char="#"):
    with open(bed_output_no_header) as read_bed:
        for line in read_bed:
            if not line.strip().startswith(comment_char):
                yield PredictedCircularRegion(parse_bed_line(line))


def get_gff_database(gtf_file):
    path_main, file_part = os.path.split(gtf_file)
    file_body_name, file_suffix = os.path.splitext(file_part)

    if ".gtf" == file_suffix:
        db_file_path = os.path.join(path_main, ".".join([file_body_name, "db"]))
        if os.path.exists(db_file_path):
            db = gffutils.FeatureDB(db_file_path)
        else:
            # @WARNING: different version of GTF will make this process time exhausting
            _logger.warning("Please check GTF file version. it is strongly recommended to build GTF database manually")
            db = gffutils.create_db(gtf_file, db_file_path)

    elif ".db" == file_suffix:
        _logger.debug("using existing database:{}".format(gtf_file))
        db = gffutils.FeatureDB(gtf_file)
    else:
        raise FileNotFoundError("Can not Get the right gffutils database file")
    return db


class SimpleMapReduce(object):
    def __init__(self, map_func, reduce_func, num_cores=None):
        self.map_func = map_func
        self.reduce_func = reduce_func
        self.pool = multiprocessing.Pool(num_cores)

    def __call__(self, inputs, chunk_size=1):
        map_response = self.pool.map(self.map_func, inputs)

        intermediate_generator = itertools.chain(*map_response)
        reduced_values = self.reduce_func(intermediate_generator)
        return reduced_values


def _get_exons(intervals_and_gff_path):
    regs, gff_path = intervals_and_gff_path
    db = get_gff_database(gff_path)
    res = []
    for reg in regs:
        res.extend(reg.arrange_exons_the_naive_way(db))
    return res


def _get_exons_with_isoforms(intervals_and_gff_path):
    regs, gff_path = intervals_and_gff_path
    gff_db = get_gff_database(gff_path)
    res = []
    for reg in regs:
        flanked_linear_entries = reg.extract_flanked_linear_entries(gff_db)
        exons_marked_circular = reg.mark_extracted_exons(flanked_linear_entries)
        res.extend(exons_marked_circular)
    return res


def _all_exons_as_list(exons):
    return list(exons)


def get_gtf_mp(path_bed, gff_db, is_structure_show=False, num_process=0):
    regions = load_region_from_file(path_bed)

    region_groups_for_each_cpu = equal_divide(regions, num_process)

    gff_db_each_cpu = [gff_db] * num_process if num_process else [gff_db]

    if is_structure_show:
        a = SimpleMapReduce(_get_exons_with_isoforms, _all_exons_as_list, num_process)
    else:
        a = SimpleMapReduce(_get_exons, _all_exons_as_list, num_process)

    result = a(zip(region_groups_for_each_cpu, gff_db_each_cpu))
    return result


def load_region_from_file(path_bed):
    func_to_get_region = parse_bed_as_region if path_bed.endswith(".bed") else parse_ciri_as_region
    regions = list(func_to_get_region(path_bed))
    return regions


def equal_divide(regs, num_process):
    if isinstance(num_process, int):

        if num_process > len(regs):
            _logger.warning("too less regions , even less than cpu core numbers")
            return [[x] for x in regs]
        else:
            size_chunk = int(round(len(regs) / num_process))
            if size_chunk * num_process < len(regs):
                size_chunk += 1

            _logger.debug("number of regions each cpu core: {}".format(size_chunk))
            _logger.debug("number of cpu core: {}".format(num_process))

            reg_parts = [regs[x: x + size_chunk] for x in
                         range(0, len(regs), size_chunk)]

            return reg_parts
    else:
        raise TypeError("Error@gtf_processing: process number should be a integer, but got : {}".format(num_process))


def do_make_gtf_for_circular_prediction_greedy(circular_candidate_regions, gff_db, output_gtf_path_name="",
                                               is_isoform_structure_shown=False):
    whole_exons = intersect_region_genome_annotation(circular_candidate_regions, gff_db, is_isoform_structure_shown)

    exons_to_gtf_file(whole_exons, output_gtf_path_name)


def exons_to_gtf_file(whole_exons, output_gtf_path_name):
    with open(output_gtf_path_name, "w") as out_gtf:
        out_gtf.write("\n".join([str(exon) for exon in whole_exons]))
        _logger.debug("GTF file already in disk : {}".format(output_gtf_path_name))


def intersect_region_genome_annotation(circular_candidate_regions, genomic_annotation,
                                       is_isoform_structure_shown=False):
    num_core = guess_num_core()
    _logger.debug("Number of CPU use during GTF generating : {}".format(num_core))
    whole_exons = get_gtf_mp(path_bed=circular_candidate_regions, gff_db=genomic_annotation,
                             is_structure_show=is_isoform_structure_shown, num_process=num_core)
    _logger.debug("number of exon in memory : {}".format(len(whole_exons)))
    return whole_exons
