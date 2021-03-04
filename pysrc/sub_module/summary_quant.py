# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#


import pysrc.body.logger
import pysrc.file_format.ciri_entry
import pysrc.file_format.gtf

CIRCULAR = "circular"

LINC = "linc"

MRNA = "mrna"

__doc__ = """
"""
__author__ = 'zerodel'

__COMMENT_CHAR = "#"

_logger = pysrc.body.logger.default_logger("SUMMARIZE_QUANTIFICATION_OUTPUT")

STRING_NA = 'n/a'


class TranscriptOwnership(object):
    """ object contains attribution of different transcript isoforms 
    members:
        self.gene_of : describes host gene of each isoform
        self.origin_of :  describes origin of each isoform
        self.transcripts_of : describes transcripts owned by each gene
        
    methods:
        merge(other_object): combine another object in the same class
        parse_gtf():parser a gtf file
        parse_ciri(): parse a ciri output 
        to_text_table_lines(): export transcript information to a text file
    """

    def __init__(self, dd_gene_of=None, dd_type_of=None, dd_transcripts_of=None):
        self.gene_of = dd_gene_of if dd_gene_of else {}
        self.type_of = dd_type_of if dd_type_of else {}
        self.transcripts_of = dd_transcripts_of if dd_transcripts_of else {}
        self.default_output_format = "{iso}\t{gene}\t{origin}\n"

    def merge(self, other):
        """use information in `other` to update 
        :param other: another object of TranscriptOwnership class
        """
        if not isinstance(other, TranscriptOwnership):
            _logger.error("try to merge a different type object")
            raise TypeError("can not merge object of different class")
        self.gene_of.update(other.gene_of)
        self.type_of.update(other.type_of)
        self.transcripts_of.update(other.transcripts_of)

    def parse_gtf(self, gtf_file):
        """parse a given annotation file in  .gtf format  
        :param gtf_file: path to a gtf file 
        """
        info_in_gtf = _load_or_update_transcript_ownership_relation(gtf_file,
                                                                    _parse_single_line_as_gtf,
                                                                    None)
        self.merge(info_in_gtf)

    def parse_ciri(self, ciri_file):
        """parse a given CIRI output file, current support CIRIv1.2 
        :param ciri_file: path to a CIRI output file
        """
        info_in_ciri = _load_or_update_transcript_ownership_relation(ciri_file,
                                                                     _parse_single_line_ciri,
                                                                     None)
        self.merge(info_in_ciri)

    def to_text_table_lines(self, include_header_line=True, format_str="", ):
        """export the host gene information and original type of each isoform.
        :param include_header_line: logical , output will contain column name if True 
        :param format_str: a string follows the str.format definition,  "iso", "gene", "origin" is allowed. 
        """
        if not format_str:
            format_str = self.default_output_format

        if include_header_line:
            yield format_str.replace("{", "").replace("}", "")

        for iso in self.__meaningful_isoforms_sorted():
            gene = self.gene_of.get(iso, "n/a")
            origin = self.type_of.get(iso, "n/a")
            yield format_str.format(**locals())

    def __meaningful_isoforms_sorted(self):
        query_gene = [iso for iso in self.gene_of]
        query_type = [iso for iso in self.type_of]
        query_gene.extend(query_type)
        all_isoforms = list(sorted(set(query_gene)))
        return all_isoforms

    def summarize_to_gene_level(self, transcript_level_quantification):

        unique_types = sorted(list(set((self.type_of[iso] for iso in self.type_of))))

        _logger.debug("the types of transcript are %s" % str(unique_types))

        all_genes = [x for x in self.transcripts_of]
        if "n/a" not in all_genes:
            all_genes.append("n/a")

        _logger.debug("there are %d genes in all " % len(all_genes))

        def slot_for_each_gene(list_of_gene):
            return dict(zip(list_of_gene, [0] * len(list_of_gene)))

#        column_types = dict(zip(unique_types, [slot_for_each_gene(all_genes)] * len(unique_types)))
        column_types = dict()
        for iso_type in unique_types:
            column_types[iso_type] = slot_for_each_gene(all_genes)

        for iso in transcript_level_quantification:
            type_this_iso = self.type_of.get(iso, "mrna")
            gene_this_iso = self.gene_of.get(iso, "n/a")
            quantification_this_iso = float(transcript_level_quantification.get(iso, 0.0))

            # if "linc" == type_this_iso:
            #     _logger.debug("add a linc : {iso}".format(iso=iso))

            column_types[type_this_iso][gene_this_iso] += quantification_this_iso

        return {"genes": all_genes, "val": column_types, "types": unique_types}


class RawRegion(object):
    def __init__(self, str_junction=""):
        if str_junction:
            self.raw_str = str_junction
            self.chr, self.start, self.end = self.parse_ciri_junction(str_junction.split(".")[0].strip())
        else:
            self.chr, self.start, self.end, self.raw_str = "", None, None, ""

    def __str__(self):
        return self.raw_str

    def is_overlap_with(self, other):
        on_the_same_chrome = self.chr == other.chr
        has_overlap = (self.start - other.end) * (self.end - other.start) < 0
        return on_the_same_chrome and has_overlap

    @staticmethod
    def parse_ciri_junction(junction_str):
        chr_name, junction_info = junction_str.strip().split(":")
        junction_start, junction_end = junction_info.strip().split("|")
        return chr_name, int(junction_start), int(junction_end)


def _make_a_cluster(list_of_junction_obj):
    starts = min([j.start for j in list_of_junction_obj])
    ends = max([j.end for j in list_of_junction_obj])
    chr_name = list_of_junction_obj[0].chr

    region_id = "Region_{name_chromosome}:{pos_start}|{pos_end}".format(name_chromosome=chr_name,
                                                                        pos_start=starts,
                                                                        pos_end=ends)
    junction_ids = [str(obj) for obj in list_of_junction_obj]
    return region_id, junction_ids


def _insert_current_cluster(cluster, dict_transcripts_of_gene):
    cluster_id, junctions_in_cluster = _make_a_cluster(cluster)
    _logger.debug("insert cluster %s : %s" % (cluster_id, ", ".join(junctions_in_cluster)))
    dict_transcripts_of_gene[cluster_id] = junctions_in_cluster


def _sort_junctions_by_chrome(na_junctions):
    chromes = {}
    for j_str in na_junctions:
        j = RawRegion(j_str)
        chromes.setdefault(j.chr, set()).add(j)
    return chromes


def replace_na_with_dis_connected_locus(dict_transcripts_of_gene):
    if STRING_NA in dict_transcripts_of_gene:
        na_junctions = dict_transcripts_of_gene.pop(STRING_NA)
        chromes = _sort_junctions_by_chrome(na_junctions)

        for chr_name in chromes:
            junction_this_chrome = sorted(chromes[chr_name], key=lambda x: x.start)
            cluster = []
            for obj_j in junction_this_chrome:
                if len(cluster) == 0:
                    cluster.append(obj_j)
                else:
                    last_obj = cluster[-1]
                    if last_obj.is_overlap_with(obj_j):
                        cluster.append(obj_j)
                    else:
                        _insert_current_cluster(cluster, dict_transcripts_of_gene)
                        cluster = [obj_j]
            _insert_current_cluster(cluster, dict_transcripts_of_gene)
    return dict_transcripts_of_gene


def _load_or_update_transcript_ownership_relation(file_contains_info, func_to_parse_single_line,
                                                  dict_transcripts_of_gene=None,
                                                  dict_type_of_transcript=None,
                                                  dict_gene_of_transcript=None,
                                                  obj_transcript_ownership=None):
    if not dict_transcripts_of_gene or not isinstance(dict_transcripts_of_gene, dict):
        dict_transcripts_of_gene = {}

    if not dict_type_of_transcript or not isinstance(dict_type_of_transcript, dict):
        dict_type_of_transcript = {}

    if not dict_gene_of_transcript or not isinstance(dict_gene_of_transcript, dict):
        dict_gene_of_transcript = {}

    if obj_transcript_ownership:
        _logger.info("loading given transcript ownership information")

        dict_transcripts_of_gene = obj_transcript_ownership.transcripts_of
        dict_type_of_transcript = obj_transcript_ownership.type_of
        dict_gene_of_transcript = obj_transcript_ownership.gene_of

    _logger.info("updating 'transcripts_under_gene' using %s" % file_contains_info)

    with open(file_contains_info) as annotation_file_obj:
        for single_line in annotation_file_obj:
            if not single_line.strip().startswith(__COMMENT_CHAR):
                gene_entry, transcript_entry, type_of_transcript = func_to_parse_single_line(single_line)
                if gene_entry and transcript_entry and type_of_transcript:
                    dict_type_of_transcript[transcript_entry] = type_of_transcript
                    dict_gene_of_transcript[transcript_entry] = gene_entry
                    dict_transcripts_of_gene.setdefault(gene_entry, set()).add(transcript_entry)
                else:
                    _logger.debug("meet a line not containing all 3 info: %s" % single_line)

    return TranscriptOwnership(dd_gene_of=dict_gene_of_transcript, dd_transcripts_of=dict_transcripts_of_gene,
                               dd_type_of=dict_type_of_transcript)


def _parse_single_line_as_gtf(single_line):
    try:
        entry = pysrc.file_format.gtf.GTFItem(single_line.strip())
        gene_id = entry.get_gene_id()
        transcript_id = entry.get_transcript_id()
    except:
        _logger.debug("in complete gtf line: %s" % single_line)
        return None, None, None
    else:

        if is_this_id_circular(transcript_id):
            transcript_type = CIRCULAR
        elif "lincRNA" == entry.get_biotype():
            # logger.debug("found linc: {iso}".format(iso=transcript_id))
            transcript_type = LINC
        else:
            transcript_type = MRNA

        return gene_id, transcript_id, transcript_type


def get_mapping_info_from_gtf(gtf_file, obj_ownership=None):
    return _load_or_update_transcript_ownership_relation(gtf_file,
                                                         _parse_single_line_as_gtf,
                                                         obj_transcript_ownership=obj_ownership)


def _parse_single_line_ciri(line):
    ciri_this_line = pysrc.file_format.ciri_entry.CIRIEntry(line.strip())
    gene_id = ciri_this_line.obj.gene_id
    transcript_id = ciri_this_line.obj.circRNA_ID
    return gene_id, transcript_id, "circular"


def get_mapping_info_from_ciri(ciri_file, obj_ownership=None):
    return _load_or_update_transcript_ownership_relation(ciri_file,
                                                         _parse_single_line_ciri,
                                                         obj_transcript_ownership=obj_ownership)


def is_this_id_circular(id_this):
    return str(id_this).startswith("chr") or (":" in str(id_this) and "|" in str(id_this))


def load_quantify_report(quantify_report_file, key_name="Name", val_name="TPM"):
    """
    main interface for load sailfish quant.sf report file. 
    
    :param quantify_report_file: path to sailfish quantification file 
    :param key_name: which column indicates the transcriptome entry
    :param val_name: which column contains the data of interest 
    :return: a dict with each line as key:value pairs . 
    """

    from collections import namedtuple

    with open(quantify_report_file) as get_quant:
        header_parts = get_quant.readline().strip().split("\t")
        # default file header "Name	Length	EffectiveLength	TPM	NumReads"
        # we assume that the header line contains valid names
        EntryISO = namedtuple("EntryISO", field_names=header_parts)

        quant_dict = {}

        for quant_line in get_quant:
            line_parts = quant_line.strip().split("\t")
            obj_this_line = EntryISO(**dict(zip(header_parts, line_parts)))
            id_this = getattr(obj_this_line, key_name)
            value_this = getattr(obj_this_line, val_name)

            quant_dict[id_this] = value_this

    return quant_dict


def dump_summarized_result(summarized_dd, output_path):
    """
    write summarize_* function result into text files 
    :param summarized_dd: output of some summarize_* function, should have three members: genes, val, types, 
    :param output_path: path to specify the tab-delimit file. 
    """
    genes = summarized_dd["genes"]
    val = summarized_dd["val"]
    value_types = summarized_dd["types"]

    with open(output_path, "w") as dump_it:
        dump_it.write("gene\t%s\n" % "\t".join(value_types))

        _logger.info(
            "start writing summarized quantification to :{out_file}".format(out_file=output_path))
        for single_gene in genes:
            line_this = [single_gene]
            for val_type in value_types:
                line_this.append(str(val[val_type][single_gene]))

            dump_it.write("%s\n" % "\t".join(line_this))

        _logger.info("finished : {out_file}".format(out_file=output_path))


def aggregate_isoform_quantify_result(quant_sf, summarized_output, gtf_annotation, ciri_output=""):
    """
    simple wrapped-up process to aggregate isoform-level quantification result into gene level
    and show different RNA source. 
    
    :param quant_sf: sailfish quantification report , isoform level 
    :param summarized_output: result tab-delimit file 
    :param gtf_annotation: gene-annotation file in .gtf format 
    :param ciri_output: optional result from CIRI (current for v1.2)
    """
    info_iso = TranscriptOwnership()

    info_iso.parse_gtf(gtf_annotation)
    #    info_iso = get_mapping_info_from_gtf(gtf_annotation, info_iso)

    if ciri_output:
        _logger.info("loading ciri report from : %s" % ciri_output)
        info_iso.parse_ciri(ciri_output)
        #        info_iso = get_mapping_info_from_ciri(ciri_output, info_iso)
        _logger.info("loading CIRI output is ok")

    # info_iso.transcripts_of = replace_na_with_dis_connected_locus(info_iso.transcripts_of)

    quantify_of_transcript_level = load_quantify_report(quant_sf)

    df_2d_res = info_iso.summarize_to_gene_level(quantify_of_transcript_level)

    _logger.debug("dump those data to %s" % summarized_output)

    dump_summarized_result(df_2d_res, summarized_output)
