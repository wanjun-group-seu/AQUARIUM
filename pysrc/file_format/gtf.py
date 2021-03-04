# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
GENE_BIOTYPE = "gene_biotype"

__doc__ = '''this module fits the GFF2 format , not GFF3
'''
__author__ = 'zerodel'


class GTFErr(Exception):
    pass


class GTFItemErr(GTFErr):
    """base class of exception of seq-file item utility
    """
    pass


class AttributionIncomplete(GTFItemErr):
    """
    this happens when some attribution segment are lost in gtf file entry.
    """
    pass


class NoGeneIDInGTF(AttributionIncomplete):
    pass


class NoTranscriptIDInGTF(AttributionIncomplete):
    pass


class NoAttributeStringInGTF(GTFItemErr):
    pass


class GTFItem(object):
    """
    read single line of .gtf file , and construct items ,
    """
    sampleAttribute = r'gene_id transcript_id exon_number gene_biotype gene_name p_id protein_id transcript_name tss_id'

    is_ensemble = True

    def __init__(self, line_in_gtf=""):
        """
        Constructor
        """
        if line_in_gtf:
            self._parse_line(line_in_gtf)
        else:  # the "null" condition
            self.__init_elements()
            self.init_null_attribute()

    def __init_elements(self):
        self._seqname = ""
        self._source = ""
        self._feature = ''
        self._start = -1
        self._end = -1
        self._score = '.'
        self._strand = "+"
        self._frame = "."

    def seqname(self):
        return self._seqname

    def set_seqname(self, seqname):
        self._seqname = seqname

    def starts(self):
        return self._start

    def ends(self):
        return self._end

    def set_start(self, start):
        self._start = start

    def set_end(self, end):
        self._end = end

    def get_gene_id(self):
        return self._attributes["gene_id"]

    def get_transcript_id(self):
        return self._attributes["transcript_id"]

    def set_gene_id(self, new_id):
        self._attributes["gene_id"] = new_id

    def set_transcript_id(self, new_transcript_id):
        self._attributes["transcript_id"] = new_transcript_id

    def get_strand(self):
        return self._strand

    def get_attribute(self):
        return self._attributes

    def set_source(self, source):
        self._source = source

    def get_source(self):
        return self._source

    def set_feature(self, feature):
        self._feature = feature

    def get_feature(self):
        return self._feature

    def set_strand(self, strand):
        self._strand = strand

    def set_frame(self, frame):
        self._frame = frame

    def get_frame(self):
        return self._frame

    def get_biotype(self):
        return self._attributes.get(GENE_BIOTYPE, "")

    # def _parse_line(self, line_in_gtf):
    #     """ parse a line in seq-file file ,
    #     only gene id and transcript id will be extracted from attribute string
    #     """
    #     if line_in_gtf.strip().startswith("#"):
    #         self.init_null_attribute()
    #         self._attributes = None
    #         raise AttributionIncomplete("This line is a comment: %s " % line_in_gtf)
    #     else:
    #         gtf_line_parts = line_in_gtf.strip().split("\t")
    #
    #         try:
    #
    #             self._seqname, self._source, self._feature, self._start, self._end, self._score, self._strand, self._frame = gtf_line_parts[:8]
    #             self._start = int(self._start)
    #             self._end = int(self._end)
    #
    #         except IndexError as e:
    #             raise e
    #
    #         try:
    #             #self._check_attribute_string(gtf_line_parts[-1])
    #             self._attributes = self.attribute2dict(gtf_line_parts[-1])
    #         except Exception as e:
    #             raise e

    def _parse_line(self, gtf_line):
        gtf_line_parts = gtf_line.strip().split()  # here use the default parameters

        if gtf_line.strip().startswith("#") or len(gtf_line_parts) < 9:
            self.init_null_attribute()
            self._attributes = None
            raise AttributionIncomplete("This line is a comment: %s " % gtf_line)
        else:
            # think about the first 8 elements
            self._seqname, self._source, self._feature, self._start, self._end, self._score, self._strand, self._frame = gtf_line_parts[:8]
            self._start = int(self._start)
            self._end = int(self._end)

            # try to handle the final attributions part .
            try:
                self._attributes = self.attribute2dict("\t".join(gtf_line_parts[8:]))
            except Exception as e:
                raise e

    @staticmethod
    def _check_attribute_string(gtf_attribute):
        #  todo: this is only for GFF2

        # gff2 has
        #         CDS
        # exon
        # start_codon
        # stop_codon
        # gff 3 CDS
        # exon
        # five_prime_utr
        # gene
        # Selenocysteine
        # start_codon
        # stop_codon
        # three_prime_utr
        # transcript

        if not gtf_attribute:
            # if nothing in attribute string
            raise NoAttributeStringInGTF

        if "gene_id" not in gtf_attribute:
            raise NoGeneIDInGTF

        if "transcript_id" not in gtf_attribute:
            raise NoTranscriptIDInGTF

    @staticmethod
    def attribute2dict(gtf_attribute):
        """
        extract information from the attribute string of gtf file.

        :param gtf_attribute:
        :return:
        """
        return dict([(item.split()[0], item.split()[-1].strip('"'))
                     for item in gtf_attribute.strip().split(";") if item])

    def __eq__(self, other_gtf_item):
        return self._seqname == other_gtf_item.seqname() \
               and self._start == other_gtf_item.starts() \
               and self._end == other_gtf_item.ends() \
               and self._strand == other_gtf_item.get_strand() \
               and self.get_gene_id() == other_gtf_item.get_gene_id() \
               and self.get_transcript_id() == other_gtf_item.get_transcript_id()

    def __len__(self):
        return self._end - self._start + 1

    def init_null_attribute(self):
        """
        two mandatory attributes : gene_id and transcript_id
        :return:
        """
        self._attributes = dict()
        self._attributes.setdefault("gene_id", "")
        self._attributes.setdefault("transcript_id", "")

    def __str__(self):
        """
        the first eight element are following GFF format, and with a description of GTF
        :return:
        """
        return "\t".join([self._seqname,
                          self._source,
                          self._feature,
                          str(self._start),
                          str(self._end),
                          self._score,
                          self._strand,
                          self._frame,
                          self._attr2str()])

    def whether_ensemble(self, yes_or_no):
        self.is_ensemble = yes_or_no

    def _attr2str(self):
        attr_rebuild = "; ".join(['%s "%s"' % (key, self._attributes.get(key))
                                  for key in self.sampleAttribute.strip().split()
                                  if key in self._attributes.keys()])

        if self.is_ensemble:
            return "%s;" % attr_rebuild

        return attr_rebuild


def filter_gff_by_source(gff, bio_type):
    res = []
    with open(gff) as obj_gff:
        for line in obj_gff:
            try:
                this_entry = GTFItem(line.strip())
            except AttributionIncomplete:
                continue
            else:
                # non-stable bio type name
                if this_entry.get_biotype() in set(bio_type):
                    res.append(line)
    return res


gene_biotype_vals = """antisense
lincRNA
miRNA
misc_RNA
Mt_rRNA
Mt_tRNA
processed_pseudogene
processed_transcript
protein_coding
pseudogene
ribozyme
rRNA
scaRNA
sense_intronic
snoRNA
snRNA
sRNA
TEC
transcribed_processed_pseudogene
transcribed_unprocessed_pseudogene
unprocessed_pseudogene
"""
