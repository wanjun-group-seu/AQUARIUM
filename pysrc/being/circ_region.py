# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import pysrc.body.logger

__doc__ = """ circ_region represent a circRNA as a region in genomic coordinate
"""

__author__ = 'zerodel'


_logger = pysrc.body.logger.default_logger("CIRC_REGION")


class RegionCirc(object):
    def __init__(self, chrom, p_start, p_end, name, count=0, host="", id_str=""):
        if id_str:
            try:
                self.parse_id_str(id_str)
            except ValueError as e:
                raise e
        else:
            self.chrom = chrom
            self.index_start = str(p_start).strip()
            self.index_end = str(p_end).strip()
            self.name = name
            self.score = str(count)
            self.host = host

    def as_bed(self):
        return "{chrom}\t{pos_start}\t{pos_end}\t{name}\t{score}".format(chrom=self.chrom,
                                                                         pos_start=self.index_start,
                                                                         pos_end=self.index_end,
                                                                         name=self.name,
                                                                         score=self.score
                                                                         )

    def parse_id_str(self, id_str):
        if "@" in id_str:
            id_segments = id_str.strip().split("@")
            region_str = id_segments[0]
            host_info = id_segments[1:]
        elif "|" in id_str and ":" in id_str:
            region_str = id_str
            host_info = ""

        else:
            _logger.error("id_str should follow such pattern chr:star|end@host_gene")
            raise ValueError("Wrong id string : {}".format(id_str))

        if not region_str:
            _logger.error("meet en empty string here: {}".format(id_str))
            raise ValueError("Wrong id string : {}".format(id_str))
        else:
            chrom, region = region_str.strip().split(":")
            pos_start, pos_end = region.strip().split("|")

            self.chrom = chrom
            self.index_start = pos_start.strip()
            self.index_end = pos_end.strip()
            self.name = id_str
            self.score = 0
            self.host = host_info

    def __eq__(self, other):
        return self.chrom == other.chrom and self.index_start == other.pos_start and self.index_end == other.pos_end

    def get_id(self):
        return "{chrom}:{start}|{end}".format(chrom=self.chrom, start=self.index_start, end=self.index_end)


if __name__ == "__main__":
    print(__doc__)
