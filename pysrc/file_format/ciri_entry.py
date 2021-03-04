# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import os
import collections

import pysrc.body.logger

JUNCTION_READS_LIMIT = 5

__doc__ = '''
'''
__author__ = 'zerodel'

_logger = pysrc.body.logger.default_logger("CIRI_ENTRY")

HEADER_v1 = """circRNA_ID	chr	circRNA_start	circRNA_end	#junction_reads	SM_MS_SMS	#non_junction_reads	junction_reads_ratio	circRNA_type	gene_id	junction_reads_ID"""

HEADER_v2 = """circRNA_ID	chr	circRNA_start	circRNA_end	#junction_reads	SM_MS_SMS	#non_junction_reads	junction_reads_ratio	circRNA_type	gene_id	strand	junction_reads_ID
"""

slots_v1 = [x.strip("#") for x in HEADER_v1.split()]


# CIRI1 = collections.namedtuple("CIRI1", slots_v1)

class CIRI1(object):
    def __init__(self, dict_ciri1):
        self.__dict__.update(dict_ciri1)


slots_v2 = [x.strip("#") for x in HEADER_v2.split()]


# CIRI2 = collections.namedtuple("CIRI2", slots_v2)
class CIRI2(object):
    def __init__(self, dict_input):
        self.__dict__.update(dict_input)


class CIRIEntry(object):
    def __init__(self, str_line):
        """ construct an empty ciri entry or from a string.
        currently , we assume the line only follows HEADER_v1 or HEADER_v2 
        
        """
        if str_line:

            line_parts = str_line.strip().split("\t")

            if 1 == len(line_parts):
                _logger.warning("this ciri line is not separated by tab : {}".format(str_line))
                line_parts = str_line.strip().split()
                _logger.debug("now it has {} parts ".format(len(line_parts)))

            if len(line_parts) == len(slots_v1):  # v1
                # self.obj = CIRI1(**dict(zip(slots_v1, line_parts)))
                self.obj = CIRI1(dict(zip(slots_v1, line_parts)))
                self.obj.strand = None

            elif len(line_parts) == len(slots_v2):  # v2
                # self.obj = CIRI2(**dict(zip(slots_v2, line_parts)))
                self.obj = CIRI2(dict(zip(slots_v2, line_parts)))
            else:
                _logger.error("this line has {} parts".format(len(line_parts)))
                raise ValueError("Error:wrong CIRI format in : {}".format(str_line))
        else:
            raise ValueError("Error: empty line ")

        # error prone:  can not set namedtuple attributes
        if self.obj.gene_id.strip().startswith("inter"):
            self.obj.gene_id = "n/a"

        self.obj.gene_id = self.obj.gene_id.strip(",")  # try to fix the tailed coma .

        if "," in self.obj.gene_id:
            _logger.debug("{circ} has multiple host : {host}, pick one ".format(circ=self.obj.circRNA_ID,
                                                                                host=self.obj.gene_id))
            # here to pick the first gene id as the only gene id
            self.obj.gene_id = self.obj.gene_id.split(",")[0].strip(",")

        self.id_show_host = "%s@%s" % (self.obj.circRNA_ID, self.obj.gene_id)

    def filter_by_junction_reads(self, num_reads_lower_limit=JUNCTION_READS_LIMIT):
        return int(self.obj.junction_reads) >= num_reads_lower_limit

    def to_bed_string(self):
        """transfer this object into a .bed file string
        """
        entry = self.obj
        strand = "0"
        if isinstance(entry, CIRI2):
            strand = entry.strand

        return "\t".join(
            [entry.chr, entry.circRNA_start, entry.circRNA_end, self.id_show_host, entry.junction_reads,
             strand]).strip()


def transform_ciri_to_bed(ciri_output_file, num_read_lower_limit=JUNCTION_READS_LIMIT, output_bed_file=""):
    abs_ciri_dir = os.path.abspath(ciri_output_file)
    main_part_ciri_path = os.path.splitext(abs_ciri_dir)[0]

    if isinstance(num_read_lower_limit, str):
        _logger.warning("doing a explict type transforming here:num_read_lower_limit")
        num_read_lower_limit = int(num_read_lower_limit)

    if not output_bed_file:
        output_bed_file = ".".join([main_part_ciri_path, "bed"])

    with open(ciri_output_file) as ciri_file:
        ciri_file.readline()  # file head should be skipped

        with open(output_bed_file, "w") as exporter:

            for line in ciri_file:
                ciri_line_entry = CIRIEntry(line.strip())

                if ciri_line_entry.filter_by_junction_reads(num_read_lower_limit):
                    new_bed_line = ciri_line_entry.to_bed_string()
                    exporter.write(new_bed_line + "\n")
                else:
                    _logger.warning("encounter a dis-qualified entry at %s" % str(ciri_line_entry))


CIRC_REPORT_TABLE_HEADER_FIRST = "circRNA_ID"


def is_ciri_file_intact(ciri_file):
    if os.path.exists(ciri_file):
        last_line_ciri_report = os.popen("tail -n 1 {}".format(ciri_file)).read().strip()
        if len(last_line_ciri_report) > 5:
            head_ciri_line = last_line_ciri_report.split()[0].strip()
            tail_ciri_line = last_line_ciri_report[-1]

            if not head_ciri_line == CIRC_REPORT_TABLE_HEADER_FIRST and tail_ciri_line == ",":
                return True

    return False