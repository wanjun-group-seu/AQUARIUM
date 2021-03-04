# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import sys
import os

upper_root = os.path.abspath(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

sys.path.append(upper_root)

import pysrc.file_format.gtf

__doc__ = '''
'''


def extract_gene_name_and_gene_id_mapping(gtf_in, table_out):
    mapping = load_mapping_info_from_gtf(gtf_in)
    with open(table_out, "w") as csv_it:
        csv_it.write("name\tid\n")
        for pair in mapping:
            csv_it.write("{}\n".format(pair))


def load_mapping_info_from_gtf(gtf_in):
    mapping = set()
    with open(gtf_in) as gtf:
        for line in gtf:
            try:
                some_entry = pysrc.file_format.gtf.GTFItem(line.strip())
            except pysrc.file_format.gtf.AttributionIncomplete:
                continue
            attr_this = some_entry.get_attribute()
            if attr_this and isinstance(attr_this, dict) and "gene_name" in attr_this and "gene_id" in attr_this:
                mapping.add("{name}\t{id}".format(name=attr_this.get("gene_name"),
                                                  id=attr_this.get("gene_id")))

    return mapping


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print(__doc__)
    else:
        extract_gene_name_and_gene_id_mapping(sys.argv[-2], sys.argv[-1])
