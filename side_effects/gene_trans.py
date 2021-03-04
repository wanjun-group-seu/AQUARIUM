# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import sys
import os

upper_root = os.path.abspath(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

sys.path.append(upper_root)

from pysrc.file_format.gtf import GTFItem

__doc__ = '''
'''


def extract_gene_name_and_gene_id_mapping(gtf_in, table_out):
    mapping = load_mapping_info_from_gtf(gtf_in)
    with open(table_out, "w") as csv_it:
        csv_it.write("transcript\tgene\n")
        for pair in mapping:
            csv_it.write("{}\n".format(pair))


def load_mapping_info_from_gtf(gtf_in):
    mapping = set()
    with open(gtf_in) as gtf:
        for line in gtf:
            if line.startswith("#"):
                continue
            else:
                some_entry = GTFItem(line.strip())
                attr_this = some_entry.get_attribute()
                if "transcript_id" in attr_this and "gene_id" in attr_this:
                    gene_id = attr_this.get("gene_id", "")
                    transcript_id = attr_this.get("transcript_id", "")
                    gene_name = attr_this.get("gene_name", "")
                    mapping.add("{transcript_id}\t{gene_id}\t{gene_name}".format(gene_id=gene_id,
                                                                                 transcript_id=transcript_id,
                                                                                 gene_name=gene_name))

    return mapping


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print(__doc__)
    else:
        extract_gene_name_and_gene_id_mapping(sys.argv[-2], sys.argv[-1])
