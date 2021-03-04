# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme: test on ciri-entry
#
__author__ = 'zerodel'

import os
import io
import pysrc.file_format.ciri_entry as ce

print(os.path.abspath(os.curdir))


test_str = """circRNA_ID      chr     circRNA_start   circRNA_end     #junction_reads SM_MS_SMS       #non_junction_reads     junction_reads_ratio    circRNA_typegene_id strand  junction_reads_ID
chr1:891303|892653      chr1    891303  892653  4       2_3_1   71      0.101   exon    ENSG00000188976,        -       ST-E00205:443:H7GLYCCXY:5:1114:32146:37665,ST-E00205:443:H7GLYCCXY:5:1106:7222:42587,ST-E00205:443:H7GLYCCXY:5:1121:30330:36698,ST-E00205:443:H7GLYCCXY:5:1121:30350:37436,
chr1:1584361|1647917    chr1    1584361 1647917 3       2_2_1   657     0.009   intergenic_region       n/a     -       ST-E00205:443:H7GLYCCXY:5:1115:3082:72156,ST-E00205:443:H7GLYCCXY:5:1103:8937:18907,ST-E00205:443:H7GLYCCXY:5:2212:10896:31441,
chr1:1647590|1647917    chr1    1647590 1647917 4       4_4_0   740     0.011   intron  ENSG00000008128,ENSG00000268575,        -       ST-E00205:443:H7GLYCCXY:5:1216:29264:47544,ST-E00205:443:H7GLYCCXY:5:1109:27925:73387,ST-E00205:443:H7GLYCCXY:5:2124:13869:32689,ST-E00205:443:H7GLYCCXY:5:2205:2909:19961,
chr1:1665886|1667437    chr1    1665886 1667437 2       2_2_0   98      0.039   intron  ENSG00000268575,ENSG00000227775,ENSG00000215790,        -  ST-E00205:443:H7GLYCCXY:5:2111:7872:54787,ST-E00205:443:H7GLYCCXY:5:2223:19207:66162,
"""



#
# obj_str = io.StringIO(test_str)
#
# lines = obj_str.readlines()[1:]
#
# objs = [ce.CIRIEntry(x) for x in lines]
#
# # =====
# import pandas as pd
# x = pd.read_table("../test/ciri.example")


with open("../test/ciri.example") as ex:
    lines = ex.readlines()[1:]
    objs = [ce.CIRIEntry(x) for x in lines]
    for x in objs:
        print(x.obj.gene_id)


