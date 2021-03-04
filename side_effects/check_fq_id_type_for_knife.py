# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import os

import sys

base_path = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]
sys.path.append(base_path)

import pysrc.wrapper.knife

__doc__ = '''
'''

__author__ = 'zerodel'

pysrc.wrapper.knife.report_id_type(sys.argv[-2], sys.argv[-1])
