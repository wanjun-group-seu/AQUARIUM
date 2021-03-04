# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import pysrc.body.cli_opts
import pysrc.sub_module.qc

__doc__ = '''
'''
__author__ = 'zerodel'

if __name__ == "__main__":
    par = pysrc.body.cli_opts.simple_arg_parser_only_config_option().parse_args()
    pysrc.sub_module.qc.main(par.config)
