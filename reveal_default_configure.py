# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#


import os.path
import sys

import pysrc.body.config

sys.path.append(os.path.dirname(__file__))

__doc__ = ''' print out where the default config file is
'''
__author__ = 'zerodel'


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print(pysrc.body.config.throw_out_where_the_default_config_is())
    else:
        print(__doc__)
