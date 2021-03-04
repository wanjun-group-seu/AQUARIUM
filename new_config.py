# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#


import os.path
import sys

import pysrc.body.config

sys.path.append(os.path.dirname(__file__))

__doc__ = ''' this command will give you a copy of default configuration file at given path,
you can change it with your editor
'''
__author__ = 'zerodel'

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print(__doc__)
    else:

        path_to_config_default = pysrc.body.config.throw_out_where_the_default_config_is()
        if not os.path.exists(path_to_config_default):
            raise FileNotFoundError("Error: unable to find default config file at {}".format(path_to_config_default))

        import shutil

        shutil.copy2(path_to_config_default, sys.argv[-1])
