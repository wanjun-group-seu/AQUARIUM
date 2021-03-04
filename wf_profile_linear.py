# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import argparse

import pysrc.body.config

import pysrc.body.worker
import pysrc.sub_module.align
from pysrc.sub_module import align, quantify

_OPT_KEY_QUANTIFIER = "quantifier"

_OPT_KEY_MAPPER = "mapper"

_QUANTIFIERS_DO_NOT_NEED_EXTERNAL_MAPPER = ["sailfish", "salmon"]

__doc__ = '''top level workflow for linear RNA profiling,
associated key: {key_in_global} in [GLOBAL],

'''.format(key_in_global=_OPT_KEY_MAPPER + "/" + _OPT_KEY_QUANTIFIER)
__author__ = 'zerodel'


def main(path_config=""):
    user_config_whole = pysrc.body.config.config(
        path_config) if path_config else pysrc.body.config.load_default_value()

    quantifier = get_value_from_global_section(user_config_whole, _OPT_KEY_QUANTIFIER)

    if _is_this_quantifier_need_external_mapper(quantifier):
        align.work(user_config_whole, get_value_from_global_section(user_config_whole, _OPT_KEY_MAPPER))

    quantify.work(user_config_whole, quantifier)


def _is_this_quantifier_need_external_mapper(name_quantifier):
    # now we simple assume that sailfish do not need external mapper
    return name_quantifier.strip() not in _QUANTIFIERS_DO_NOT_NEED_EXTERNAL_MAPPER


def get_value_from_global_section(user_config_whole, key_name):
    if key_name in user_config_whole[pysrc.body.config.SECTION_GLOBAL]:
        return user_config_whole[pysrc.body.config.SECTION_GLOBAL][key_name]
    else:
        raise KeyError("Error@config_global_section: {} must be in GLOBAL section".format(key_name))


def __cli_arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cfg_file", help="file path to a configuration file of isoform detection")
    parser.add_argument("-l", "--log_file", help="logging file path", default="")
    return parser


if __name__ == "__main__":
    arg_parser = __cli_arg_parser()
    args = arg_parser.parse_args()
    main(args.cfg_file)
