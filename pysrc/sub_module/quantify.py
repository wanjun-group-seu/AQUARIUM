# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import pysrc.wrapper.rsem
import pysrc.wrapper.sailfish

import pysrc.body.cli_opts
import pysrc.body.config
import pysrc.wrapper.salmon

__doc__ = '''
'''
__author__ = 'zerodel'

available_tools = {
    "sailfish": pysrc.wrapper.sailfish,
    "rsem": pysrc.wrapper.rsem,
    "salmon": pysrc.wrapper.salmon
}

_OPT_KEY_NAME_QUANTIFIER = "quantifier"


def work(whole_config_content=None, quantifier_name=None, inputs=None):
    quantifier = available_tools[quantifier_name]
    setting_index = dict(whole_config_content[quantifier.SECTION_INDEX])
    setting_quantify = dict(whole_config_content[quantifier.SECTION_QUANTIFY])

    index_quantify = _prepare_index(setting_index, quantifier)
    result = _do_quantify(setting_quantify, quantifier, index_quantify, inputs)

    return result


def _prepare_index(index_config, quantifier):
    index_path = quantifier.get_index_path(index_config)
    try:
        if quantifier.is_path_contain_index(index_path):
            return index_path
        else:
            quantifier.index(index_config)
    except FileNotFoundError as e:
        raise e


def _do_quantify(quantify_config, quantifier, quantify_ref, inputs):
    config_dict_quantify = pysrc.body.cli_opts.chain_map(quantify_config,
                                                         quantifier.interpret_index_path(quantify_ref),
                                                         quantifier.interpret_seq_files(inputs))

    return quantifier.quantify(config_dict_quantify)


def main(config_file):
    config_latter = pysrc.body.config.config(config_file)
    quantifier_name = config_latter[pysrc.body.config.SECTION_GLOBAL][_OPT_KEY_NAME_QUANTIFIER]
    work(config_latter, quantifier_name)


if __name__ == "__main__":
    arg_parser = pysrc.body.cli_opts.simple_arg_parser_only_config_option()
    par = arg_parser.parse_args()
    main(par.config)
