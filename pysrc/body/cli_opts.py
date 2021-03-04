# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import argparse
import os

import pysrc.body.config
import pysrc.body.logger

__doc__ = '''
'''
__author__ = 'zerodel'

_module_logger = pysrc.body.logger.default_logger("commandline_option")


def drop1(dict_para, key, template):
    return template.format(dict_para.pop(key))


def option_and_value(key, dict_para):
    if dict_para[key]:
        return "{} {}".format(key, dict_para[key])
    else:
        return key


def drop_key(key, dict_para):
    if key in dict_para:
        if dict_para[key]:
            return "{} {}".format(key, dict_para.pop(key))
        else:
            return key
    else:
        return ""


def cat_options_no_replace(options, opts_value_dict):
    return " ".join([drop_key(key, opts_value_dict)
                     for key in options
                     if key in opts_value_dict])


def cat_options(options, opts_value_dict):
    return " ".join([option_and_value(key, opts_value_dict)
                     for key in options
                     if key in opts_value_dict])


def all_options(options):
    return " ".join([option_and_value(key, options) for key in options])


def enum_all_opts(opts_value_dict):
    return " ".join([option_and_value(key, opts_value_dict) for key in opts_value_dict])


def update_parameters(para_default, para_cli, para_conf):
    import copy
    tmp_para = copy.copy(para_default)
    tmp_para.update(para_cli)
    tmp_para.update(para_conf)
    return tmp_para


# todo: [2017_11_06 15:45] seems this will cause an bug :
# when your want to update the parameters in a new config file , it wont work as you expected
# it ignored META and GLOBAL in the newer config file
def merge_parameters(kwargs, parameters_custom, section_name_in_config):
    config_default = pysrc.body.config.load_default_value()
    # parameters_custom = config_custom[section_name_in_config] if section_name_in_config in config_custom else {}
    parameters_default = dict(config_default[section_name_in_config]) if section_name_in_config in config_default else {}
    updated_para = update_parameters(parameters_default,
                                     parameters_custom,
                                     kwargs) if parameters_custom else parameters_default
    return updated_para


def chain_map(*args):
    result = {}
    for single_map in args:
        result.update(single_map)
    return result


def check_if_these_files_exist(filename_str):
    single_paths = [x for x in filename_str.split() if x]
    for fa in single_paths:
        if not os.path.exists(fa):
            return False
    return True


def is_suitable_path_with_prefix(path_prefix):
    raw_path, prefix = os.path.split(path_prefix)
    return os.path.exists(raw_path) and os.path.isdir(raw_path)


def transform_input_general(input_files):
    if isinstance(input_files, list):
        paths = " ".join(input_files)
    elif isinstance(input_files, str):
        paths = input_files
    else:
        raise TypeError("Error: only list and string are allowed in assign input_files")
    return paths


def simple_arg_parser_only_config_option():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="specify the config file for this align job")

    return parser


def catch_one(opts_dict, *args):
    for arg in args:
        if arg in opts_dict:
            return opts_dict.get(arg)
    else:
        raise KeyError("Error: no such key %s in option %s" % ("/".join(args), opts_dict))


def extract_one(opts, key, default=""):
    val = opts.pop(key, default)
    if key in opts:
        raise KeyError(
            "ERROR@OPTION: unable to remove {this_entry} in a dict using dict.pop".format(this_entry=key))
    return val


def extract_entries_from_value_str(tools_str_in_config):
    str_tools = tools_str_in_config.strip()
    if "," in str_tools:
        tz = [x.strip() for x in str_tools.split(",")]
    else:
        tz = [x.strip() for x in str_tools.split()]
    return tz
