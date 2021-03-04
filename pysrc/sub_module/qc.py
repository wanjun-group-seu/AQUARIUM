# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import copy


import pysrc.body.cli_opts
import pysrc.body.config
import pysrc.body.logger
import pysrc.wrapper.fastqc

__doc__ = '''
'''
__author__ = 'zerodel'

available_tools = {
    "fastqc": pysrc.wrapper.fastqc
}

_OPT_KEY_NAME_QC_TOOL = "qc_tool"

_logger = pysrc.body.logger.default_logger("QC")


def work(whole_config_content=None, qc_tool_name=None, output_path=None, input_files=None):
    qc_tool = available_tools[qc_tool_name]
    setting_qc = dict(whole_config_content[qc_tool.SECTION_QC_SETTING])
    _logger.debug("setting for {tool} is {setting}".format(tool=qc_tool_name,
                                                           setting=str(setting_qc)))
    qc_tool.run_qc(copy.copy(setting_qc))


def main(config_file):
    config_latter = pysrc.body.config.config(config_file)
    qc_tool_name = config_latter[pysrc.body.config.SECTION_GLOBAL][_OPT_KEY_NAME_QC_TOOL]
    _logger.info("use %s as QC tool " % qc_tool_name)

    work(config_latter, qc_tool_name)


