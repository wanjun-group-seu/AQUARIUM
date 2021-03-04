# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
import pysrc.body.config
import pysrc.wrapper.bwa
import pysrc.wrapper.ciri
import pysrc.wrapper.ciri2
import pysrc.wrapper.ciri_full
import pysrc.wrapper.knife
import pysrc.wrapper.rsem
import pysrc.wrapper.sailfish
import pysrc.wrapper.salmon
import pysrc.wrapper.star
import pysrc.wrapper.fastqc

import pysrc.wrapper.ciri_as
import wf_profile_circRNA
import wf_detect_circRNA

__doc__ = '''
'''
__author__ = 'zerodel'

_tool_wrapper = {}
_tool_description = {}


def show_it(tool_name):
    try:
        wrapper_tool = _tool_wrapper[tool_name]
    except KeyError as e:
        print("no such tool as %s" % tool_name)
        raise e

    module_desc = wrapper_tool.__doc__
    module_desc = "\n".join(["# %s" % line for line in module_desc.strip().split('\n')])
    print(module_desc)

    default_config = pysrc.body.config.load_default_value()
    for checker in wrapper_tool.OPTION_CHECKERS:
        default_this_section = dict(default_config[checker.name]) if checker.name in default_config else {}
        checker.dict_opt = default_this_section
        print(checker)


def reveal():
    print("\nthis pipeline now has following wrappers, choose one to see more information \n")
    keys = sorted([k for k in _tool_description])
    for key in keys:
        if key in _tool_description:
            print("%s \t-\t %s " % (key, _tool_description[key]))

    print("\n")
    print(
        r""" you could use 'grep -v "^#\ " | grep -v "^$"' to filter out the options and use them in config file """)


def _add_tool(name, wrapper, description):
    _tool_wrapper[name] = wrapper
    _tool_description[name] = description


_add_tool("fastqc", pysrc.wrapper.fastqc, "fastqc , a QC tool")
_add_tool("bwa", pysrc.wrapper.bwa, "BWA and BWA MEM ")
_add_tool("ciri", pysrc.wrapper.ciri, "CIRI : a circular RNA detection tool ")
_add_tool("ciri2", pysrc.wrapper.ciri2, "CIRI2 : a circular RNA detection tool, multiple core empowered ")
_add_tool("ciri_as", pysrc.wrapper.ciri_as, "CIRI-AS: circular RNA Alternative Splicing Event detection tool")
_add_tool("ciri_full", pysrc.wrapper.ciri_full, "CIRI-FULL: a powerful circular RNA detection and rebuilding tool, "
                                                "which combines CIRI and CIRI-AS")
_add_tool("knife", pysrc.wrapper.knife, "KNIFE: a circular RNA detection tool ")
_add_tool("rsem", pysrc.wrapper.rsem, "RSEM : a RNA seq quantification tool")
_add_tool("sailfish", pysrc.wrapper.sailfish, "Sailfish: a RNA-seq quantification tool based on k-mer")
_add_tool("salmon", pysrc.wrapper.salmon, "Salmon: a RNA-seq quantification tool based on fragment ")
_add_tool("star", pysrc.wrapper.star, "STAR : a junction sensitive aligner")
_add_tool("profile_circRNA", wf_profile_circRNA, "home made pipeline for profiling the circRNA ")
_add_tool("detect_circRNA", wf_detect_circRNA, "home made pipeline for detecting circRNA")

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        reveal()
    else:
        show_it(sys.argv[-1])
