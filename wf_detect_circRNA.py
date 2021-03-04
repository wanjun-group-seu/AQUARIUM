# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import argparse
import copy

import pysrc.being.circ_region
import pysrc.body.cli_opts
import pysrc.body.config

import pysrc.body.logger
import pysrc.body.option_check
import pysrc.body.utilities
import pysrc.wrapper.bwa
import pysrc.wrapper.ciri
import pysrc.wrapper.ciri2
import pysrc.wrapper.ciri_as
import pysrc.wrapper.knife
import pysrc.wrapper.ciri_full

_OPT_FINAL_BED_OUT = "bed_output"

_OPT_DETECTOR = "detector"

__TOOL_KNIFE = "knife"

__TOOL_CIRI_AS = "ciri_as"

__TOOL_CIRI = "ciri"

__TOOL_CIRI_FULL = "ciri_full"

__TOOL_BWA = "bwa"

SECTION_CIRC_DETECTION = 'CIRC_DETECTION'

available_tools = {
    __TOOL_BWA: pysrc.wrapper.bwa,
    __TOOL_CIRI: pysrc.wrapper.ciri2,
    __TOOL_CIRI_AS: pysrc.wrapper.ciri_as,
    __TOOL_KNIFE: pysrc.wrapper.knife,
    __TOOL_CIRI_FULL: pysrc.wrapper.ciri_full
}

__tools_for_detection = [__TOOL_CIRI, __TOOL_CIRI_AS, __TOOL_KNIFE, __TOOL_CIRI_FULL]

__doc__ = ''' top level interface of circRNA detection workflow.\n
choose detector:  '{key_global}' in section [{section_name}]\n
usable value : {usable_values}  or combined 
'''.format(key_global=_OPT_DETECTOR,
           section_name=SECTION_CIRC_DETECTION,
           usable_values=", ".join([x for x in available_tools]))

__author__ = 'zerodel'

_logger = pysrc.body.logger.default_logger(SECTION_CIRC_DETECTION)


def __get_option_checker(opts=None):
    oc = pysrc.body.option_check.OptionChecker(opts, name=SECTION_CIRC_DETECTION)

    def __check_tools_string(tools_str_in_config):
        return all([x in __tools_for_detection for x in
                    (pysrc.body.cli_opts.extract_entries_from_value_str(tools_str_in_config))])

    def __is_bed_file_name_usable(path_bed):
        return pysrc.body.utilities.is_path_creatable(path_bed) and path_bed.strip().endswith(".bed")

    oc.must_have("detector",
                 check_fun=__check_tools_string,
                 exception_on_error=LookupError("unable to find this circ-detection tool"),
                 des_str="specify the circ_rna detection tool, current support {tz}".format(
                     tz=",".join(__tools_for_detection)))

    oc.may_need("log_file", check_fun=pysrc.body.utilities.is_path_creatable,
                exception_on_error=FileNotFoundError("unable to set logfile for {}".format(SECTION_CIRC_DETECTION)),
                des_str="path to log file of {}".format(SECTION_CIRC_DETECTION))

    oc.may_need(_OPT_FINAL_BED_OUT, check_fun=__is_bed_file_name_usable,
                exception_on_error=FileNotFoundError("bed output file name not suitable"),
                des_str="path to final output in .bed format")

    return oc


opt_checker = __get_option_checker()
OPTION_CHECKERS = [opt_checker]


def __cli_arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cfg_file", help="file path to a configuration file of detection job")
    return parser


def main(cfg):
    _logger.info("circRNA Detecting starting\n" + "configure is a {} , and content is {}".format(type(cfg), str(cfg)))

    detect_option_section = pysrc.body.config.load_or_update_option_section(SECTION_CIRC_DETECTION, cfg)
    opt_checker.check(copy.copy(detect_option_section))

    detector_names = pysrc.body.cli_opts.extract_entries_from_value_str(detect_option_section[_OPT_DETECTOR])
    _logger.info("using {num_detector} detectors : {detector_str} ".format(num_detector=len(detector_names),
                                                                           detector_str=", ".join(detector_names)))

    bed_file_paths = {}
    for single_detector in detector_names:
        _logger.info("using %s as detector" % single_detector)
        opts_this = _detect_by(cfg, single_detector)
        bed_result = opts_this.get("bed_out", "")
        if bed_result:
            bed_file_paths[single_detector] = bed_result

    combined_bed_path = pysrc.body.cli_opts.extract_one(detect_option_section,
                                                        _OPT_FINAL_BED_OUT,
                                                        "")
    if combined_bed_path:
        merge_all_bed(bed_file_paths, combined_bed_path)


def merge_all_bed(dd_path_bed, out_bed):
    whole_regions = []
    for single_detector in dd_path_bed:
        bed_this = dd_path_bed.get(single_detector, "")
        if bed_this:
            with open(bed_this) as read_bed:
                region_ids = [bed_line_to_region(line) for line in read_bed]
                whole_regions.extend(region_ids)
    whole_regions = sorted(list(set(whole_regions)))

    with open(out_bed, "w") as flush_it:
        for region_circ in whole_regions:
            tmp_rc = pysrc.being.circ_region.RegionCirc(id_str=region_circ)
            flush_it.write("{}\n".format(tmp_rc.as_bed()))


def bed_line_to_region(line):
    parts = line.strip().split()[:3]
    return "{chr}:{start}|{end}".format(chr=parts[0], start=parts[1], end=parts[2])


def _detect_by(cfg, name_of_detector):
    if name_of_detector not in available_tools:
        raise KeyError("Error: no such circular RNA detection tool : {}".format(name_of_detector))
    detector = available_tools[name_of_detector]
    config_detector = pysrc.body.config.load_or_update_option_section(detector.SECTION_DETECT, cfg)

    _logger.info("content of detector is =====\n{}\n".format(str(config_detector)))

    return detector.detect(config_detector)


if __name__ == "__main__":
    arg_parser = __cli_arg_parser()
    args = arg_parser.parse_args()
    main(args.cfg_file)
