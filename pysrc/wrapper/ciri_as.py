# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import copy
import os
import os.path

import pysrc.body.cli_opts
import pysrc.body.logger
import pysrc.body.option_check
import pysrc.body.worker
import pysrc.wrapper.ciri

_OPT_CIRI_AS_PATH = "ciri_as_path"

__doc__ = '''this is the wrapper of CIRI-AS, attention: this wrapper will output all information by default .
'''
__author__ = 'zerodel'

SECTION_DETECT = "CIRI_AS"

# _ESSENTIAL_ARGUMENTS = ["--sam", "--ciri", "--out", "--ref_dir", "--ref_file", "--anno", "--log"]

_ARGUMENT_ORDER = ["--sam", "--ciri", "--out", "--ref_dir", "--ref_file", "--anno", "--output_all", "--log"]

_logger = pysrc.body.logger.default_logger(SECTION_DETECT)


def _is_a_suitable_file_path(path_given):
    p_dir, p_file = os.path.split(path_given)
    return os.path.exists(p_dir) and os.path.isdir(p_dir)


def _check_opts(opts=None):
    check_your_option = pysrc.body.option_check.OptionChecker(opts, name=SECTION_DETECT)

    check_your_option.must_have(_OPT_CIRI_AS_PATH, os.path.exists,
                                FileNotFoundError("Error: can not find CIRI-AS script"),
                                "path to CIRI-AS script file ")

    check_your_option.must_have("--sam", os.path.exists,
                                FileNotFoundError(
                                    "Error: unable to find CIRI-AS input sam file"),
                                "input sam file , should be the same as the CIRI used")

    check_your_option.must_have("--ciri", os.path.exists,
                                FileNotFoundError("Error@CIRI-AS : unable to find CIRI output file"),
                                "CIRI output file , should be the same version with CIRI AS")

    check_your_option.must_have("--out", pysrc.body.cli_opts.is_suitable_path_with_prefix,
                                FileNotFoundError("Error@CIRI-AS: incorrect output file for CIRI-AS"),
                                "output file path for CIRI-AS")

    check_your_option.one_and_only_one(["--ref_dir", "--ref_file"], pysrc.wrapper.ciri.check_ref,
                                       FileNotFoundError("Error@CIRI-AS: unable to find ref-file for CIRI-AS"),
                                       "genomic reference file, should be the same as CIRI")

    check_your_option.may_need("--anno", os.path.exists,
                               FileNotFoundError("Error@CIRI-AS: incorrect annotation file provide for CIRI-AS"),
                               "genomic annotation , should be the same as CIRI")

    check_your_option.may_need("--log", _is_a_suitable_file_path,
                               FileNotFoundError("Error@CIRI-AS: incorrect path for a log file "),
                               "output log file name (optional)")

    check_your_option.forbid_these_args("--help", "-H")
    return check_your_option


opt_checker = _check_opts()  # set up the opt_checker

OPTION_CHECKERS = [opt_checker]


def detect(para_config=None, **kwargs):
    opts_raw = pysrc.body.cli_opts.merge_parameters(kwargs, para_config, SECTION_DETECT)

    _logger.debug("ciri-as args: %s" % str(opts_raw))

    opt_checker.check(copy.copy(opts_raw))

    opts_full_output = copy.copy(opts_raw)
    opts_full_output["--output_all"] = "yes"  # force to output all
    cmd_detect = _get_detect_cmd(opts_full_output)

    _logger.debug("ciri-as command is : %s" % cmd_detect)

    pysrc.body.worker.run(cmd_detect)

    return opts_raw


def _get_detect_cmd(opts):
    cmd_as = "perl {ciri_as_path}".format(ciri_as_path=opts.pop(_OPT_CIRI_AS_PATH))
    cmd_latter = " ".join([pysrc.body.cli_opts.drop_key(key, opts) for key in _ARGUMENT_ORDER if key in opts])

    cmd_other = pysrc.body.cli_opts.enum_all_opts(opts)
    return " ".join([cmd_as, cmd_latter, cmd_other])


def to_bed():
    pass


def interpret_seq_files():
    pass


if __name__ == "__main__":
    print(__doc__)
    print(opt_checker)
