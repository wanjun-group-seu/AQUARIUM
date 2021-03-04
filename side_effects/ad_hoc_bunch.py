# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#

import argparse
import copy
import os

import pysrc.body.cli_opts
import pysrc.body.config
import pysrc.body.logger
import pysrc.body.worker
import pysrc.sub_module.align
import wf_detect_circRNA
import wf_profile_circRNA
import wf_profile_linear
import wf_reveal_circ_isoform
from pysrc.sub_module import align, quantify

__doc__ = ''' python ad_hoc_bunch.pysrc result_folder_path sra_files_pat
'''
__author__ = 'zerodel'

own_cfg = {
    "CIRI": {
        # "-T": "10",
        "--in": "",
        "--seqs": "",
        "--out": "",
        "-P": None,
    },
    "GLOBAL": {
        "detector": "ciri",
    },
    "CIRC_PROFILE": {
        "quantifier": "salmon",
        #        "--annotation": "",
        #        "--genomic_seqs_fasta": "",
        "-o": "",
        "-k": "",
        #        "--mll": ""
    },
    "KNIFE": {
        "dataset_name": "knife",
        "read_directory": "",
        "alignment_parent_directory": ""
    }

}
#
# [GLOBAL]
# detector = knife
# [KNIFE]
# read_directory = /home/zerodel/zData/parts_blood/profile_sra/ERR335312/fq
# alignment_parent_directory = /home/zerodel/zData/parts_blood/profile_sra/ERR335312/detection_report
# dataset_name = knife

WORKING_FLOWS_OF = {"align": align,
                    "quantify": quantify,
                    "profile_linear": wf_profile_linear,
                    "BSJ_detection": wf_detect_circRNA,
                    "AS_detection": wf_reveal_circ_isoform,
                    "profile_circRNA": wf_profile_circRNA}

JOB_ID_SECTION = "jobs"
WORKING_PATH_SECTION = "working_path"

_logger = pysrc.body.logger.default_logger("BUNCH_FOR_SRA")


def make_parameters_for_this_job(para_dict, job_id, tap_root):
    dict_args = copy.copy(para_dict)
    abs_path_up_level = os.path.abspath(tap_root)

    def _put_it_under(x):
        return os.path.join(abs_path_up_level, x(job_id))

    dict_args["CIRI"]["--in"] = os.path.join(_put_it_under(_bwa_sam_path), "pe.sam")
    dict_args["CIRI"]["--out"] = os.path.join(_put_it_under(_detection_report_path), "ciri.out")

    path_to_fq_files = os.listdir(_put_it_under(_fq_path))
    fq1 = [f for f in path_to_fq_files if f.endswith("1.fastq") or f.endswith("1.fq")][0]
    fq2 = [f for f in path_to_fq_files if f.endswith("2.fastq") or f.endswith("2.fq")][0]

    dict_args["CIRI"]["--seqs"] = " ".join([os.path.join(_put_it_under(_fq_path), fq1),
                                            os.path.join(_put_it_under(_fq_path), fq2)])

    dict_args["CIRC_PROFILE"]["-o"] = _put_it_under(_quantify_result_path)
    dict_args["CIRC_PROFILE"]["-k"] = 31

    dict_args["KNIFE"]["read_directory"] = _put_it_under(_fq_path)
    dict_args["KNIFE"]["alignment_parent_directory"] = _put_it_under(_detection_report_path)

    _logger.debug("parameters: %s" % str(dict_args))

    return dict_args


def main(folder_of_result_file, folder_of_sra_files, work_flow, cfg_file):
    if not (folder_of_result_file and folder_of_result_file and work_flow and cfg_file):
        raise KeyError("ERROR: you should give some arguments")

    folder_of_sra_files = os.path.abspath(folder_of_sra_files)
    folder_of_result_file = os.path.abspath(folder_of_result_file)

    _logger.debug(
        "samples in {sra}, result will be in {out} ".format(sra=folder_of_sra_files, out=folder_of_result_file))

    sra_ids = guess_sra_ids(folder_of_sra_files)

    _logger.debug("there are {} samples".format(str(len(sra_ids))))

    temp_dict = own_cfg if not cfg_file else pysrc.body.config.cfg2dict(pysrc.body.config.config(cfg_file))

    for sra_id in sra_ids:
        _logger.debug("staring %s" % str(sra_id))

        prepare_sub_path_for(sra_id, folder_of_result_file)

        _logger.debug("extract sra %s ...." % sra_id)

        extract_fq(sra_id, folder_of_result_file, folder_of_sra_files)

        para_this_job = make_parameters_for_this_job(temp_dict, sra_id, folder_of_result_file)

        _logger.debug("final parameter is : %s" % str(para_this_job))

        work_flow.main(para_this_job)


def guess_sra_ids(sra_path):
    return [x.split(".")[0] for x in os.listdir(sra_path) if x.endswith(".sra")]


def _get_sra_file_path(sra, sra_root):
    return os.path.join(os.path.abspath(sra_root), sra + ".sra")


def extract_fq(sra_id, up_level_path, sra_root):
    fq_path = os.path.join(up_level_path, _fq_path(sra_id))
    sra_path = _get_sra_file_path(sra_id, sra_root)
    # WARNING : this code snippet need a specific default config file
    default_setting = pysrc.body.config.load_default_value()
    fastq_dump_bin = default_setting["META"]["fastq_dump_bin"] if "fastq_dump_bin" in default_setting[
        "META"] else "fastq-dump"

    if not _is_paired_end_fq_extracted(sra_id, fq_path):
        pysrc.body.worker.run("{fq_dump_bin} --split-files {sra} -O {fq}".format(fq_dump_bin=fastq_dump_bin,
                                                                                 sra=sra_path,
                                                                                 fq=fq_path))


def _is_paired_end_fq_extracted(sra_id, path_fq):
    fqs = [f for f in os.listdir(path_fq) if f.endswith(".fq") or f.endswith(".fastq")]
    this_sra_fq = [f.strip().split(".")[0] for f in fqs if f.startswith(sra_id)]

    mate1ok = any([f.endswith("1") for f in this_sra_fq])
    mate2ok = any([f.endswith("2") for f in this_sra_fq])
    return mate1ok and mate2ok


def prepare_sub_path_for(job_id, up_level_path):
    func_list = [_sra_path, _fq_path, _bwa_sam_path, _detection_report_path, _quantify_result_path]
    abs_path_up_level = os.path.abspath(up_level_path)

    def prepare_folder(x):
        _make_it_exist(os.path.join(abs_path_up_level, x))

    for func in func_list:
        prepare_folder(func(job_id))


def _make_it_exist(path_to):
    if not os.path.exists(path_to):
        os.mkdir(path_to)


def _sra_path(sra):
    return sra.strip()


def _fq_path(sra):
    return os.path.join(sra, "fq")


def _bwa_sam_path(sra):
    return os.path.join(sra, "sam")


def _detection_report_path(sra):
    return os.path.join(sra, "detection_report")


def _quantify_result_path(sra):
    return os.path.join(sra, "quantify_result")


def get_sra_id_in_config(cfg):
    if pysrc.body.config.SECTION_GLOBAL in cfg:
        job_ids = [x.strip() for x in cfg[pysrc.body.config.SECTION_GLOBAL][JOB_ID_SECTION].strip().split(" ")
                   if x]
        return job_ids
    else:
        raise KeyError("Error@bunch_work: no GLOBAL section in configuration file")


def _get_tap_root_from_cfg(cfg):
    if pysrc.body.config.SECTION_GLOBAL in cfg:

        part_global = cfg[pysrc.body.config.SECTION_GLOBAL]
        if WORKING_PATH_SECTION in part_global:
            return part_global[WORKING_PATH_SECTION]
        else:
            raise KeyError("Error@'GLOBAL' in configure : WORKING_PATH not found ")
    else:
        raise KeyError("Error@ConfigureFIle: no 'GLOBAL' part in configure file. ")


def _make_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sra", help="path to your sra samples folder", default="")
    parser.add_argument("--out", help="path to your output path, each sra will have a sub-dir under this path",
                        default="")
    parser.add_argument("--config", help="path to your config file", default="")
    parser.add_argument("--work_flow",
                        help="choose your work flow : align , quantify , BSJ_detection, AS_detection, profile_linear, "
                             "profile_circRNA",
                        choices=[x for x in WORKING_FLOWS_OF], default="")

    return parser


if __name__ == "__main__":
    parser_this_bunch = _make_arg_parser()
    par = parser_this_bunch.parse_args()
    # if par.work_flow not in WORKING_FLOWS_OF:
    #     raise KeyError("ERROR: @ parameter , incorrect workflow name, allowed workflow names: {}".format(
    #         ", ".join([x for x in WORKING_FLOWS_OF])))
    main(par.out, par.sra, WORKING_FLOWS_OF[par.work_flow], par.config)
