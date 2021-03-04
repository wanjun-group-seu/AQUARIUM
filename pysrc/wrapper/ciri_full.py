# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
#
from pysrc.body.utilities import safe_rename_with_postfix


__author__ = 'zerodel'

import os
import os.path
import glob
import copy

import pysrc.body.config
import pysrc.body.cli_opts
import pysrc.body.logger
import pysrc.body.option_check
import pysrc.body.worker
import pysrc.body.utilities as ut
import pysrc.file_format.ciri_entry
import pysrc.file_format.fa
import pysrc.file_format.bsj_gtf
import pysrc.file_format.gtf

import pysrc.wrapper.ciri2
import pysrc.wrapper.ciri_as
import pysrc.wrapper.bwa

_OPT_CIRI_AS_PATH = "ciri_as_path"
_OPT_CIRI_PATH = "ciri_path"
_OPT_BWA_BIN = "bwa_bin"

_OPT_JAR_FULL_PATH = "jar_full"
_OPT_JAR_VIS_PATH = "jar_vis"

BASE_NAME_R_SCRIPT_TRANSFORM_LIST = "transform_ciri_vis_list_exon.R"
BASE_NAME_R_SCRIPT_FILTER_BSJ = "filter_untouched_bsj.R"

__doc__ = '''this is the wrapper of CIRI-AS, attention: this wrapper will output all information by default .
'''

SECTION_DETECT = "CIRI_FULL"

# _ESSENTIAL_ARGUMENTS = ["--sam", "--ciri", "--out", "--ref_dir", "--ref_file", "--anno", "--log"]

_ARGUMENT_ORDER = ["--sam", "--ciri", "--out", "--ref_dir",
                   "--ref_file", "--anno", "--output_all", "--log"]

_logger = pysrc.body.logger.default_logger(SECTION_DETECT)

_OPT_SHOW_ALL = "--no_strigency"


def _is_a_suitable_file_path(path_given):
    p_dir, p_file = os.path.split(path_given)
    return os.path.exists(p_dir) and os.path.isdir(p_dir)


def _check_opts(opts=None):
    check_your_option = pysrc.body.option_check.OptionChecker(
        opts, name=SECTION_DETECT)

    check_your_option.must_have(_OPT_CIRI_AS_PATH, os.path.exists,
                                FileNotFoundError(
                                    "Error: can not find CIRI-AS script"),
                                "path to CIRI-AS script file ")

    check_your_option.must_have(_OPT_CIRI_PATH, os.path.exists,
                                FileNotFoundError(
                                    "Error: can not find CIRI script"),
                                "path to CIRI script file ")

    check_your_option.must_have(_OPT_BWA_BIN, os.path.exists,
                                FileNotFoundError(
                                    "Error: can not find bwa executive file"),
                                "path to bwa binary")

    check_your_option.must_have(_OPT_JAR_VIS_PATH, os.path.exists,
                                FileNotFoundError(
                                    "Error: can not find CIRI-Vis jar file"),
                                "path to CIRI_Vis.jar ")
    check_your_option.must_have(_OPT_JAR_FULL_PATH, os.path.exists,
                                FileNotFoundError(
                                    "Error: can not find CIRI-full jar file"),
                                "path to CIRI-full.jar")

    check_your_option.must_have("-1", os.path.exists,
                                FileNotFoundError(
                                    "Error: unable to find pair-end reads file 1 "),
                                "path to PE reads 1")

    check_your_option.must_have("-2", os.path.exists,
                                FileNotFoundError(
                                    "Error: unable to find pair-end reads file 2 "),
                                "path to PE reads 2")

    check_your_option.must_have("-a", os.path.exists,
                                FileNotFoundError(
                                    "Error@CIRI-AS: incorrect annotation file provide for CIRI-full"),
                                "genomic annotation file")

    check_your_option.must_have("-r", os.path.exists,
                                FileNotFoundError(
                                    "Error@CIRI-full: genomic reference .fa file not found"),
                                "genomic reference file")

    check_your_option.must_have("-d", ut.is_path_a_legal_dir,
                                NotADirectoryError(
                                    "Error@CIRI-full: improper output directory given"),
                                "path to output files")

    check_your_option.must_have("-l", lambda x: isinstance(int(x), int),
                                TypeError(
                                    "Error@CIRI-full , improper read length given"),
                                "read length")

    check_your_option.may_need("-t", lambda x: isinstance(int(x), int),
                               TypeError(
                                   "Error@CIRI-full , improper thread number given"),
                               "number of threads")

    check_your_option.may_need("-T", lambda x: isinstance(int(x), int),
                               TypeError(
                                   "Error@CIRI-full , improper mapping threshold given"),
                               "bwa mem mapping threshold, for human, 19 is ok, more sensitive for lower threshold")

    check_your_option.may_need("circ_fa", ut.is_path_creatable,
                               FileNotFoundError("Error@CIRI-full: given path is not valid for reconstructed fasta "
                                                 "file"), "path for reconstructed fa file ")

    check_your_option.may_need(_OPT_SHOW_ALL, lambda x: True,
                               ValueError("Error@CIRIFULL, you need a flag"),
                               """a flag whether to show all possible BSJ"""
                               )

    check_your_option.may_need("-min", lambda x: int(x) > 0,
                               ValueError(
                                   "Error@CIRIFULL, minimal expression should be a integer larger than 1"),
                               """ the minimal expression in step ciri-vis""")
    # check_your_option.must_have("--sam", os.path.exists,
    #                             FileNotFoundError(
    #                                 "Error: unable to find CIRI-full input sam file"),
    #                             "input sam file , should be the same as the CIRI used")
    #
    # check_your_option.must_have("--ciri", os.path.exists,
    #                             FileNotFoundError("Error@CIRI-AS : unable to find CIRI output file"),
    #                             "CIRI output file , should be the same version with CIRI AS")

    check_your_option.forbid_these_args("--help", "-H")
    return check_your_option


opt_checker = _check_opts()  # set up the opt_checker

OPTION_CHECKERS = [opt_checker]


def detect(para_config=None, **kwargs):
    opts_raw = pysrc.body.cli_opts.merge_parameters(
        kwargs, para_config, SECTION_DETECT)

    _logger.debug("ciri-full args: %s" % str(opts_raw))

    opt_checker.check(copy.copy(opts_raw))

    # prepare common arguments
    r1, r2, ref, anno, = opts_raw.get(
        "-1"), opts_raw.get("-2"), opts_raw.get("-r"), opts_raw.get("-a")

    read_length = opts_raw.get("-l")

    bwa_score = opts_raw.get("-T")

    min_expr = opts_raw.get("-min", 1)

    ciri_path, ciri_as_path, bwa_path = opts_raw.get(_OPT_CIRI_PATH), opts_raw.get(_OPT_CIRI_AS_PATH), opts_raw.get(
        _OPT_BWA_BIN)

    jar_full, jar_vis = opts_raw.get(
        _OPT_JAR_FULL_PATH), opts_raw.get(_OPT_JAR_VIS_PATH)
    dir_out = opts_raw.get("-d")
    num_thread = opts_raw.get("-t", 1)
    fa_target = opts_raw.get("circ_fa")

    sub_dir_detection, sub_dir_vis = detection_sub_dir_under(
        dir_out), vis_sub_dir_under(dir_out)

    ut.do_make_dir(dir_out)
    ut.do_make_dir(sub_dir_detection)
    _logger.debug(
        "ciri-full detection report will be put in {}".format(sub_dir_detection))
    _logger.debug(
        "ciri-full visualization report will be put in {}".format(sub_dir_vis))

    # ut.do_make_dir(sub_dir_vis) # ciri-vis will make a empty dir itself, otherwise this will cause an error.

    ciri_file = ciri_report_under(dir_out)
    sam_file = os.path.join(sub_dir_detection, "align.sam")

    ciri_as_prefix = os.path.join(sub_dir_detection, "as")
    ciri_full_prefix = os.path.join(sub_dir_detection, "full")

    # run ciri2
    # get default setting dict
    ciri_args_dict = dict(pysrc.body.config.load_or_update_option_section(
        pysrc.wrapper.ciri2.SECTION_DETECT))

    # substitute those ciri arguments
    ciri_args_dict["bwa_bin"] = bwa_path
    ciri_args_dict["bwa_index"] = ref
    ciri_args_dict["ciri_path"] = ciri_path
    ciri_args_dict["--in"] = sam_file
    ciri_args_dict["--seqs"] = " ".join([r1, r2])
    ciri_args_dict["--out"] = ciri_file
    ciri_args_dict["--ref_file"] = ref
    ciri_args_dict["--anno"] = anno
    ciri_args_dict["--thread_num"] = num_thread
    ciri_args_dict["bwa_score"] = bwa_score

    if _OPT_SHOW_ALL in opts_raw:
        ciri_args_dict[_OPT_SHOW_ALL] = ""

    # now run it
    _logger.debug("now detecting circular RNA using CIRI")

    _logger.debug("and ciri running under the following  parameters: {}".format(
        str(ciri_args_dict)))

    pysrc.wrapper.ciri2.detect(ciri_args_dict)

    _logger.debug("CIRI part finished ....")

    if not pysrc.file_format.ciri_entry.is_ciri_file_intact(ciri_file):
        _logger.warning(
            "WARNING: CIRI may meet some ERROR , report file : {}".format(ciri_file))
    else:
        _logger.debug("seems ciri works well ")

    # run ciri-as
    ciri_as_args_dict = dict(pysrc.body.config.load_or_update_option_section(
        pysrc.wrapper.ciri_as.SECTION_DETECT))

    ciri_as_args_dict["ciri_as_path"] = ciri_as_path
    ciri_as_args_dict["--sam"] = sam_file
    ciri_as_args_dict["--ciri"] = ciri_file

    ciri_as_args_dict["--ref_file"] = ref
    ciri_as_args_dict["--anno"] = anno
    ciri_as_args_dict["--out"] = ciri_as_prefix

    _logger.debug(
        "CIRI-AS running under the parameters: {}".format(str(ciri_as_args_dict)))

    pysrc.wrapper.ciri_as.detect(ciri_as_args_dict)

    _logger.debug("starting RO1 phase")
    cmd_ro1 = _get_cmd_ro1(jar_full, r1, r2, ciri_full_prefix, num_thread)
    _logger.debug("command RO1 of ciri-full is : %s" % cmd_ro1)
    pysrc.body.worker.run(cmd_ro1)

    # phase RO2
    fq_ro1 = ciri_full_prefix + "_ro1.fq"
    sam_ro1 = ciri_full_prefix + "_ro1.sam"

    mapping_job_setting = {
        "bwa_bin": bwa_path,
        "bwa_index": ref,
        "read_file": fq_ro1,
        "sam": sam_ro1,
        "thread_num": num_thread
    }

    _logger.info("bwa mapping for RO2 with parameter: %s" %
                 str(mapping_job_setting))

    pysrc.wrapper.bwa.bwa_mapping_ciri_only(
        mapping_job_setting)  # perform the mapping using bwa
    _logger.debug("mapping completed")

    cmd_ro2 = _get_cmd_ro2(jar_full, ref, sam_ro1,
                           read_length, ciri_full_prefix)
    _logger.debug(
        "command RO2 of ciri-full is {cmd_ro2}".format(cmd_ro2=cmd_ro2))
    pysrc.body.worker.run(cmd_ro2)

    # phase Merge
    cmd_merge = _get_cmd_merge(
        jar_full, ref, anno, ciri_file, ciri_as_prefix, ciri_full_prefix)
    _logger.debug(
        "command merge of ciri-full is {cmd_merge}".format(cmd_merge=cmd_merge))
    pysrc.body.worker.run(cmd_merge)

    # phase Vis
    cmd_vis = _get_cmd_vis(
        jar_vis, ref, ciri_full_prefix, sub_dir_vis, min_expr)
    _logger.debug(
        "command vis of ciri-full is {cmd_vis}".format(cmd_vis=cmd_vis))
    pysrc.body.worker.run(cmd_vis)

    # after ciri-vis, try to extract the reconstructed FA
    if fa_target:
        fa_reconstructed = glob.glob(os.path.join(sub_dir_vis, "*.fa"))
        _logger.debug("after phase vis , we got fa files as : {}".format(
            str(fa_reconstructed)))
        _logger.debug(
            "try to combine reconstructed fa file into one file: {}....".format(fa_target))
        try:
            pysrc.file_format.fa.incremental_updating(fa_target, fa_reconstructed)
        except Exception as e:
            _logger.error("error happens while extracting circular RNA at {}".format(fa_target))
            
    return opts_raw


def ciri_report_under(dir_out):
    return os.path.join(detection_sub_dir_under(dir_out), "ciri.report")


def vis_sub_dir_under(dir_out):
    return os.path.join(dir_out, "vis")


def detection_sub_dir_under(dir_out):
    return os.path.join(dir_out, "detection")


def rebuild_fa_path_under(dir_out):
    your_sub_dir = vis_sub_dir_under(dir_out)

    fa_path = [x for x in os.listdir(your_sub_dir) if x.endswith(
        ".fa") or x.endswith(".fasta")]

    if fa_path:
        return os.path.join(your_sub_dir, fa_path[0])
    else:
        _logger.warning("""NO FA File Found in ciri-full result : {} , 
        we hope the stdout.list file of ciri-full do not mess it up """.format(your_sub_dir))
        return None


def vis_list_path_under(dir_out):
    your_sub_dir = vis_sub_dir_under(dir_out)
    list_path = [x for x in os.listdir(your_sub_dir) if x.endswith(".list")]
    return os.path.join(your_sub_dir, list_path[0])


def _get_cmd_ro1(jar_full, r1, r2, prefix, num_threads=None, minM=None, minI=None):
    pattern_cmd = "java -jar {jar_full} RO1 -1 {r1} -2 {r2} -o {prefix}"
    cmd = pattern_cmd.format(jar_full=jar_full,
                             r1=r1, r2=r2, prefix=prefix)

    if num_threads:
        cmd = " ".join([cmd, "-t", str(int(num_threads))])
    if minM:
        cmd = " ".join([cmd, "-minM", str(int(minM))])
    if minI:
        cmd = " ".join([cmd, "-minI", str(int(minI))])

    return cmd


def _get_cmd_ro2(jar_full, ref, sam, read_length, prefix, range_circ=None):
    cmd = "java -jar {jar_full} RO2 -r {ref} -s {sam} -l {read_length} -o {prefix}".format(jar_full=jar_full,
                                                                                           ref=ref,
                                                                                           sam=sam,
                                                                                           read_length=read_length,
                                                                                           prefix=prefix)
    if range_circ:
        cmd = " ".join([cmd, "-range", str(int(range_circ))])
    return cmd


def _get_cmd_merge(jar_full, ref, anno, ciri, ciri_as_prefix, ciri_full_prefix):
    ro2 = "".join([ciri_full_prefix, "_ro2_info.list"])
    as_out = "".join([ciri_as_prefix, "_jav.list"])
    cmd = "java -jar {jar_full} Merge -r {ref} -a {anno} -c {ciri} -as {ciri_as} -ro {ro2} -o {prefix}".format(
        jar_full=jar_full, ref=ref, anno=anno, ciri=ciri, ciri_as=as_out, ro2=ro2, prefix=ciri_full_prefix)
    return cmd


def _get_cmd_vis(jar_vis, ref, ciri_full_prefix, output_dir, min_expr):
    detail_anno = "".join([ciri_full_prefix, "_merge_circRNA_detail.anno"])
    cmd = "java -jar {jar_vis} -i {merged} -d {out_dir} -r {ref} -min {min_expr}".format(jar_vis=jar_vis,
                                                                                         merged=detail_anno,
                                                                                         out_dir=output_dir,
                                                                                         ref=ref,
                                                                                         min_expr=min_expr)
    return cmd


def to_bed():
    _logger.debug("this function : to_bed@ciri_full is not implemented yet")
    pass


def interpret_seq_files():
    _logger.debug(
        "this function : interpret_seq_files@ciri_full is not implemented yet")
    pass


def translate_vis_list(path_vis_list, tmp_dir=""):
    _logger.debug("be careful, this process will invoke en external R script")

    path_r_script = _find_that_r_script(BASE_NAME_R_SCRIPT_TRANSFORM_LIST)

    _logger.debug(
        "current path of external R script is {}".format(path_r_script))

    dir_par = tmp_dir if tmp_dir else os.path.dirname(path_vis_list)

    dir_parpar = os.path.dirname(dir_par)

    path_bed_circexon, path_bed_blank = os.path.join(
        dir_parpar, "circexon.bed"), os.path.join(dir_parpar, "blank.bed")

    # what if those two file already exists ? replace it with a new name.
    if os.path.exists(path_bed_circexon) or os.path.exists(path_bed_blank):
        time_stamp_now = pysrc.body.logger.raw_timestamp()

        _logger.debug("found bed file of previous run, rename them :  {}".format("[" + ". ".join([path_bed_blank,
                                                                                                  path_bed_circexon
                                                                                                  ]) + "]"))
        safe_rename_with_postfix(path_bed_circexon, time_stamp_now)
        safe_rename_with_postfix(path_bed_blank, time_stamp_now)

    _logger.debug("circexon bed file will be : {circexon},\\"
                  " blank bed will be : {blank}".format(circexon=path_bed_circexon,
                                                        blank=path_bed_blank))

    cmd_filter_r = " ".join(
        ["Rscript", path_r_script, path_vis_list, path_bed_circexon, path_bed_blank])

    _logger.debug(
        """external cmd for filter ciri-full vis list : \n {}""".format(cmd_filter_r))

    pysrc.body.worker.run(cmd_filter_r)

    return path_bed_circexon, path_bed_blank


def filter_out_un_touched_circular_rna(path_bed_ciri, path_vis_list, tmp_dir="", exon_only=True):
    _logger.debug(
        "start filter out those BSJ which contains partial inner structure information")

    path_r_script = _find_that_r_script(BASE_NAME_R_SCRIPT_FILTER_BSJ)

    _logger.debug(
        "external R script will be invoked : {}".format(path_r_script))

    dir_par = tmp_dir if tmp_dir else os.path.dirname(path_bed_ciri)

    path_bed_out = os.path.join(dir_par, "untouched_bsj.bed")

    _logger.debug("output bed file will in : {}".format(path_bed_out))

    r_args = ["Rscript", path_r_script, path_bed_ciri,
                      path_vis_list, path_bed_out]
    if exon_only:
        r_args.append("T")

    cmd_r = " ".join(r_args)

    _logger.debug("""external cmd for filter ciri bsj is:  {}""".format(cmd_r))
    pysrc.body.worker.run(cmd_r)

    return path_bed_out


def _find_that_r_script(basename_r):
    path_current_script = os.path.abspath(os.path.realpath(__file__))
    dir_r = os.path.dirname(os.path.dirname(
        os.path.dirname(path_current_script)))

    path_r_script = os.path.join(dir_r, "R", basename_r)

    return path_r_script


def exons_within_blank_region(bed_file_show_blank_region, genomic_annotation):
    bed_objs = load_bed_as_tuple(bed_file_show_blank_region)
    host_gene_of = dict([(str(x.name), str(x.itemRGB)) for x in bed_objs])

    exons_inside_blank_region = pysrc.file_format.bsj_gtf.intersect_region_genome_annotation(bed_file_show_blank_region,
                                                                                             genomic_annotation)

    for single_exon in exons_inside_blank_region:
        single_exon.set_gene_id(host_gene_of.get(
            single_exon.get_transcript_id()))

    return exons_inside_blank_region


def exons_directly_from_bed_file(bed_file):
    res = load_bed_as_tuple(bed_file)

    return [_bed2gtf(x) for x in res]


def load_bed_as_tuple(bed_file):
    import collections
    bed_obj = collections.namedtuple("bed", ["chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart",
                                             "thickEnd", "itemRGB", "blockCount", "blockSize", "blockStarts"])
    res = []
    with open(bed_file) as your_bed_file:
        for line in your_bed_file:
            parts = line.strip().split()
            len_should_pad = 12 - len(parts)
            if len_should_pad > 0:
                parts.extend([None] * len_should_pad)
            if len_should_pad < 0:
                parts = parts[:12]

            bed_obj_this_line = bed_obj._make(parts)
            res.append(bed_obj_this_line)
    return res


def _bed2gtf(bed):
    return pysrc.file_format.bsj_gtf.PredictedCircularRegion.generate_exon_for_circular_isoform(
        host_seqname=str(bed.chrom),
        start=int(bed.chromStart),
        end=int(bed.chromEnd),
        host_gene_id=str(bed.itemRGB),
        host_tran_id=str(bed.name),
        strand=bed.strand
    )


def summarize_circ_isoform_structure_marked_break(path_vis_list, genomic_gtf, summarized_gtf, tmp_dir=""):

    path_bed_circ_exon, path_bed_blank = translate_vis_list(
        path_vis_list, tmp_dir)

    import itertools
    whole_exon_generator = itertools.chain(exons_within_blank_region(path_bed_blank, genomic_gtf),
                                           exons_directly_from_bed_file(
                                               path_bed_circ_exon))

    pysrc.file_format.bsj_gtf.exons_to_gtf_file(
        whole_exon_generator, summarized_gtf)


def _rename_rebuild_fa(fa_file_in, use_suffix=True):
    import Bio.SeqIO
    fa_in = Bio.SeqIO.parse(fa_file_in, "fasta")
    __suffix_rebuild_fa = ".f"

    dict_fa_id = {}
    for seq in fa_in:
        short_id_isoform = seq.id.strip().split("#")[-1]
        if use_suffix:  # 19_11_16 using .f to mark rebuild fa.
            short_id_isoform = short_id_isoform + __suffix_rebuild_fa

        # in case there are multiple isoform under one BSJ
        if short_id_isoform not in dict_fa_id:
            dict_fa_id[short_id_isoform] = 1
        else:
            dict_fa_id[short_id_isoform] += 1
            short_id_isoform = ".".join(
                [short_id_isoform, str(dict_fa_id[short_id_isoform])])

        seq.id = short_id_isoform

        yield seq


def summarize_rebuild_fa(fa_in, fa_out, use_suffix=True):
    with open(fa_out, "w") as dump_it:
        dump_it.write("\n".join([">{fa_id}\n{fa_seq}".format(fa_id=seq.id.strip(),
                                                             fa_seq=seq.seq.strip())
                                 for seq in _rename_rebuild_fa(fa_in, use_suffix)]))
    return fa_out


if __name__ == "__main__":
    print(__doc__)
    print(opt_checker)
