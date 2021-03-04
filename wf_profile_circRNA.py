# !/usr/bin/env python
# -*- coding:utf-8 -*-
# author : zerodel
# Readme:
# changelog: 2020-01-09 : move from sailfish to salmon
#
import argparse
import copy
import os
import shutil

import pysrc.being.linc
import pysrc.body.cli_opts
import pysrc.body.config
import pysrc.body.logger
import pysrc.body.option_check
import pysrc.body.utilities
import pysrc.body.worker
import pysrc.file_format.bsj_gtf
import pysrc.file_format.ciri_as_to_gtf
import pysrc.file_format.ciri_entry
import pysrc.file_format.fa
import pysrc.sub_module.summary_quant
import pysrc.wrapper.gffread
import pysrc.wrapper.sailfish
import pysrc.wrapper.salmon
import pysrc.wrapper.ciri_full
from pysrc.body.cli_opts import catch_one

_SUB_DIR_PROFILE_RESULT = "profile_result"

_SUB_DIR_INDEX_FINAL = "index_final"

_OPT_CIRI_AS_OUTPUT_PREFIX = "ciri_as_prefix"

_OPT_VALUE_SAILFISH = "sailfish"

_OPT_KEY_QUANTIFIER = "quantifier"

_OPT_KEY_ADDITIONAL_CIRC_REF = "additional_circ_ref"

_OPT_KEY_ADDITIONAL_LINEAR_REF = "additional_linear_ref"

_OPT_KEY_ADDITIONAL_ANNOTATION = "additional_annotation"

_OPT_KEY_PRESERVE = "preserved_id_list"

_OPT_KEY_USE_LINC_EXPLICITLY = "flag_use_linc_explicitly"

_OPT_KEY_REJECT_LINEAR = "flag_reject_linear"

_QUANTIFIER_BACKEND_OF = {"sailfish": pysrc.wrapper.sailfish,
                          "salmon": pysrc.wrapper.salmon}

_OPT_KEY_DECOY_FILES = "--decoys"

__doc__ = '''
'''

__author__ = 'zerodel'

SECTION_PROFILE_CIRCULAR_RNA = "CIRC_PROFILE"

_logger = pysrc.body.logger.default_logger(SECTION_PROFILE_CIRCULAR_RNA)


def __cli_arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "cfg_file", help="file path to a configuration file of detection job")
    parser.add_argument("-l", "--log_file",
                        help="logging file path", default="")
    parser.add_argument(
        "-f", "--force", help="forced refresh", action="store_true")
    return parser


def _option_check_main_interface(opts=None):
    oc = pysrc.body.option_check.OptionChecker(
        opts, name=SECTION_PROFILE_CIRCULAR_RNA)

    oc.must_have(_OPT_KEY_QUANTIFIER, lambda x: x in ["sailfish", "salmon"],
                 KeyError(
                     "Error@circular_RNA_profiling: incorrect quantifier back-end"),
                 "the back end quantifier: sailfish or salmon")

    oc.one_and_only_one(["-g", "--genomic_seqs_fasta"], os.path.exists,
                        FileNotFoundError(
                            "Error@circular_RNA_profiling: incorrect genomic reference fa file"),
                        "path to genomic sequence fasta file")

    oc.one_and_only_one(["-a", "--annotation"], os.path.exists,
                        FileNotFoundError(
                            "Error@circular_RNA_profiling: incorrect genome annotation file"),
                        "path to gene annotation file, ie, .gtf or .gff files")

    oc.one_and_only_one(["-c", "--ciri_bsj", "--bed"], os.path.exists,
                        FileNotFoundError(
                            "Error@circular_RNA_profiling: can not find circular detection report"),
                        "path to  circRNA detection (file for ciri, folder for ciri-full) to specify circular RNA")

    oc.may_need(_OPT_CIRI_AS_OUTPUT_PREFIX, pysrc.body.cli_opts.is_suitable_path_with_prefix,
                FileNotFoundError(
                    "Error@circular_RNA_profiling: incorrect circular Alternative Splice file prefix"),
                "path prefix to CIRI-AS report of circular exons")

    # make sure the path should be a folder
    oc.may_need(_OPT_KEY_DECOY_FILES, lambda x: all([os.path.exists(y) for y in x.split()]),
                FileNotFoundError("ERROR@circular_RNA_profiling: unable to find decoy sequences"), "path to decoy sequences for salmon")

    oc.may_need(_OPT_KEY_PRESERVE, lambda x: all([os.path.exists(y) for y in x.split()]),
                FileNotFoundError("ERROR@circular_RNA_profiling: unable to find preserved id list"), "path to preserved id list ")

    oc.must_have("-o", pysrc.body.utilities.make_sure_there_is_a_folder,
                 FileNotFoundError(
                     "Error@circular_RNA_profiling: no place for output"),
                 "output folder that contains the index built by sailfish and quantification results")

    oc.may_need("--mll", lambda x: x.isdigit(),
                TypeError(
                    "Error@circular_RNA_profiling: mean library length should be a integer"),
                "mean library length, this option is to fix up the effective length.")

    oc.may_need("-k", lambda x: x.isdigit(),
                TypeError(
                    "Error@circular_RNA_profiling: k in k mer should be integer"),
                "k-mer size used by sailfish to built index. default is 21")

    oc.may_need("-1", os.path.exists,
                FileNotFoundError(
                    "Error@circular_RNA_profiling: no mate1 input seq file"),
                "path to raw pair-end reads, mate 1")
    oc.may_need("-2", os.path.exists,
                FileNotFoundError(
                    "Error@circular_RNA_profiling: no mate2 input seq file"),
                "path to raw pair-end reads, mate 2")
    oc.may_need("-r", os.path.exists,
                FileNotFoundError(
                    "Error@circular_RNA_profiling: no single end sequence input file"),
                "path to single-end raw sequencing reads file.")

    oc.may_need(_OPT_KEY_ADDITIONAL_CIRC_REF, os.path.exists,
                FileNotFoundError(
                    "Error@circular_RNA_profiling: additional circular reference file not exist"),
                "path to additional circular RNA reference file (.fa), ")

    oc.may_need(_OPT_KEY_ADDITIONAL_ANNOTATION, os.path.exists,
                FileNotFoundError(
                    "Error@circular_RNA_profiling: additional annotation not exist"),
                "path to additional circular RNA annotation file in gtf format")

    oc.may_need(_OPT_KEY_ADDITIONAL_LINEAR_REF, os.path.exists,
                FileNotFoundError(
                    "Error@circular_RNA_profiling: additional linear reference file not exist"),
                "path to additional linear RNA reference file(.fa)")

    oc.may_need(_OPT_KEY_USE_LINC_EXPLICITLY, lambda x: x in ("T", "F", "True", "False", ""),
                KeyError(
                    "Error@circular_RNA_profiling: incorrect flag to specify whether linc should be explicit"),
                """flag to specify whether linc RNA should be include in quantification result, commenting it out to 
                set False
                """)

    oc.may_need(_OPT_KEY_REJECT_LINEAR, lambda x: x in ("T", "F", "True", "False", ""),
                KeyError("""Error@circular_RNA_profiling: incorrect flag to specify whether index should reject 
                linear RNA"""),
                """flag to specify whether to reject linear RNA during quantification, for example for a RNase R  
                treated sample, comment it out to set False"""
                )

    oc.forbid_these_args("-h", "--help")
    return oc


# assign package-scale global objects
option_checker = _option_check_main_interface()
OPTION_CHECKERS = [option_checker]

# dummy implement
_seq_extractor = pysrc.wrapper.gffread
_gtf_operator = pysrc.file_format.bsj_gtf


def __determine_kmer_length(obj_circ_profile):
    if "-k" in obj_circ_profile and obj_circ_profile['-k']:
        kmer_length = int(obj_circ_profile["-k"])
    else:
        kmer_length = 21
    return kmer_length


def _filter_gtf(gtf_in, by_what, gtf_out):
    id_pool = []
    with open(by_what) as perserve_it:
        for line in perserve_it:
            id_pool.append(line.strip())

    lines_output = []
    with open(gtf_in) as input_gtf:
        for line in input_gtf:
            if line.startswith("#"):
                continue
            hit = any([x in line for x in id_pool])
            is_exon = line.strip().split()[2] == "exon"
            if hit and is_exon:
                lines_output.append(line)

    with open(gtf_out, "w") as output_gtf:
        output_gtf.writelines(lines_output)


def main(path_config, forced_refresh=False):
    whole_config = pysrc.body.config.config(path_config)
    circ_profile_config = _load_to_update_default_options(path_config)

    _logger.debug("profile config dict is : %s" % str(circ_profile_config))

    option_checker.check(copy.copy(circ_profile_config))  # check your options

    quantifier = _confirm_quantifier(circ_profile_config)

    genomic_annotation = catch_one(circ_profile_config, "-a", "--annotation")
    genome_fa = catch_one(circ_profile_config, "-g", "--genomic_seqs_fasta")
    circ_detection_report = catch_one(
        circ_profile_config, "-c", "--bed", "--ciri_bsj")

    output_path = catch_one(circ_profile_config, "-o")

    # additional options
    additional_circ_ref = circ_profile_config.get(_OPT_KEY_ADDITIONAL_CIRC_REF)
    additional_annotation = circ_profile_config.get(
        _OPT_KEY_ADDITIONAL_ANNOTATION)
    additional_linear_ref = circ_profile_config.get(
        _OPT_KEY_ADDITIONAL_LINEAR_REF)

    decoy_file_for_salmon = circ_profile_config.get(_OPT_KEY_DECOY_FILES)
    use_linc = _OPT_KEY_USE_LINC_EXPLICITLY in circ_profile_config
    reject_linear = _OPT_KEY_REJECT_LINEAR in circ_profile_config

    preserved_id_file = circ_profile_config.get(_OPT_KEY_PRESERVE)

    # assign file path
    spliced_linear_reference = os.path.join(output_path, "ref_linear.fa")
    circular_rna_gtf = os.path.join(output_path, "circ_only.gtf")
    circ_reference_seq = os.path.join(output_path, "circ_only.fa")

    linc_rna_gtf = os.path.join(output_path, "linc_only.gtf")
    linc_reference_seq = os.path.join(output_path, "linc_only.fa")

    k = __determine_kmer_length(circ_profile_config)

    # 1st, extract the linear sequence
    if not reject_linear:
        if not os.path.exists(spliced_linear_reference) or forced_refresh:
            _prepare_linear_transcriptome(
                genome_fa, genomic_annotation, spliced_linear_reference)
            _logger.debug("linear RNA from Genomic Sequence is included")

    # 2nd, get the circular RNA gtf sequences
    circular_detection_not_only_bsj = os.path.isdir(circ_detection_report)
    if not os.path.exists(circular_rna_gtf) or forced_refresh:
        _prepare_circular_rna_annotation(circ_detection_report, circular_rna_gtf, genomic_annotation,
                                         circular_detection_not_only_bsj)

    # 3rd, extracts circular RNA sequence
    _seq_extractor.do_extract_non_coding_transcript(gff=circular_rna_gtf,
                                                    path_ref_sequence_file=genome_fa,
                                                    output=circ_reference_seq)

    # 4th, do operations on circular RNA reference .lincRNA are treated as linear mRNA

    # mean of effective length use 150 as default
    mean_library_length = int(circ_profile_config["--mll"]) if "--mll" in circ_profile_config and \
                                                               circ_profile_config[{"--mll"}] else 150

    _add_adapter_k_mll(circ_reference_seq,
                       circ_reference_seq, k, mean_library_length)

    # ###### ========================================================
    # process additional reference . including linc and custom circular RNA 18-12-21

    lst_reference_fa = [circ_reference_seq] if reject_linear else [
        spliced_linear_reference, circ_reference_seq]
    lst_annotation = [circular_rna_gtf] if reject_linear else [
        genomic_annotation, circular_rna_gtf]

    if additional_linear_ref:
        lst_reference_fa.extend([single_fa for single_fa in additional_linear_ref.strip().split() if os.path.exists(
            single_fa)])

    if circular_detection_not_only_bsj:
        _logger.debug("transform ciri_full rebuild fa file")

        ciri_full_rebuild_fa_file = pysrc.wrapper.ciri_full.rebuild_fa_path_under(
            circ_detection_report)

        if ciri_full_rebuild_fa_file:
            _logger.debug(
                "ciri-full rebuild circular RNA file found in {}".format(ciri_full_rebuild_fa_file))
            ciri_full_rebuild_fa_with_short_name = pysrc.wrapper.ciri_full.summarize_rebuild_fa(
                ciri_full_rebuild_fa_file, os.path.join(output_path, "ciri_full_renamed.fa"))

            # output of this step.
            rebuild_fa_encoded = os.path.join(output_path, "rebuild_seq.fa")

            # need adapter-k

            _add_adapter_k_mll(ciri_full_rebuild_fa_with_short_name,
                               rebuild_fa_encoded, k, mean_library_length)

            lst_reference_fa.append(rebuild_fa_encoded)

        else:
            _logger.warning(
                " NO CIRI-FULL rebuild file found under {}".format(circ_detection_report))

    if additional_circ_ref:
        additional_circ_ref_decoded = os.path.join(
            output_path, "additional_circ_ref.fa")

        # here also need adapter-k

        _add_adapter_k_mll(additional_circ_ref,
                           additional_circ_ref_decoded, k, mean_library_length)

        lst_reference_fa.append(additional_circ_ref_decoded)

    if additional_annotation:
        lst_annotation.append(additional_annotation)

    # Wednesday, 5 April 2017: add same procedure for lincRNA
    if use_linc:
        if not os.path.exists(linc_rna_gtf) or forced_refresh:
            pysrc.being.linc.prepare_linc_annotation(original_gff=genomic_annotation,
                                                     target_linc_annotation=linc_rna_gtf)

        if not os.path.exists(linc_reference_seq) or forced_refresh:
            pysrc.being.linc.prepare_linc_transcriptome_seq(linc_annotation=linc_rna_gtf,
                                                            genomic_seq=genome_fa,
                                                            target_fa=linc_reference_seq)
        lst_reference_fa.append(linc_reference_seq)

    if preserved_id_file:
        # prepare gtf file
        preserved_annotation = os.path.join(output_path, "preserved.gtf")
        _logger.debug("extract annotation for preserved, and store in {}".format(
            preserved_annotation))
        _filter_gtf(genomic_annotation, preserved_id_file,
                    preserved_annotation)

        preserved_fa = os.path.join(output_path, "preserved.fa")
        pysrc.wrapper.gffread.do_extract_non_coding_transcript(
            preserved_annotation, genome_fa, preserved_fa)
        _logger.debug(
            "sequence for preserved will be stored in {}".format(preserved_fa))

        lst_annotation.append(preserved_annotation)
        lst_reference_fa.append(preserved_fa)

    # 5th , combined those fa files
    final_refer = os.path.join(output_path, "final.fa")
    # pysrc.body.utilities.do_merge_files(final_refer, lst_reference_fa)
    pysrc.file_format.fa.incremental_updating(final_refer, lst_reference_fa)

    # linc RNA is already in original gtf file
    final_annotation = os.path.join(output_path, "final.gtf")
    pysrc.body.utilities.do_merge_files(final_annotation, lst_annotation)

    # 6th , make index for quantifier

    path_to_quantifier_index = os.path.join(output_path, _SUB_DIR_INDEX_FINAL)

    index_parameters = {"sailfish_bin": whole_config["META"]["sailfish_bin"],
                        "--kmerSize": str(k),
                        "--transcripts": final_refer,
                        "--out": path_to_quantifier_index
                        } if quantifier is pysrc.wrapper.sailfish else {
        "salmon_bin": whole_config["META"]["salmon_bin"],
        "--kmerLen": str(k),
        "--transcripts": final_refer,
        "--index": path_to_quantifier_index,
    }

    if quantifier is pysrc.wrapper.salmon and decoy_file_for_salmon:
        index_parameters[_OPT_KEY_DECOY_FILES] = decoy_file_for_salmon

    quantifier.index(para_config=index_parameters)

    # 7th , do quantification!
    path_to_quantify_result = os.path.join(
        output_path, _SUB_DIR_PROFILE_RESULT)
    if not os.path.exists(path_to_quantify_result):
        os.mkdir(path_to_quantify_result)

    opts_quantifier = {"salmon_bin": whole_config["META"]["salmon_bin"]} if quantifier is pysrc.wrapper.salmon else {
        "sailfish_bin": whole_config["META"]["sailfish_bin"]}

    opts_quantifier["--index"] = path_to_quantifier_index
    if "-1" and "-2" in circ_profile_config:
        opts_quantifier["--mates1"] = circ_profile_config["-1"]
        opts_quantifier["--mates2"] = circ_profile_config["-2"]
        opts_quantifier["--libType"] = 'IU'
    elif "-r" in circ_profile_config:
        opts_quantifier["--unmatedReads"] = circ_profile_config["-r"]
        opts_quantifier["--libType"] = 'U'

    opts_quantifier["--output"] = path_to_quantify_result
    opts_quantifier["--geneMap"] = final_annotation

    # # on salmon's bias model
    if quantifier is pysrc.wrapper.salmon:
        opts_quantifier["--seqBias"] = None
        opts_quantifier["--gcBias"] = None
        opts_quantifier["--validateMappings"] = ""

    quantifier.quantify(para_config=opts_quantifier)  # do qunatification

    _logger.debug("All finished....... exiting ...")
    # pysrc.sub_module.summary_quant.aggregate_isoform_quantify_result(
    #     quant_sf=os.path.join(path_to_quantify_result, "quant.sf"),
    #     summarized_output=os.path.join(path_to_quantify_result, "summarized.quant"),
    #     gtf_annotation=final_annotation)

    _post_quantify(sf_in=os.path.join(opts_quantifier["--output"], "quant.sf"),
                   anno=opts_quantifier["geneMap"],
                   tab_out=os.path.join(opts_quantifier["--output"], "post_ratio.tab"))


def _add_adapter_k_mll(raw_circ_seq, decorated_seq, k, mean_library_length):
    seq_backup_dest = raw_circ_seq + ".raw_seq"
    shutil.copyfile(raw_circ_seq, seq_backup_dest)
    if os.path.exists(seq_backup_dest):
        _logger.debug("sequence file {src}  has been backed up in {dest}".format(src=raw_circ_seq,
                                                                                 dest=seq_backup_dest))
    else:
        _logger.error("backup failed for {}".format(raw_circ_seq))

    pysrc.file_format.fa.convert_all_entries_in_fasta(fa_in=raw_circ_seq,
                                                      fa_out=decorated_seq,
                                                      convert_fun=pysrc.file_format.fa.make_adapter(k))

    pysrc.file_format.fa.convert_all_entries_in_fasta(fa_in=decorated_seq,
                                                      fa_out=decorated_seq,
                                                      convert_fun=pysrc.file_format.fa.pad_for_effective_length(
                                                          mean_library_length))
    _logger.debug("{} has been decorated".format(raw_circ_seq))


def _prepare_circular_rna_annotation(circ_detection_report, circular_rna_gtf, genomic_annotation,
                                     detection_not_only_bsj):
    # this means use ciri-full as circRNA detection source.
    if detection_not_only_bsj:

        dir_par = os.path.dirname(circular_rna_gtf)

        _logger.debug(
            "building circRNA using comprehensive information (ciri-full) from :{}".join(circular_rna_gtf))

        ciri_full_list_file = pysrc.wrapper.ciri_full.vis_list_path_under(
            circ_detection_report)

        ciri_report_path = pysrc.wrapper.ciri_full.ciri_report_under(
            circ_detection_report)

        _logger.debug(
            "first, we need BSJ information in {}".format(ciri_report_path))
        bed_bsj_only = pysrc.wrapper.ciri_full.filter_out_un_touched_circular_rna(ciri_report_path,
                                                                                  ciri_full_list_file,
                                                                                  tmp_dir=dir_par,
                                                                                  exon_only=True)

        bsj_only_gtf = os.path.join(dir_par, "bsj_only.gtf")
        _logger.debug(
            "annotation for BSJ only will be put in : {}".format(bsj_only_gtf))
        _logger.debug(
            "genomic annotation source is from : {} . ".format(genomic_annotation))
        _gtf_operator.do_make_gtf_for_circular_prediction_greedy(
            bed_bsj_only, genomic_annotation, bsj_only_gtf)

        partial_gtf = os.path.join(dir_par, "partial_structure.gtf")
        _logger.debug(
            "circular RNA with partial structure information is in {}".join(partial_gtf))
        pysrc.wrapper.ciri_full.summarize_circ_isoform_structure_marked_break(ciri_full_list_file, genomic_annotation,
                                                                              partial_gtf, dir_par)

        pysrc.body.utilities.do_merge_files(
            circular_rna_gtf, [bsj_only_gtf, partial_gtf])
        _logger.debug("the final circular exclusive annotation file is : {} ".format(
            circular_rna_gtf))

    else:
        _logger.debug("building circRNA GTF using BSJ information only")
        _gtf_operator.do_make_gtf_for_circular_prediction_greedy(circular_candidate_regions=circ_detection_report,
                                                                 gff_db=genomic_annotation,
                                                                 output_gtf_path_name=circular_rna_gtf)


def _prepare_linear_transcriptome(genome_fa, genomic_annotation, spliced_linear_reference):
    _seq_extractor.do_extract_classic_message_transcript(gff=genomic_annotation,
                                                         path_ref_sequence_file=genome_fa,
                                                         output=spliced_linear_reference)


def _load_to_update_default_options(path_config):
    user_config = pysrc.body.config.config(
        path_config) if path_config else pysrc.body.config.load_default_value()

    if SECTION_PROFILE_CIRCULAR_RNA not in user_config:
        raise KeyError(
            "ERROR@Config_file: should have a section with the name :{}".format(SECTION_PROFILE_CIRCULAR_RNA))

    user_option_section = dict(user_config[SECTION_PROFILE_CIRCULAR_RNA])
    default_config = pysrc.body.config.load_default_value()

    if SECTION_PROFILE_CIRCULAR_RNA in default_config:
        default_option_section = dict(
            default_config[SECTION_PROFILE_CIRCULAR_RNA])
        default_option_section.update(user_option_section)
        user_option_section = default_option_section

    return user_option_section


def _confirm_quantifier(circ_profile_config):
    str_quantifier = circ_profile_config.get(
        _OPT_KEY_QUANTIFIER, _OPT_VALUE_SAILFISH)
    quantifier = _QUANTIFIER_BACKEND_OF[str_quantifier]
    _logger.debug("using %s as quantification backend" % str_quantifier)
    return quantifier


def _find__r_script(basename_r):
    path_current_script = os.path.abspath(os.path.realpath(__file__))
    dir_r = os.path.dirname(path_current_script)
    path_r_script = os.path.join(dir_r, "R", basename_r)
    return path_r_script


def _post_quantify(sf_in, anno, tab_out):
    # process the isoform level sf file ,and calculate the ratio for DE
    path_r_script = _find__r_script("post_ratio.R")
    # Rscript post_ratio.R gtf sf ratio_out
    
    if not os.path.exists(sf_in):
        _logger.error("NO sf_in file in {} ! unable to calculate the ratio".format(sf_in))
        return None
    if not os.path.exists(anno):
        _logger.error("NO genomic annotation  in {} ! unable to calculate the ratio".format(anno))
        return None

    cmd = " ".join(["Rscript", path_r_script, anno, sf_in, tab_out])

    _logger.debug("start post quantification step , calculate the ratio ")
    pysrc.body.worker.run(cmd)


if __name__ == "__main__":
    arg_parser = __cli_arg_parser()
    cl_args = arg_parser.parse_args()
    _logger = pysrc.body.logger.set_logger_file(_logger, cl_args.log_file)
    main(cl_args.cfg_file, forced_refresh=cl_args.force)
