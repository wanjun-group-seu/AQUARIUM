#! /usr/bin/env Rscript

base::rm(list = base::ls())

library(dplyr)
library(magrittr)
library(rtracklayer)
library(GenomicRanges)

# functions ---------------------------------------------------------------

#' Load CIRI report file
#'
#' @param path.ciri path to your ciri report file
#' @param ...  additional parameters
#'
#' @return a dataframe contains ciri report information
#' @export
#'
#' @examples
load_ciri_report <- function(path.ciri, ...) {
  if (require(data.table)) {
    tmp <-
      data.table::fread(path.ciri,
        header = T,
        stringsAsFactors = F,
        ...
      )
  } else {
    tmp <-
      utils::read.table(path.ciri,
        header = T,
        stringsAsFactors = F,
        ...
      )
  }

  name.tmp <- names(tmp)
  name.tmp <- gsub("^#", "", name.tmp)
  names(tmp) <- name.tmp
  return(tmp)
} # end of function: load_ciri_report


#' read gtf as granges
#'
#' @param path_gtf path to a gtf file
#'
#' @return a granges object contains the information of gtf file
#' @export
#'
#' @examples
load_gtf_as_gr <- function(path_gtf, ...) {
  linear_anno <-
    rtracklayer::import(path_gtf, ...)

  GenomicRanges::GRanges(
    as.character(seqnames(linear_anno)),
    ranges = ranges(linear_anno),
    strand = strand(linear_anno),
    mcols(linear_anno)
  )
} # end of function: load_gtf_as_gr


load_untouched_bed_as_gr <- function(path_untouched_bed) {
  x <- data.table::fread(path_untouched_bed)
  names(x) <- c(
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "itemRGB"
  )

  GenomicRanges::makeGRangesFromDataFrame(
    x,
    seqnames.field = "chrom",
    start.field = "chromStart",
    end.field = "chromEnd",
    strand.field = "strand",
    keep.extra.columns = T
  )
} # end of function: load_untouched_bed_as_gr

#' 总体中哪些门类在采样中被全取出来了
#'
#' @param sample_part  采样部分
#' @param population  总体数据
#'
#' @return 被采样取干净的门类名称
#' @export
#'
#' @examples
which_sample_is_total_picked_out <-
  function(sample_part, population) {
    df_sample_cnt <- as.data.frame(table(sample_part))
    names(df_sample_cnt) <- c("name", "cnt_sample")

    df_population_cnt <- as.data.frame(table(population))
    names(df_population_cnt) <- c("name", "cnt_population")

    df_which_is_total_picked_out <-
      dplyr::left_join(df_sample_cnt, df_population_cnt, by = c("name" = "name")) %>% dplyr::filter(cnt_sample == cnt_population)

    as.character(df_which_is_total_picked_out$name)
  } # end of function: which_sample_is_total_picked_out


#' print something into log
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
say_log <- function(...) {
  cat("\n")
  cat(paste(..., sep = " ", collapse = "\n "))
  cat("\n")
} # end of function: say_log


helper_pick_first_occur <- function(some_gr, name) {
  some_gr[match(name, some_gr$name)]
}



helper_pick_by_name <- function(some_gr, name) {
  some_gr[some_gr$name == name]
}
# load the data ----------------------------------------------------------


cat(paste0("\n", "Start finding circRNA that hit known database"))

given_args <- commandArgs(trailingOnly = T)

if (length(given_args) < 4) {
  stop(
    paste0(
      "essential argument lost, need 4 argument:\n  1. path_to_bsj \n 2. path_to_partial_rebuild \n 3. annotation \n 4. path_to_id_hit_annotation \nbut only",
      length(given_args),
      " given"
    )
  )
}


path_untouch_bed <- given_args[1]
path_partial_bed <- given_args[2]
path_anno <- given_args[3]
path_output <- given_args[4]


say_log("start loading known database GTF file ....")
known_gr <- load_gtf_as_gr(path_anno, format = "GFF")

if (file.exists(path_partial_bed)) {
  partial_gr <- load_untouched_bed_as_gr(path_partial_bed)
} else {
  say_log("no partial circRNA ....")
  partial_gr <- NULL
}

if (file.exists(path_untouch_bed)) {
  if (length(grep(".bed$", path_untouch_bed))>0) {
    say_log("loading untouched bsj...:", path_untouch_bed)
    untouched_gr <- load_untouched_bed_as_gr(path_untouch_bed)
  } else {
    say_log("loading untouched bsj from ciri report", path_untouch_bed)
    ciri_report <- load_ciri_report(path_untouch_bed)
    untouched_gr <-
      GenomicRanges::GRanges(
        seqnames = ciri_report$chr,
        ranges = IRanges(
          start = ciri_report$circRNA_start,
          end = ciri_report$circRNA_end
        ),
        strand = ciri_report$strand
      )

    untouched_gr$name <- ciri_report$circRNA_ID
  }
} else {
  untouched_gr <- NULL
}


# check partial part ------------------------------------------------------

check_partial_part <- function(partial_gr, known_gr, verbose=F) {
  if (is.null(partial_gr)) {
    say_log("no partial circRNA, so no circRNA hit database")
    return(NULL)
  }


  partial_that_hits <- partial_gr[partial_gr %in% known_gr]

  partial_id_all_exon_hits <-
    which_sample_is_total_picked_out(partial_that_hits$name, partial_gr$name)

  partial_all_hits <-
    partial_that_hits[partial_that_hits$name %in% partial_id_all_exon_hits]

  exon_hits_pair <- findMatches(partial_all_hits, known_gr)

  df_exon_hits_pair <-
    data.frame(
      "name_from" = partial_all_hits$name[exon_hits_pair@from],
      "name_to" = known_gr$name[exon_hits_pair@to],
      "from" = exon_hits_pair@from,
      "to" = exon_hits_pair@to
    )


  circ_need_check <- unique(partial_all_hits$name)
  tx_cover_this_circ <- sapply(circ_need_check, function(x) {
    df_tmp <- df_exon_hits_pair %>% dplyr::filter(name_from == x)
    count_tmp <- table(df_tmp$name_to, df_tmp$from) %>% as.matrix()

    is_tx_cover_circ <- apply(count_tmp, 1, function(x) {
      all(x > 0)
    })

    if (any(is_tx_cover_circ)) {
      rownames(count_tmp)[which(is_tx_cover_circ)[1]]
    } else{
      NA
    }
  })

  if (verbose) {
    data.frame("circ" = circ_need_check,
               "known_tx" = tx_cover_this_circ)

  } else{
    return(circ_need_check[!is.na(tx_cover_this_circ)])
  }
} # end of function: check_partial_part

partial_circ_hit_known <- check_partial_part(partial_gr, known_gr)
say_log("partial circRNA checked , total ", length(partial_circ_hit_known), "found ")

# check only bsj part -----------------------------------------------------

check_untouch_bsj <- function(untouched_gr, known_gr, verbose = F) {
  if (is.null(untouched_gr)) {
    say_log("no untouched bsj ....")
    return(NULL)
  }

  untouch5 <-
    resize(untouched_gr,
           width = 2,
           fix = "end",
           use.names = T)
  untouch5$side <- "p5"

  untouch3 <-
    resize(untouched_gr,
           width = 2,
           fix = "start",
           use.names = T)
  untouch3$side <- "p3"

  untouch_end <- c(untouch5, untouch3)

  hit_untouch_end <- findOverlaps(untouch_end, known_gr)

  df_hit_untouch <-
    data.frame(
      "name_from" = untouch_end$name[hit_untouch_end@from],
      "side_from" = untouch_end$side[hit_untouch_end@from],
      "name_to" = known_gr$transcript_id[hit_untouch_end@to],
      "from" = hit_untouch_end@from,
      "to" = hit_untouch_end@to
    )


#' 移除只命中单边的.
  df_side_align <-
    df_hit_untouch %>%
    dplyr::group_by(name_from, name_to) %>%
    dplyr::summarise(both_side = length(unique(side_from)))

  df_have_both_side <-
    df_side_align %>% dplyr::filter(both_side == 2)

  df_hit_untouch_both_side <-
    df_hit_untouch %>% dplyr::filter(
      name_from %in% df_have_both_side$name_from,
      name_to %in% df_have_both_side$name_to
    )


  untouch_hit_both <-
    untouched_gr[untouched_gr$name %in% df_hit_untouch_both_side$name_from]

  #' 检查数据库里面的isoform 的起止

  known_gr_hits_both_side <-
    known_gr[known_gr$transcript_id %in% df_hit_untouch_both_side$name_to]

  lst_known_gr <-
    split(known_gr_hits_both_side, f = known_gr_hits_both_side$transcript_id)

  range_gr <- range(GRangesList(lst_known_gr)) %>% stack()

  match_of_bsj <- match(untouch_hit_both, range_gr)

  if (verbose) {
    df_res <- data.frame("circ" = untouch_hit_both$name,
                         "known_tx" = match_of_bsj) %>% dplyr::filter(!is.na(known_tx))    %>% dplyr::mutate(known_tx = as.character(range_gr$name)[known_tx])
    return(df_res)
  } else{
    return(untouch_hit_both$name[!is.na(match_of_bsj)])
  }
} # end of function: check_untouch_bsj

untouch_circ_hit_known <- check_untouch_bsj(untouched_gr, known_gr)
say_log("BSJ checked, total ", length(untouch_circ_hit_known), "found.")

# combine two parts -------------------------------------------------------

id_hit_known <- c(partial_circ_hit_known, untouch_circ_hit_known)

writeLines(id_hit_known, path_output)

say_log("circRNA ID exported to ", path_output)
say_log("total ", length(id_hit_known), " inside. ")