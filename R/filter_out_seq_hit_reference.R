#! /usr/bin/env Rscript

base::rm(list = base::ls())


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
    df_sample_cnt <- table(sample_part) %>% as.data.frame()
    names(df_sample_cnt) <- c("name", "cnt_sample")

    df_population_cnt <- table(population) %>% as.data.frame()
    names(df_population_cnt) <- c("name", "cnt_population")

    df_which_is_total_picked_out <-
      dplyr::left_join(df_sample_cnt, df_population_cnt, by = c("name" = "name")) %>% dplyr::filter(cnt_sample == cnt_population)

    as.character(df_which_is_total_picked_out$name)
  } # end of function: which_sample_is_total_picked_out

# load the data ----------------------------------------------------------


cat(paste0("\n", "Find circRNA already in database"))

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

stopifnot()

known_gr <- load_gtf_as_gr(path_anno, format = "GFF3")

if (file.exists(path_partial_bed)) {
  partial_gr <- load_untouched_bed_as_gr(path_partial_bed)
} else {
  partial_gr <- NULL
}

if (file.exists(path_untouch_bed)) {
  if (endsWith(path_untouch_bed, ".bed")) {
    untouched_gr <- load_untouched_bed_as_gr(path_untouch_bed)
  } else {
    cat("\n loading untouched bsj from ciri report \n")
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

check_partial_part <- function(partial_gr, known_gr) {
  if (is.null(partial_gr)) {
    return(NULL)
  }

  hitpart <- partial_gr[partial_gr %in% known_gr]
  bsj_all_exon_hit <-
    which_sample_is_total_picked_out(hitpart$name, partial_gr$name)

  hitpart_likely <- hitpart[hitpart$name %in% bsj_all_exon_hit]

  hits_pair <- findMatches(hitpart_likely, known_gr)

  df_hits_pair <-
    data.frame(
      "name_from" = hitpart_likely$name[hits_pair@from],
      "name_to" = known_gr$name[hits_pair@to],
      "from" = hits_pair@from,
      "to" = hits_pair@to
    )

  circ_need_check <- unique(hitpart_likely$name)

  #' 如果所有exon都在某个已知的isoform中, 那么就是命中了
  circ_all_in_known_isoform <- sapply(circ_need_check, function(x) {
    df_tmp <- df_hits_pair %>% dplyr::filter(name_from == x)
    count_tmp <- table(df_tmp$name_to, df_tmp$from) %>% as.matrix()
    any(apply(count_tmp, 1, function(x) {
      all(x > 0)
    }))
  })

  circ_need_check[circ_all_in_known_isoform]
}
partial_circ_hit_known <- check_partial_part(partial_gr, known_gr)
# check only bsj part -----------------------------------------------------

check_untouch_bsj <- function(untouched_gr, known_gr) {
  if (is.null(untouched_gr)) {
    return(NULL)
  }

  untouch5 <-
    resize(untouched_gr,
      width = 2,
      fix = "end",
      use.names = T
    )
  untouch5$side <- "p5"

  untouch3 <-
    resize(untouched_gr,
      width = 2,
      fix = "start",
      use.names = T
    )
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
  #' 移除只命中一边的.
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

  untouch_hit_both$name[!is.na(match_of_bsj)]
}

untouch_circ_hit_known <- check_untouch_bsj(untouched_gr, known_gr)


# combine two parts -------------------------------------------------------

id_hit_known <- c(partial_circ_hit_known, untouch_circ_hit_known)

writeLines(id_hit_known, path_output)
