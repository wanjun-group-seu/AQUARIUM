#! /usr/bin/env Rscript

base::rm(list = base::ls())


# process the arguments
args_given <- commandArgs(trailingOnly = T)

if (length(args_given) < 3) {
    cat("\nUsage: \n
Rscript post_ratio.R gtf sf ratio_out\n

gtf: path to gtf file \n
sf: path to quant.sf file \n
ratio_out: path to circRNA ratio output \n
    ")

    quit(save = "no", status = 1, runLast = F)
}



path_gtf <- args_given[1]
path_sf <- args_given[2]
path_ratio <- args_given[3]



# if(FLAG_TEST){
# path_gtf <-  "./out/quant2/final.gtf"
# path_sf <-  "./out/quant2/profile_result/quant.sf"
# path_ratio <-  "./out/some_ratio.tab"

# }

# load librarys -----------------------------------------------------------


library(tidyverse)
library(rtracklayer)
library(GenomicFeatures)


#' Load Quant.sf file of salmon/sailfish
#'
#' @param path_to_sf path to quant.sf file
#' @param sample.id (optional) sample name
#'
#' @return data.frame of a single quant.sf file,
#' @export
#'
#' @examples
load_quant <- function(path_to_sf, sample.id = NULL) {
    if (require(data.table)) {
        x <- data.table::fread(path_to_sf, stringsAsFactors = F)
    } else {
        x <- read.table(path_to_sf, stringsAsFactors = F)
    }

    if (!is.null(sample.id)) {
        x$sample_id <- sample.id
    }
    if (require(tibble)) {
        x <- tibble::as_tibble(x)
    }
    return(x)
} # end of function: load_quant

#' 去掉带尾巴的id, 换成ciri风格的id
#'
#' @param bsj_with_tail
#'
#' @return
#' @export
#'
#' @examples
extract_bsj <- function(bsj_with_tail) {
    sapply(bsj_with_tail, function(x) {
        strsplit(x, "\\.")[[1]][[1]]
    })
} # end of function: extract_bsj

#' 把ciri风格的bsj id 转换成一个数据框
#'  seq start end bsj
#' @param bsj_raw
#' @param cut_tail
#'
#' @return
#' @export
#'
#' @examples
bsj2df <- function(bsj_raw, cut_tail = F) {
    if (cut_tail) {
        bsj <- extract_bsj(bsj_raw)
    } else {
        bsj <- bsj_raw
    }

    df_gr <-
        as.data.frame(t(sapply(bsj, function(x) {
            strsplit(gsub("\\|", ":", x), ":")[[1]]
        })))
    colnames(df_gr) <- c("seq", "start", "end")
    df_gr$bsj <- bsj_raw
    df_gr
} # end of function: bsj2df


# end of functions start process ------------------------------------------



gtf <- rtracklayer::import(path_gtf, format = "GFF")
df_quant <- load_quant(path_sf)


# todo: 需要有个方法来得到线性RNA与环状RNA的区分
non_bsj_rna <- stringr::str_subset(df_quant$Name, ":\\d", negate = T)


#' 筛出 非环状RNA的部分.
gtf_exon_expr <-
    gtf[gtf$transcript_id %in% non_bsj_rna & !is.na(gtf$exon_id)]

#' 找到环状RNA部分
bsj_rna <- stringr::str_subset(df_quant$Name, ":\\d")

# 对每个bsj找到一个起止点.
df_bsj <- bsj2df(bsj_rna, cut_tail = T) %>% dplyr::distinct()

gr_bsj_5 <-
    df_bsj %>%
    dplyr::mutate(end = start) %>%
    dplyr::distinct() %>%
    GenomicRanges::makeGRangesFromDataFrame(
        seqnames.field = "seq",
        start.field = "start",
        end.field = "end",
        keep.extra.columns = T
    )

gr_bsj_3 <-
    df_bsj %>%
    dplyr::mutate(start = end) %>%
    dplyr::distinct() %>%
    GenomicRanges::makeGRangesFromDataFrame(
        seqnames.field = "seq",
        start.field = "start",
        end.field = "end",
        keep.extra.columns = T
    )


#' 建立环状RNA首尾两个位置的1nt长的区域

#' 5-primer
hit5 <- GenomicRanges::findOverlaps(gr_bsj_5, gtf_exon_expr)

linear_exon_hit5 <- gtf_exon_expr[hit5@to]

#' 3-primer
hit3 <- GenomicRanges::findOverlaps(gr_bsj_3, gtf_exon_expr)
linear_exon_hit3 <- gtf_exon_expr[hit3@to]


tx_id_needed <-
    gtf_exon_expr[c(hit3@to, hit5@to) %>%
        unique()]$transcript_id %>% unique()

exon_tx_need <- gtf_exon_expr[gtf_exon_expr$transcript_id %in% tx_id_needed]

gr_tx_needed <- split(exon_tx_need, exon_tx_need$transcript_id)


# calculate the start and end point of transcript
start_of_tx <- min(start(gr_tx_needed))
end_of_tx <- max(end(gr_tx_needed))

#' 分别找到与环状RNA首尾相交的线性RNA.
#' 这里考虑了一些特殊情况, 就是环状RNA是第一个或者最后exon.
#' 5-primer端
df_hit_5 <- data.frame(
    "bsj" = gr_bsj_5$bsj[hit5@from],
    "tx" = gtf_exon_expr$transcript_id[hit5@to],
    "primer5" = start(gr_bsj_5)[hit5@from],
    "loc" = "primer5"
) %>%
    dplyr::mutate(
        "start_tx" = start_of_tx[tx],
        "is_overlap" = start_tx < primer5,
        "is_flush_with" = start_tx == primer5,
        "weired" = start_tx > primer5
    )


df_hit_3 <- data.frame(
    "bsj" = gr_bsj_3$bsj[hit3@from],
    "tx" = gtf_exon_expr$transcript_id[hit3@to],
    "primer3" = end(gr_bsj_3)[hit3@from],
    "loc" = "primer3"
) %>%
    dplyr::mutate(
        "end_tx" = end_of_tx[tx],
        "is_overlap" = end_tx > primer3,
        "is_flush_with" = end_tx == primer3,
        "weired" = end_tx < primer3
    )


#' 对于每个BSJ 计算ratio

df_quant_slim <-
    df_quant %>%
    dplyr::select(Name, TPM) %>%
    dplyr::rename("expr" = "TPM")

df_hit_tx_bsj <- dplyr::bind_rows(df_hit_5, df_hit_3) %>%
    dplyr::select(bsj, tx, loc)

df_hit_tx_bsj <- df_hit_tx_bsj %>%
    dplyr::left_join(df_quant_slim, by = c("bsj" = "Name")) %>%
    dplyr::rename("circ" = "expr")

df_hit_tx_bsj <- df_hit_tx_bsj %>%
    dplyr::left_join(df_quant_slim, by = c("tx" = "Name")) %>%
    dplyr::rename("linear" = "expr")

# # 这里可以用来修改数据来测试.
# if (FLAG_TEST) {
#   df_hit_tx_bsj$circ <- 0
# }

df_ratio_hit_tx_bsj <- df_hit_tx_bsj %>%
    dplyr::group_by(bsj) %>%
    dplyr::summarise(
        circ_final = mean(circ),
        linear_final = sum(linear) / 2,
        ratio = ifelse(circ_final + linear_final == 0,
            NA,
            circ_final / (linear_final + circ_final)
        )
    )

#' 输出相关的ratio到文件

write.table(df_ratio_hit_tx_bsj %>% dplyr::select(bsj, ratio),
    file = path_ratio,
    row.names = F, quote = F
)
