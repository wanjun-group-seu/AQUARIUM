#! /usr/bin/env Rscript

# description of this script ----------------------------------------------
#' FUNCTIONALITY

#' USAGE

#' ARGUMENTS


# check dependency  and try install missing package -----------------------


# clean environment and load packages --------------------------------------------------------
base::rm(list = base::ls())

#' function definitions ------------------------------------
adhoc_load_list <- function(path_list) {
    col_name_list <- c(
        "image_id",
        "bsj",
        "chr",
        "start",
        "end",
        "exp",
        "isoform_num",
        "isoform_exp",
        "isoform_len",
        "isoform_state",
        "strand",
        "origin_gene",
        "isoform_circexon"
    )
    df_list <- read.table(
        path_list,
        header = F,
        sep = "\t",
        stringsAsFactors = F,
        row.names = NULL,
        skip = 1,
        col.names = col_name_list
    )
    colnames(df_list)[1] <- "id"
    return(df_list)
}

load_ciri_report <- function(path.ciri, ...) {
    tmp <-
        read.table(
            path.ciri,
            header = T,
            stringsAsFactors = F,
            sep = "\t",
            comment.char = ""
        )
    name.tmp <- names(tmp)
    name.tmp <- gsub("^#", "", name.tmp)
    names(tmp) <- name.tmp
    return(tmp)
}

# parse arguments ---------------------------------------------------------



cat(paste0("\n", "trim CIRI BSJ into bed file "))

given_args <- commandArgs(trailingOnly = T)

if (length(given_args) < 3) {
    stop(
        paste0(
            "essential argument lost, need 3 argument:\n  1. path_to_ciri_file \n 2. path_to_cirifull_list_file \n 3. where_to_put_bed_file \nbut only",
            length(given_args),
            " given"
        )
    )
}


ciri_report_path <- given_args[1]
list_path <- given_args[2]
output_bed_path <- given_args[3]

only_exon_circ_remain <- length(given_args) > 3

# loading data -----------------------------------------------------------
#
# ciri_report_path <- "./res/detection/ciri.report"
# list_path <- "./res/vis/stout.list"
#



if (!file.exists(ciri_report_path)) {
    stop(
        paste0("ERROR: NO ciri report file in ", ":", ciri_report_path),
        "\n"
    )
}

if (!file.exists(list_path)) {
    stop(paste0("ERROR: NO ciri-full list file in ", ":", list_path, "\n"))
}

raw_ciri <- load_ciri_report(ciri_report_path)
raw_list <- adhoc_load_list(list_path)

cat(paste0("\n", dim(raw_ciri)[1], " circRNA detected by BSJ information \n"))

cat(paste0("\n", dim(raw_list)[1], " circRNA with inner structure detected\n "))


filtered_ciri <-
    base::subset(raw_ciri, !raw_ciri$circRNA_ID %in% raw_list$bsj)

if (only_exon_circ_remain) {
    filtered_ciri <- base::subset(filtered_ciri, filtered_ciri$circRNA_type == "exon")
}

bed_out <- with(
    filtered_ciri,
    data.frame(
        "chrom" = chr,
        "chromStart" = circRNA_start,
        "chromEnd" = circRNA_end,
        "name" = circRNA_ID,
        "score" = X.junction_reads,
        "strand" = strand,
        "thickStart" = X.junction_reads,
        "thickEnd" = X.non_junction_reads,
        "itemRGB" = gsub(",$", "", gene_id)
    )
)

cat(paste0("\n", "there is ", dim(bed_out)[1], " untouched bsjs\n"))


write.table(
    bed_out,
    file = output_bed_path,
    col.names = F,
    row.names = F,
    quote = F,
    sep = "\t",
)

cat(paste0("\n", "output file in : ", output_bed_path, "\n"))