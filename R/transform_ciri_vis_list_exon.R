#! /usr/bin/env Rscript

# description of this script ----------------------------------------------
#' FUNCTIONALITY
#' transform ciri-vis *.list file into exon bed files

#' USAGE
#' Rscript this.script path_to_ciri_vis.list path_to_exon_bed path_to_empty_bed
#'
#' ARGUMENTS
#' path_to_ciri_vis.list : path to intermedian file under ciri-vis report folder
#' path_to_exon_bed : path to bed file contain KNOWN exons
#' path_to_empty_bed: path to bed file contain "empty" regions


# check dependency  and try install missing package -----------------------
#' use base R , no additional package needed .

# clean environment and load packages --------------------------------------------------------
base::rm(list = base::ls())

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
} # end of adhoc_load_list



#' count the occurence of each element in a vector
#'
#' @param x an vector ,
#' elements in x should be as.charact
#'
#' @return a vect showing how many time each element occurs
#' @export
#'
#' @examples
#' x <- c(1,1,2,3,4,5,5,4,5)
#' print(count_of_occur(x))
#' [1] 1 2 1 1 1 1 2 2 3
count_of_occur <- function(x) {
    tmpx <- seq_along(x)
    names(tmpx) <- as.character(x)
    count_x <- rep(NA, length(x))
    count_now <- 1
    
    while (any(is.na(count_x))) {
        if (count_now > length(x)) {
            stop("error here, count is larger than whole length")
        }
        # check duplicated
        is_dup <- duplicated(names(tmpx))
        #mark count
        count_x[tmpx[!is_dup]] <- count_now
        # update tmpx
        tmpx <- tmpx[is_dup]
        # increase count num
        count_now = count_now + 1
    }
    return(count_x)
}  # end of count_of_occur


#' mark which the current element is duplicated
#'
#' @param x input vector
#'
#' @return bool vector
#' T: current element occurs multiple times
#' F: current element occurs only once
#' @export
#'
#' @examples
#' x <- c(1,2,3,3,3)
#' > is_duplicate(x)
#' [1] FALSE FALSE  TRUE  TRUE  TRUE
#'
is_duplicate <- function(x) {
    x %in% base::unique(x[base::duplicated(x)])
    # count_x <- table(x)
    # has_dup <- names(count_x)[count_x > 1]
    # return(x %in% has_dup)
}

start_of_region <- function(s) {
    as.numeric(base::regmatches(s, base::regexpr("^\\d+", s)))
}

end_of_region <- function(s) {
    as.numeric(base::regmatches(s, base::regexpr("\\d+$", s)))
}

start_of_bsj <- function(s) {
    as.numeric(base::regmatches(s, base::regexpr("(?<=\\:)\\d+(?=|)", s, perl = T)))
}

end_of_bsj <- function(s) {
    as.numeric(base::regmatches(s, base::regexpr("(?<=\\|)\\d+(?=\\.*)", s, perl = T)))
}

chr_of_bsj <- function(s) {
    base::regmatches(s, base::regexpr("^.*(?=\\:)", s, perl = T))
}

count_suffix_of <- function(x) {
    has_dup <- is_duplicate(x)
    count <- count_of_occur(x)
    suffix <- ifelse(has_dup, paste0(".", count), "")
    return(suffix)
}

adhoc_unique_id_of <- function(x) {
    suffix <- count_suffix_of(x)
    uid <- base::paste0(x, suffix)
    return(uid)
} # end of adhoc_unique_id_of

add_unique_id <- function(df_lst) {
    isoform_state_suffix <-
        ifelse(df_lst$isoform_state == "Full", ".f", ".b")
    count_suffix <- count_suffix_of(df_lst$bsj)
    
    df_lst$uid <-
        base::paste0(df_lst$bsj, isoform_state_suffix, count_suffix)
    
    return(df_lst)
}

adhoc_lst_to_df <- function(bsj, circexons) {
    lst_circexon <- strsplit(circexons, split = ",")
    names(lst_circexon) <- bsj
    
    vec_circexon <- NULL
    for (x in names(lst_circexon)) {
        region_exon <- lst_circexon[[x]]
        names(region_exon) <- base::rep(x, length(region_exon))
        vec_circexon <- c(vec_circexon, region_exon)
    }
    
    df_circexon <-
        data.frame("bsj" = names(vec_circexon), "circexon" = vec_circexon)
    
    df_circexon$start <- start_of_region(df_circexon$circexon)
    df_circexon$end <- end_of_region(df_circexon$circexon)
    df_circexon$start_bsj <- start_of_bsj(df_circexon$bsj)
    df_circexon$end_bsj <- end_of_bsj(df_circexon$bsj)
    return(df_circexon)
} # end of adhoc_lst_to_df

adhoc_deco_lst_df_with_anno_info <-
    function(df_raw_region, df_list) {
        df_raw_region$strand <-
            df_list$strand[match(df_raw_region$bsj, df_list$uid)]
        df_raw_region$gene <-
            df_list$origin_gene[match(df_raw_region$bsj, df_list$uid)]
        df_raw_region$chr <- chr_of_bsj(df_raw_region$bsj)
        return(df_raw_region)
    } # end of adhoc_deco_lst_df_with_anno_info


confirm_empty_region <- function(tmp_df) {
    index_empty <- which(tmp_df$start == 0 & tmp_df$end == 0)
    
    if (1 == index_empty) {
        cat(paste0("\n 0-0 is in the head of exon region: ", tmp_df$bsj[1]))
        return(c(tmp_df$start_bsj[1]), tmp_df$start[2])
    }
    if (dim(tmp_df)[1]  == index_empty) {
        cat(paste0("\n 0-0 is in the end of exon region: ", tmp_df$bsj[1]))
        return(c(tmp_df$end[index_empty - 1], tmp_df$end_bsj[index_empty]))
    }
    return(c(tmp_df$end[index_empty - 1], tmp_df$start[index_empty + 1]))
} # end of confirm_empty_region

summarize_region_empty <- function(df_circexon) {
    cat("here, we assume all exon are arranged by position 5' -> 3'")
    
    lst_region <-
        lapply(split(df_circexon, df_circexon$bsj),
               confirm_empty_region)
    
    df_res <- data.frame(
        "chr" = chr_of_bsj(names(lst_region)),
        "start" = sapply(lst_region, function(x) {
            x[1]
        }),
        "end" = sapply(lst_region, function(x) {
            x[2]
        }),
        "bsj" = names(lst_region),
        stringsAsFactors = F
    )
    
    num_bsj_abnormal <- sum(!(df_res$start < df_res$end))
    
    if (num_bsj_abnormal > 0) {
        cat(
            paste0(
                "\nthere is " ,
                num_bsj_abnormal,
                " empty region with length less than 0, should be removed  "
            )
        )
    }
    
    return(df_res)
} # end of summarize_region_empty

# need a filter process
filter_region_empty <- function(df_region_empty) {
    subset(df_region_empty,
           df_region_empty$start < df_region_empty$end)
}

# parse arguments ---------------------------------------------------------


given_args <- commandArgs(trailingOnly = T)
if (length(given_args) < 3) {
    stop(
        paste0(
            "essential argument lost, need 3 argument:\n  1. path_to_raw_list_file\n 2. where_to_put_circexon_bed_file\n 3. where_to_put_empty_bed_file \nbut only ",
            length(given_args),
            " given"
        )
    )
}


raw_lst_path <- given_args[1]
exon_bed_path <- given_args[2]
empty_bed_path <- given_args[3]


# loading files -----------------------------------------------------------

cat(paste0("reading circ-full list file from ", raw_lst_path, "\n"))

df_lst <- adhoc_load_list(raw_lst_path)

cat(paste0(
    "\n there is ",
    sum(df_lst$isoform_state == "Full"),
    " complete sequence in ",
    length(df_lst$isoform_state),
    " BSJs\n"
))

# parsing ciri_full list file  ====

#df_lst$uid <- adhoc_unique_id_of(df_lst$bsj)
df_lst <- add_unique_id(df_lst)

df_partial_isoform <-
    base::subset(df_lst, df_lst$isoform_state == "Break")


if (dim(df_partial_isoform)[1] < 1) {
    file.create(exon_bed_path)
    cat("\n all your bsj sequence is intact, nothing to do here, quitting \n")
    
} else{
    # [2019_10_15 18:28] we should make unique id before division
    #  df_partial_isoform$uid <- adhoc_unique_id_of(df_partial_isoform$bsj)
    cat(paste0(
        "\n",
        "we need to serparate the verified exons and  0-0 empty region"
    ))
    
    df_raw_region <-
        adhoc_lst_to_df(df_partial_isoform$uid,
                        df_partial_isoform$isoform_circexon)
    
    df_raw_region <-
        adhoc_deco_lst_df_with_anno_info(df_raw_region, df_partial_isoform)
    
    df_region_circexon <-
        base::subset(df_raw_region,
                     !(df_raw_region$start == 0 |
                           df_raw_region$end == 0))
    
    df_region_empty <- summarize_region_empty(df_raw_region)
    df_region_empty <-
        adhoc_deco_lst_df_with_anno_info(df_region_empty, df_partial_isoform)
    df_region_empty$start_bsj <-
        df_region_circexon$start_bsj[base::match(df_region_empty$bsj, df_region_circexon$bsj)]
    df_region_empty$end_bsj <-
        df_region_circexon$end_bsj[base::match(df_region_empty$bsj, df_region_circexon$bsj)]
    
    bed_export_circexon <- with(
        df_region_circexon,
        data.frame(
            "chr" = chr,
            "start" = start,
            "end" = end,
            "name" = bsj,
            "score" = end - start + 1,
            "strand" = strand,
            thickStart = start_bsj,
            thickEnd = end_bsj,
            itemRGB = gene
        )
    )
    
    
    
    bed_export_empty <- with(
        df_region_empty,
        data.frame(
            "chr" = chr,
            "start" = start,
            "end" = end,
            "name" = bsj,
            "score" = end - start + 1,
            "strand" = strand,
            thickStart = start_bsj,
            thickEnd = end_bsj,
            itemRGB = gene
        )
    )
    
    
    write.table(
        bed_export_circexon,
        file = exon_bed_path,
        row.names = F,
        col.names = F,
        sep = "\t",
        quote = F
    )
    cat(paste0("\n circexon bed file created under ", exon_bed_path, "   \n"))
    
    write.table(
        bed_export_empty,
        file = empty_bed_path,
        row.names = F,
        col.names = F,
        sep = "\t",
        quote = F
    )
    cat(paste0("\n", "empty region bed file located in : ", exon_bed_path, " \n"))
    
}

cat(paste0("\n", "ciri-vis list file dissert, R part is OK\n"))
