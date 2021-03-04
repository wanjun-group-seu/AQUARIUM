#! /usr/bin/env Rscript

# description of this script ----------------------------------------------
#' FUNCTIONALITY

#' USAGE

#' ARGUMENTS


# check dependency  and try install missing package -----------------------
CheckDependencies <- function(package_list) {
    for (p in package_list) {
        if (!suppressWarnings(suppressMessages(require(
            p,
            character.only = TRUE,
            quietly = TRUE,
            warn.conflicts = FALSE
        )))) {
            #' try use cran as package source.
            install.packages(p, repos = "http://cran.r-project.org")
            suppressWarnings(suppressMessages(library(
                p,
                character.only = TRUE,
                quietly = TRUE,
                warn.conflicts = FALSE
            )))
        }
    }
}


CheckDependencies(c("DESeq2", "tidyverse", "optparse"))


# clean enviroment and load packages --------------------------------------------------------
base::rm(list = base::ls())

source('./helper_functions.R')

# parse arguments ---------------------------------------------------------

adhocParseArgs <- function(verbose.print = T) {
    library(optparse)
    option_list <- list(
        make_option(
            c("-i", "--input"),
            type = "character",
            default = NULL,
            help = "directory of your quantification files, we assume each sample has a sub-directory under this folder"
        ),
        make_option(
            c("-o", "--output"),
            type = "character",
            default = NULL,
            help = "output directory or prefix, "
        ),
        make_option(
            c("-d", "--design"),
            type = "character",
            default = NULL,
            help = "design of sample series"
        ),
        make_option(
            c("-a", "--anno"),
            type = "character",
            default = NULL,
            help = "host information for each transcript"
        )
    )
    
    opts <- parse_args(OptionParser(option_list = option_list))
    
    if (verbose.print) {
        print(paste("The input file is ", opts$input,  sep = ""))
        print(paste("The output file prefix is ", opts$output, sep = ""))
    }
    return(opts)
} # end of adhocParseArgs


#' parse and assign arguments
args <- adhocParseArgs()
path.quant <- args$input
tab.design <- args$design

if (!is.null(tab.design)) {
    stopifnot(file.exists(tab.design))
}



# load Quant files --------------------------------------------------------

Load_Quant_Directory <- function(dir.quant) {
  lst.sample.quant.files <- Serial_Sample_Sub_Path(SUB_PATH_QUANT)
  
  lst.quant.df <- list()
  
  for (sample_name in names(lst.sample.quant.files)) {
      lst.quant.df[[sample_name]] <- LoadQuant(lst.sample.quant.files[[sample_name]], sample_name = sample_name)
  }
  
  res.df <- dplyr::bind_rows(lst.quant.df)
  return(res.df)
}


if (file_test("-d", path.quant)) {
    #' if this is a folder , we take it as the quantification result folder
    df.count <- Load_Quant_Directory(path.quant)
    
} else {
    #' or given it a tabulate file
    df.count <- data.table::fread(path.quant)
}


# load desgin tab file  ---------------------------------------------------

exp.design <- data.table::fread(tab.design) %>% data.frame(stringsAsFactors = F)

#' check data integrity

# start DE ----------------------------------------------------------------
df.DE <- DESeq2::DESeqDataSetFromMatrix(countData = as.matrix(df.count[, -1]),
                                        colData = exp.design %>% dplyr::select(-sample, -formula),
                                        design = as.formula(as.character(exp.design$formula[1])))


obj.de <- DESeq2::DESeq(df.DE)

res.de <- DESeq2::results(obj.de)

# export the result  ------------------------------------------------------

