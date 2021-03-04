#! /usr/bin/env Rscript
# load packages -----------------------------------------------------------

library(tidyverse)
library(Biostrings)

helper_CountTableFaList <- function(biostring.fa.obj) {
    Biostrings::alphabetFrequency(biostring.fa.obj) %>% as_tibble() %>% mutate("name" = names(biostring.fa.obj)) %>% select(one_of(c("name", "A", "C", "G", "T"))) %>% mutate("len" = A + C + G + T)
}


helper_ExtractInfoFromFaHeader <- function(biostrings.fa.obj) {
    fa.name.trimed <-
        names(biostrings.fa.obj) %>% stringr::str_replace_all(c("gene\\=" = "","CDS\\=.*$" = "", " $" = ""))
    
    tibble::tibble("raw" = fa.name.trimed) %>% tidyr::separate(
        col = "raw",
        into = c("transcript", "gene"),
        sep = "\\ "
    ) %>% mutate('name' = names(biostrings.fa.obj))
}

GetInfoFa <- function(fa.path, sample_name) {
    lst.fa.obj <-
        Biostrings::readDNAStringSet(filepath = fa.path, format = "fasta")
    count.fa <- helper_CountTableFaList(lst.fa.obj)
    info.fa <- helper_ExtractInfoFromFaHeader(lst.fa.obj)
    res <- info.fa %>% left_join(count.fa, by = "name")
    if (!is.null(sample_name))
    {
        res <- res %>% mutate(sample_name = sample_name)
    }
    return(res)
}


# end of function ---------------------------------------------------------



ExtractInfoFromSampleSerial <- function(lst.quant.sf) {
    res <- list()
    for (x in names(lst.quant.sf)) {
        tmp <- GetInfoFa(lst.quant.sf[[x]], x)
        res[[x]] <- tmp
    }
    
    df.whole <- do.call(dplyr::bind_rows, res)
    
    df.total <-
        df.whole %>% select(-sample_name, -gene, -name) %>% dplyr::distinct()
    return(df.total)
}


LoadSampleSerial <- function(dir.quant) {
    safe.dir <- sub("\\/$", "", dir.quant)
    samples.raw <-
        list.dirs(dir.quant, full.names = F, recursive = F)
    path.sf <-
        paste(safe.dir, "/", samples.raw, "/", "final.fa", sep = "")
    res <- list()
    for (x in seq_along(samples.raw)) {
        res[[samples.raw[x]]] <- path.sf[x]
    }
    
    return(res)
}


# invoke in shell  --------------------------------------------------------

args <- commandArgs(trailingOnly = T)
dir.in <- args[1]
file.out <- args[2]

dir.in  %>% LoadSampleSerial() %>% ExtractInfoFromSampleSerial() %>% write.table(file = file.out,
                                                                                quote = F,
                                                                                row.names = F)
