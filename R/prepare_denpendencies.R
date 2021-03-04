#! /usr/bin/env Rscript

BasicInstall <- function(p) {
    install.packages(p, repos = "http://cran.r-project.org")
}

BioconductorInstall <- function(p) {
    source("https://bioconductor.org/biocLite.R")
    biocLite(p)
}

CheckDependencies <- function(package_list, func_install=BasicInstall) {
    for (p in package_list) {
        if (!suppressWarnings(suppressMessages(require(
            p,
            character.only = TRUE,
            quietly = TRUE,
            warn.conflicts = FALSE
        )))) {
            #' try use cran as package source.
            func_install(p)
            suppressWarnings(suppressMessages(library(
                p,
                character.only = TRUE,
                quietly = TRUE,
                warn.conflicts = FALSE
            )))
        }
    }
}


CheckDependencies(c("tidyverse", "optparse", "data.table"))
CheckDependencies(c("DESeq2", "GenomicFeatures"), func_install = BioconductorInstall)

