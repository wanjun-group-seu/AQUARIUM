#! /usr/bin/env Rscript

rm(list = ls())
#source("./sum_up.R")

library(dplyr)
library(data.table)

#' Title load a salfish quant.sf file
#'
#' @param sf
#' @param ...
#'
#' @return : a dataframe
#' @export
#'
#' @examples
LoadQuant <- function(sf, ...) {
    dd.sf <- data.table::fread(
        sf,
        header = T,
        sep = "\t",
        stringsAsFactors = F,
        ...
    )
    return(
        data.frame(
            "Name" = dd.sf$Name,
            "Length" = dd.sf$EffectiveLength,
            "TPM" = dd.sf$TPM,
            "NumReads" = dd.sf$NumReads,
            stringsAsFactors = F
        )
    )
}


#' Title load a salfish quantification result with a sample name
#'
#' @param path.sf : path to the .sf file
#' @param sample.name : a name for this sample
#'
#' @return : a data.frame with a sample column
#' @export
#'
#' @examples
GetSampleExpress <- function(path.sf, sample.name) {
    cbind(LoadQuant(path.sf),
          # here, chose your origin of quant format.
          "sample" = sample.name,
          stringsAsFactors = F)
}


#' Title remove the ".r" surfix made by gffread
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
TrimStrand <- function(x) {
    if (base::endsWith(x, ".r")) {
        x <- substr(x, start = 1, stop = nchar(x) - 2)
    }
    return(x)
}

RowNameAsFirstCol <- function(iso.level.exp, name.this.column) {
    entry <- row.names(iso.level.exp)
    name.col <- colnames(iso.level.exp)
    iso.level.exp  <-
        as.data.frame(iso.level.exp,
                      col.names = name.col,
                      stringAsFactors = F)
    iso.level.exp <- cbind(entry, iso.level.exp)
    names(iso.level.exp) <-
        c(name.this.column , names(iso.level.exp)[-1])
    iso.level.exp[, 1] <-
        as.character(levels(iso.level.exp[, 1]))[iso.level.exp[, 1]]
    return(iso.level.exp)
}



SumItUp <-
    function(df.whole.data,
             origin.type,
             to.be.sumed,
             by.what) {
        tmp <- df.whole.data %>% filter(origin == origin.type)
        tapply(tmp[[to.be.sumed]], list(tmp[[by.what]], tmp$sample), function(x) {
            sum(x, na.rm = T)
        })
    }

NA20 <- function(x) {
    x[is.na(x)] <- 0
    return(x)
}

MakeUpTable <- function(sum.up.data,
                        whole.entry.list,
                        entry.name) {
    
    x.out <-
        sum.up.data[match(whole.entry.list, sum.up.data[, 1]), -1]
    x.out[is.na(x.out)] <- 0
    if (is.null(dim(x.out)))
    {
        
        return(data.frame("name" = whole.entry.list, "val" = x.out, stringsAsFactors = F))
    } else{
        row.names(x.out) <- whole.entry.list
        x.out <- RowNameAsFirstCol(x.out, entry.name)
        return(x.out)
    }
}



str.mrna <- "mrna"
str.linc <- "linc"
str.circ <- "circular"
str.gene <- "gene"
str.iso <- "transcript"
str.count <- "NumReads"
str.tpm <- "TPM"

#' this should be IMPORTANT!!!

#' @param df.total : a data.frame contains the whole information of a series of expriments
#' @param circ.linear.linc : a string shows which kind of transcript you want , you can use str.mrna/str.linc/str.circ
#' @param tpm.count : what kind of measurement , you can use str.tpm/str.count
#' @param gene.transcript : what unit you want , you can use str.iso/str.gene

adhocSumUp <-
    function(df.total,
             circ.linear.linc,
             tpm.count,
             gene.transcript) {
        
        df.total %>% SumItUp(circ.linear.linc, tpm.count, gene.transcript) %>% RowNameAsFirstCol(gene.transcript) %>% NA20() %>% MakeUpTable(sort(unique(df.total[[gene.transcript]])), gene.transcript)
    } # end of adhocSumUp


#

#



# end of function ====

adhocLoadListMappingFile <- function(lst.path.mapping.info, ...) {
    path.all.files <- c(lst.path.mapping.info, ...)
    info.mapping.list <- lapply(path.all.files, function(path.x) {
        data.table::fread(path.x,
                          header = T,
                          stringsAsFactors = F)
    })
    
    dplyr::bind_rows(info.mapping.list) %>% dplyr::distinct()
} # end of adhocLoadListMappingFile


adhocDecoratedMappingInfo <-
    function(mapping.of.transcript.and.gene,
             mapping.of.gene.id.and.name = NULL) {
        df.tmp <-
            mapping.of.transcript.and.gene %>% mutate(transcript = sapply(iso, TrimStrand))
        
        if (!is.null(mapping.of.gene.id.and.name)) {
            df.tmp <-
                df.tmp %>% mutate(gene_name = mapping.of.gene.id.and.name$name[match(gene, mapping.of.gene.id.and.name$id)])
        }
        
        return(df.tmp)
    } # end of adhocDecoratedMappingInfo


adhocDecoratedExpressionWithGeneInfo <-
    function(df.expression.only, df.mapping.info) {
        cbind(df.expression.only, df.mapping.info[match(df.expression.only$Name, df.mapping.info$iso),]) %>% select(-iso)
    } # end of adhocDecoratedExpressionWithGeneInfo


# example ====







############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

#all.your.args <- c("./summarization/host.info.tab" , "./summarization"," p1h0", "/home/wanjun/timecourse/p1_0h/quant_ciri/profile_result/quant.sf", "p1h12","/home/wanjun/timecourse/p1_12h/quant_ciri/profile_result/quant.sf", "p1h24", "/home/wanjun/timecourse/p1_24h/quant_ciri/profile_result/quant.sf" )

#all.your.args <- c("./test/host.info", "./test/sum/", "p1h0", "./test/p1h0/profile_result/quant.sf")


all.your.args <- commandArgs(trailingOnly = T)

if (length(all.your.args) <= 2) {
    stop("=========================================================
functions to sum up your quantification result by type and source, 

usage:

Rscript this.R host_info output_path sample_name_1 sample_quant_path_1 [sample_name_2 sample_quant_path_2] ... ..."
         )    
}
stopifnot(length(all.your.args) > 2)


path.to.mapping.info <- all.your.args[1]

report.path <- all.your.args[2]
report.path <- gsub("\\/$", "", report.path, perl = T)

args.left <- all.your.args[-c(1, 2)]

stopifnot(length(args.left) %% 2 == 0)

num.pairs <- (length(args.left)) / 2

ind.names <- 2 * (1:num.pairs) - 1
ind.path <- ind.names + 1

sample.info <- data.frame("name" = args.left[ind.names],
                          "path" = args.left[ind.path],
                          stringsAsFactors = F)





#' 读取基因组的结构信息, 转录本的归属等.
total.mapping.info <-
    adhocDecoratedMappingInfo(adhocLoadListMappingFile(path.to.mapping.info))


#' 读取多个定量结果文件.
# total.expression <- do.call(rbind,
#                          list(GetSampleExpress("./lnc_level/li.sf", "li"),
#                               GetSampleExpress("./lnc_level/lu.sf", "lu"),
#                               GetSampleExpress("./lnc_level/sun.sf", "sun"),
#                               GetSampleExpress("./lnc_level/zhao.sf", "zhao"))
# )
total.expression <-
    do.call(rbind, lapply(sample.info$name, function(x, df) {
        path.x <- df$path[match(x, df$name)]
        GetSampleExpress(path.sf = path.x, sample.name = x)
    }, df = sample.info))


#' 给定量结果添加基因组信息.
total.expression <-
    adhocDecoratedExpressionWithGeneInfo(total.expression, total.mapping.info)

# 汇总数据 ====

#' 调用 adhocSumUp 来汇总数据.  三个参数 的含义可见函数定义.
# tpm.circ.gene <-
#     total.expression %>% adhocSumUp(str.circ, str.tpm, str.gene)
# tpm.circ.iso <-
#     total.expression %>% adhocSumUp(str.circ, str.tpm, str.iso)
# 
# tpm.mrna.gene <-
#     total.expression %>% adhocSumUp(str.mrna, str.tpm, str.gene)
# tpm.mrna.iso <-
#     total.expression %>% adhocSumUp(str.mrna, str.tpm, str.iso)
# 
# count.linc.iso <-
#     total.expression %>% adhocSumUp(str.linc, str.count, str.iso)
# count.linc.iso <-
#     total.expression %>% adhocSumUp(str.linc, str.count, str.iso)

export.pattern <- expand.grid(c(str.circ, str.linc, str.mrna), c(str.tpm, str.count), c(str.iso, str.gene))

on_each_pattern <- function(source.entry.and.unit.exp.and.unit.entry, df.exp.data){
    source.entry = source.entry.and.unit.exp.and.unit.entry[1]
    unit.exp = source.entry.and.unit.exp.and.unit.entry[2]
    unit.entry = source.entry.and.unit.exp.and.unit.entry[3]
    path.full <- paste(report.path, "/", paste(source.entry, "_", unit.exp, "_", unit.entry, ".tab", sep = ""), sep = "")
    
    
    
    print(path.full)
    
    df.exp.data %>% adhocSumUp(source.entry, unit.exp, unit.entry) %>% data.table::fwrite(file = path.full, sep = "\t", row.names = F)
    
}


apply(export.pattern, 1, on_each_pattern, df.exp.data = total.expression)

print("all is done! ")

#' 输出分类的汇总结果.
#data.table::fwrite(tpm.circ.gene, file = "./lnc_level/report/gene.circ.csv")
