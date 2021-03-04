#! /usr/bin/env Rscript

# description of this script ----------------------------------------------
#' FUNCTIONALITY

#' USAGE

#' ARGUMENTS


# check dependency  and try install missing package -----------------------


# clean environment and load packages --------------------------------------------------------
base::rm(list = base::ls())





#! /usr/bin/env Rscript

# description of this script ----------------------------------------------
#' FUNCTIONALITY
#' a little tweak from https://github.com/dvdscripter/ini/blob/master/R/ini.R

#' USAGE

#' ARGUMENTS


.trim <- function(x) {
    sub('^\\s*(.*?)\\s*$', '\\1', x)
}


#'  a inner helper function for detecting status of each ini line
#'
#'
#' @param line.str a raw string from ini file
#'
#' @return one of 3 statust : "ignore" "section" "key"
#' @export
#'
#' @examples
.guess_line_status <- function(line.str) {
    sectionREGEXP <- '^\\s*\\[\\s*(.+?)\\s*]'
    # match section and capture section name
    
    keyValueREGEXP <- '^\\s*[^=]+=?.?'
    # match "key = value" pattern
    
    ignoreREGEXP <- '^\\s*[;#]'
    # match lines with ; or # at start
    if (nchar(line.str) < 1 ||
        grepl(ignoreREGEXP, line.str)) {
        return("ignore")
    }
    if (grepl(sectionREGEXP, line.str)) {
        return("section")
    }
    if (grepl(keyValueREGEXP, line.str)) {
        return("key")
    }
    
}


#'  reads a ini file and return a list
#'
#' @param path.ini path to a ini file
#'
#' @return
#' @export
#'
#' @examples
.load_ini_as_list <- function(path.ini) {
    lines.ini <- readLines(path.ini)
    
    try(if (length(lines.ini) < 1) {
        #" exit program
        stop("file is empty...")
    })
    
    
    HelperExtractSection <- function(line.ini.section) {
        sub("\\[(.+)\\]", "\\1", line.ini.section)
    } # end of ExtractSection
    
    
    HelperExtractKey <- function(line.ini) {
        tmp <- sub("=.*", "", line.ini)
        .trim(tmp)
    } # end of ExtractKey
    
    
    HelperExtractValue <- function(line.ini) {
        tmp <- sub(".*=", "", line.ini)
        .trim(tmp)
    } # end of ExtractValue
    
    
    
    section.name <- ""
    res <- list()
    for (li in seq_along(lines.ini))
    {
        line.this <- .trim(lines.ini[li])
        status.line <- .guess_line_status(line.this)
        
        key.this <- NULL
        value.this <- NULL
        if (status.line == "ignore") {
            next
        }
        
        if (status.line == "section") {
            section.name <- HelperExtractSection(line.this)
            res[[section.name]] <- list()
        }
        if (status.line == "key") {
            key.this <- HelperExtractKey(line.this)
            value.this <- HelperExtractValue(line.this)
            res[[section.name]][[key.this]] <- value.this
        }
    }
    return(res)
}

# step1 ,  list -> table
.table_from_raw_list <- function(lst.ini) {
    k_section <- NULL
    k_key <- NULL
    val <- NULL
    for (section in names(lst.ini)) {
        this.section <- lst.ini[[section]]
        for (key in names(this.section)) {
            k_section <- c(k_section, section)
            k_key <- c(k_key, key)
            val <- c(val, this.section[[key]])
        }
    }
    
    uid <- paste(k_section , ":", k_key, sep = "")
    data.frame(
        "uid" = uid,
        "section" = k_section,
        "key" = k_key,
        "value" = val,
        stringsAsFactors = F
    )
    
} # end of table_from_raw_list

#' step2, list -> table -> adjace_list
.adjace_list_from_table <- function(table.form.ini) {
    stopifnot(all(
        c("uid", "value", "section", "key") %in% names(table.form.ini)
    ))
    
    regexpINFER <- "\\$\\{.+?\\:.+?\\}"
    
    helper_detect_infer <- function(chr.val) {
        grepl(regexpINFER, chr.val)
    }
    
    helper_extract_infer <- function(chr.val) {
        raw.lst <- regmatches(chr.val, gregexpr(regexpINFER, chr.val))
        sapply(raw.lst, function(x) {
            sub("\\$\\{(.+?\\:.+?)\\}", "\\1", x)
        })
    }
    
    
    res <- list()
    for (index.row in seq_along(table.form.ini$uid)) {
        val.this <- table.form.ini$value[index.row]
        if (helper_detect_infer(val.this)) {
            uid.this <- table.form.ini$uid[index.row]
            res[[uid.this]] <- helper_extract_infer(val.this)
        }
        
    }
    return(res)
}


#' step4, detect whether there is a cycle in adjace_list
.detect_circle_in_adjacent_list <- function(adjace.list.ini) {
    visit.status.lst <- list()
    for (x in names(adjace.list.ini)) {
        visit.status.lst[[x]] <- 0
    }
    
    graph.al <- adjace.list.ini
    has_circ <- F
    
    DFS <- function(uid, visited.status) {
        visited.status[[uid]] <- -1
        while (length(graph.al[[uid]]) > 0) {
            neighbor.this <- graph.al[[uid]][1]
            graph.al[[uid]] <- graph.al[[uid]][-1]
            
            
            if (is.null(visited.status[[neighbor.this]])) {
                #            cat(uid, neighbor.this)
                next
            }
            if (visited.status[[neighbor.this]] == 0) {
                DFS(neighbor.this, visited.status)
                visited.status[[neighbor.this]] = 1
                
            } else {
                if (visited.status[[neighbor.this]] == -1) {
                    has_circ <<- T
                }
            }
            
        }
    }
    for (x in names(graph.al)) {
        if (visit.status.lst[[x]] == 0) {
            DFS(x, visit.status.lst)
        }
    }
    return(has_circ)
}

#' another step4, detect whether there is a NULL inference
.broken_reference <- function(adjace.list.ini, tab.ini) {
    entry.missing <- NULL
    graph.al <- adjace.list.ini
    
    for (uid in names(graph.al)) {
        if (!uid %in% tab.ini$uid) {
            entry.missing <- c(uid, entry.missing)
        }
        for (neigh in graph.al[[uid]]) {
            if (!neigh %in% tab.ini$uid) {
                entry.missing <- c(neigh, entry.missing)
            }
        }
        
    }
    return(entry.missing)
}


.value_from_table <- function(tab.ini) {
    stopifnot(all(c("uid", "value") %in% names(tab.ini)))
    ZipToList <- function(vec.section.key, vec.value) {
        res <- list()
        for (x in seq_along(vec.section.key)) {
            res[[vec.section.key[x]]] <- vec.value[x]
        }
        return(res)
    }
    ZipToList(tab.ini$uid, tab.ini$value)
}


.abs_value <- function(uid,
                       raw.value.lst.from.table,
                       adjact.list,
                       pool.abs.val = list()) {
    if (uid %in% names(pool.abs.val)) {
        return(pool.abs.val[[uid]])
    }
    
    cur.str <- raw.value.lst.from.table[[uid]]
    if (uid %in% names(adjact.list)) {
        lst.neigh <- adjact.list[[uid]]
        for (neigh in lst.neigh) {
            value.neigh <-
                .abs_value(neigh,
                           raw.value.lst.from.table,
                           adjact.list,
                           pool.abs.val)
            pattern.neigh <- paste0("${", neigh, "}")
            cur.str <-
                gsub(pattern.neigh, value.neigh, cur.str, fixed = T)
        }
    }
    pool.abs.val[[uid]] <- cur.str
    return(pool.abs.val[[uid]])
}


.rebuild_abs_list <- function(tab.ini,  adj.ini) {
    raw.value.lst <- .value_from_table(tab.ini)
    pool.abs.val <- list()
    res <- list()
    for (section.key in tab.ini$uid) {
        section <- gsub("\\:.+$", "", section.key)
        key <- gsub("^.+\\:", "", section.key)
        #print(paste0(section, " --- ", key))
        res[[section]][[key]] <-
            .abs_value(section.key, raw.value.lst, adj.ini, pool.abs.val)
        
    }
    return(res)
}

.rebuild_list <- function(tab.ini) {
    vec.section.key <- tab.ini$uid
    vec.val <- tab.ini$value
    res <- list()
    for (index.vec in seq_along(vec.section.key)) {
        section <- tab.ini$section[index.vec]
        key <- tab.ini$key[index.vec]
        res[[section]][[key]] <- vec.val[index.vec]
    }
    return(list)
}


.write_down_list <- function(obj.ini.lst, file.path) {
    obj.file <- file(file.path, open = "w")
    
    for (sect in names(obj.ini.lst)) {
        writeLines(paste0("[", sect, "]"), obj.file)
        section <- obj.ini.lst[[sect]]
        for (key in names(section)) {
            writeLines(paste0(key, " = " , section[[key]]), obj.file)
        }
        writeLines("", obj.file)
    }
    
    close(obj.file)
}


.expand_raw_ini_list <- function(raw.list.ini) {
    res <- list()
    res$raw <- raw.list.ini
    res$ini <- res$raw
    res$tab <- .table_from_raw_list(raw.list.ini)
    res$adj.lst <- .adjace_list_from_table(res$tab)
    class(res) <- c("iniparser", "list")
    return(res)
}

.self_update <- function(object.ini) {
    if (identical(object.ini$raw, object.ini$ini)) {
        object.ini
    } else {
        print("ok , something has been changed")
        .expand_raw_ini_list(object.ini$ini)
    }
}



.make_sure <- function(ini) {
    if (is.character(ini)) {
        .load_ini_as_list(ini)
    } else {
        ini
    }
}



# interface to outside --------------------------------------------------

#'  Read INI file
#'
#' @param path.ini path to INI file
#'
#' @return list of list
#' @export
#'
#' @examples
read_ini <- function(path.ini) {
    .load_ini_as_list(path.ini)
}


#'  check whether there exists a cyclic reference
#'
#' @param ini string contains the path to INI file , or a iniparser object returned by ReadINI
#'
#' @return bool variable indicates whether a loop exists
#' @export
#'
#' @examples
has_cycle_ref <- function(lst.ini) {
    .detect_circle_in_adjacent_list(.adjace_list_from_table(.table_from_raw_list(lst.ini)))
}


#'  check whether there exists unassigned reference
#'
#' @param ini string contains the path to INI file , or a iniparser object returned by ReadINI
#' @param verbose bool, if True this function will return the list of null refers.
#'
#' @return return True is there exists null reference, if verbose is True, return the list of null reference
#' @export
#'
#' @examples
has_null_ref <- function(ini, verbose = F) {
    lst.ini <- .make_sure(ini)
    tab.lst <- .table_from_raw_list(lst.ini)
    adj.lst <- .adjace_list_from_table(tab.lst)
    ref.null <- .broken_reference(adj.lst, tab.lst)
    if (verbose) {
        return(ref.null)
    } else {
        return(length(ref.null) > 0)
    }
}

#'  Export iniparser object into other format.
#'
#' @param ini.object a 'iniparse' object returned by ReadINI
#' @param use.abs if True, all relative reference will be replaced by absolute value
#' @param tabulate if True, this function will return a data.frame instead a nested list
#'
#' @return a data.frame if tabulate is True, else a list .
#' @export
#'
#' @examples
export_ini <- function(ini.object,
                       use.abs = F,
                       tabulate = F) {
    lst.ini <- .make_sure(ini.object)
    tab.lst <- .table_from_raw_list(lst.ini)
    adj.lst <- .adjace_list_from_table(tab.lst)
    ref.null <- .broken_reference(adj.lst, tab.lst)
    
    if (.detect_circle_in_adjacent_list(adj.lst)) {
        stop("Found cycle reference in your ini")
    }
    
    if (length(ref.null) > 0) {
        stop("Has NULL ref in your ini")
    }
    
    if (use.abs) {
        lst.final <- .rebuild_abs_list(tab.lst, adj.lst)
    } else {
        lst.final <- lst.ini
    }
    
    if (tabulate) {
        return(.table_from_raw_list(lst.final))
        
    } else {
        return(lst.final)
    }
}




#'  write down iniparser object in INI format
#'
#' @param obj.ini an iniparser object returned by ReadINI
#' @param path.to.export path to the result file
#' @param use.abs relative reference will be replaced by absolute value if True
#'
#' @return nothing. a INI format file will be created.
#' @export
#'
#' @examples
write_ini <- function(obj.ini, path.to.export, use.abs = F) {
    .write_down_list(export_ini(obj.ini, use.abs, tabulate = F),
                     path.to.export)
}

#'
#'
#' @param obj.ini an iniparser object returned by ReadINI
#' @param path.to.export path to the result file
#' @param use.abs relative reference will be replaced by absolute value if True
#'
#' @return nothing. a CSV format file will be created.
#' @export
#'
#' @examples
write_csv <- function(obj.ini, path.to.export, use.abs = F) {
    write.table(
        export_ini(obj.ini, use.abs, tabulate = T),
        path.to.export,
        quote = F,
        sep = "," ,
        col.names = T,
        row.names = F
    )
}



# end of ini part ---------------------------------------------------------

FQ_FILENAME_PATTERN_GZ <- "(\\.fq$)|(\\.fastq$)|(\\.fq.gz$)|(\\.fastq.gz$)"

FQ_FILENAME_PATTERN <- "(\\.fq$)|(\\.fastq$)"


guess_sample_info <- function(abs.dir.fq.raw, pattern_fq_filename=FQ_FILENAME_PATTERN) {
    if (is.null(abs.dir.fq.raw)) {
        return(NULL)
    }
    abs.dir.fq <- gsub("\\/$", "", abs.dir.fq.raw)
    lst.fq <-
        list.files(abs.dir.fq, pattern = pattern_fq_filename)
    strip.fq.name <- gsub(pattern_fq_filename, "", lst.fq)
    
    pe_num <-
        regmatches(strip.fq.name, regexec("\\d$", strip.fq.name))
    
    pe_num_simple <-
        as.numeric(unlist(ifelse(nchar(pe_num) == 1 , pe_num, NA)))
    
    info.fq <- data.frame(
        "fq" = strip.fq.name,
        "pe" = pe_num_simple,
        "sample" = gsub("\\_\\d$", "", strip.fq.name),
        "abs.path" = paste0(abs.dir.fq, "/", lst.fq),
        stringsAsFactors = F
    )
    return(info.fq)
} # end of Guess_Sample_Info


extract_seq_info <- function(dir_seq) {
    information.seqs <- guess_sample_info(dir_seq)
    base::split(information.seqs$abs.path, information.seqs$sample)
    
}

dump_cfg_files <- function(lst.of.sample.ini, cfg.dir) {
    for (sample.id in names(lst.of.sample.ini)) {
        path.ini.this.sample <- paste0(cfg.dir, "/", sample.id, ".cfg")
        write_ini(lst.of.sample.ini[[sample.id]], path.ini.this.sample)
    }
}


mix_up_shell_content_this_sample <-
    function(path.cfg,
             path.workflow.detection,
             path.workflow.quant) {
        c(
            "#! /usr/bin/env bash",
            paste("python3", path.workflow.detection, path.cfg, sep = "\t"),
            paste("python3", path.workflow.quant, path.cfg, sep = "\t")
        )
    }

mix_up_shell_content_main <- function(lst.path.shell.script) {
    tmp <- paste("bash", lst.path.shell.script, sep = " ")
    c("#! /usr/bin/env bash", tmp)
}


write_the_shell_script <-
    function(path.shell.should.be,
             lst.str.shell.content) {
        file.shell.op <- file(path.shell.should.be)
        writeLines(lst.str.shell.content, file.shell.op)
        close(file.shell.op)
    }

# start the command args part -----------------------------------------------------

STR_USAGE <-             "

usage:
          Rscript make_serial_config.R template.path fq.dir cfg.dir sh.dir

          template.path: file path of your config file template
          
          fq.dir: folder of fq files
          
          cfg.dir: path to a folder where to store generated config files
          
          sh.dir: optional , a folder to put you PBS shell scripts

Notes:

you need to specify 'pipeline_script_detection' and 'pipeline_script_profile' in META section of template file. 

use abstract path to avoid potential error. 

        
            "


main <- function(template.path, fq.dir, cfg.dir, sh.dir = NULL) {
    #' initialize the ini object and read the correct information from the
    
    if (is.na(template.path)) {
        cat(STR_USAGE)
        quit(save = "no")
    }
    
    if (!dir.exists(fq.dir)) {
        stop("NO FQ DIR!")
    }
    
    if (!dir.exists(cfg.dir)) {
        dir.create(cfg.dir)
    }
    
    if (is.null(sh.dir) || is.na(sh.dir)) {
        sh.dir <- cfg.dir
        cat("sh file per sample would be in the same folder as cfg files\n")
    }
    
    if (!dir.exists(sh.dir)) {
        dir.create(sh.dir)
    }
    
    
    ini.template <- read_ini(template.path)
    
    detection.dir <-  ini.template$GLOBAL$detection_dir
    if (!dir.exists(detection.dir)) {
        cat(paste("\n create folder for detection : ", detection.dir, '\n'))
        dir.create(detection.dir)
    }
    
    quantification.dir <- ini.template$GLOBAL$quant_root_dir
    if (!dir.exists(quantification.dir)) {
        cat(paste("\n create folder for quantification: ", quantification.dir), "\n")
        dir.create(quantification.dir)
    }
    
    #' get script path from the cfg file
    path.work.flow.detection <-
        ini.template$META$pipeline_script_detection
    path.work.flow.profile <-
        ini.template$META$pipeline_script_profile
    
    
    if (!file.exists(path.work.flow.detection)) {
        stop("No valid detection workflow py file in template!")
    }
    
    if (!file.exists(path.work.flow.profile)) {
        stop("No valid profile workflow py file in template!")
    }
    
    lst.fq <- extract_seq_info(fq.dir)
    lst.shell.path <- NULL   #' for the main shell script
    
    for (sample.id in names(lst.fq)) {
        #' for each sample , produce a cfg from the template file
        fq.sorted <- sort(lst.fq[[sample.id]])
        num.fq.this.sample <- length(fq.sorted)
        stopifnot(num.fq.this.sample %in% c(1, 2))
        
        #' make sure the quantification root for this sample exists
        dir.this.sample.quant <-
            paste0(quantification.dir, "/", sample.id)
        if (!dir.exists(dir.this.sample.quant)) {
            dir.create(dir.this.sample.quant)
        }
        
        #' copy the ini here
        ini.this <- ini.template
        ini.this$GLOBAL$sample_id <- sample.id
        ini.this$CUSTOM$quant_dir <- dir.this.sample.quant
        ini.this$CIRI[["--in"]] <-
            paste0(detection.dir, "/", sample.id, ".sam")
        
        if (num.fq.this.sample == 2) {
            #' PE
            
            ini.this$CUSTOM$r1 = fq.sorted[1]
            ini.this$CUSTOM$r2 = fq.sorted[2]
            ini.this$CIRI[["--seqs"]] = "${CUSTOM:r1} ${CUSTOM:r2}"
            ini.this$CIRC_PROFILE[["-1"]] = "${CUSTOM:r1}"
            ini.this$CIRC_PROFILE[["-2"]] = "${CUSTOM:r2}"
            ini.this$CIRC_FULL[["-1"]] = "${CUSTOM:r1}"
            ini.this$CIRC_FULL[["-2"]] = "${CUSTOM:r2}"
            
        } else {
            #' SE
            ini.this$CUSTOM$r1 = fq.sorted[1]
            ini.this$CIRI[["--seqs"]] = "${CUSTOM:r1}"
            ini.this$CIRC_PROFILE[["-r"]] = "${CUSTOM:r1}"
            ini.this$CIRC_FULL[["-r"]] = "${CUSTOM:r1}"
        }
        
        
        path.ini.this.sample <-
            paste0(cfg.dir, "/", sample.id, ".cfg")
        write_ini(ini.this, path.ini.this.sample) #' this is where we produce the ini files
        
        #" make a shell script for each sample
        shell.content.this.sample <-
            mix_up_shell_content_this_sample(path.ini.this.sample,
                                             path.work.flow.detection,
                                             path.work.flow.profile)
        path.shell.script.this.sample <-
            paste0(sh.dir, "/", sample.id, ".sh")
        
        lst.shell.path <-
            c(lst.shell.path, path.shell.script.this.sample)
        write_the_shell_script(path.shell.script.this.sample,
                               shell.content.this.sample)
        
        
    } # end of loop on each sample
    
    
    #' write the main script for this serial samples.
    path.to.main.shell <- paste0(dirname(sh.dir), "/", "batch.sh")
    
    cat(paste("\n", "path to batch shell: ", path.to.main.shell, "\n"))
    write_the_shell_script(path.to.main.shell,
                           mix_up_shell_content_main(lst.shell.path))
    
    
}


# the real runing part -----------------------------------------------------------------

cli_args = commandArgs(trailingOnly = T)


template.path <- cli_args[1]
fq.dir <- cli_args[2]
cfg.dir <- cli_args[3]
sh.dir <- cli_args[4]


# be careful of user given arguments
main(template.path, fq.dir, cfg.dir, sh.dir)
