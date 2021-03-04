#! /usr/bin/env Rscript

# description of this script ----------------------------------------------
#' FUNCTIONALITY
#' https://github.com/dvdscripter/ini/blob/master/R/ini.R

#' USAGE

#' ARGUMENTS


.trim <- function(x) {
    sub('^\\s*(.*?)\\s*$', '\\1', x)
}


#' Title a inner helper function for detecting status of each ini line
#'
#' @param line.str a raw string from ini file
#'
#' @return one of 3 statust : "ignore" "section" "key"
#' @export
#'
#' @examples
.GuessLineStatus <- function(line.str) {
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


#' Title reads a ini file and return a list
#'
#' @param path.ini path to a ini file
#'
#' @return
#' @export
#'
#' @examples
.LoadINIRaw <- function(path.ini) {
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
        status.line <- .GuessLineStatus(line.this)
        
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


.TableFromRawList <- function(lst.ini) {
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
    
} # end of Lst2Table


.AdjaceListFromTable <- function(tab.ini) {
    stopifnot(all(c("uid", "value", "section", "key") %in% names(tab.ini)))
    
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
    for (index.row in seq_along(tab.ini$uid)) {
        val.this <- tab.ini$value[index.row]
        if (helper_detect_infer(val.this)) {
            uid.this <- tab.ini$uid[index.row]
            res[[uid.this]] <- helper_extract_infer(val.this)
        }
        
    }
    return(res)
}

#' detect whether there is a cycle in this ini file
.DetectCircleInAdjacentList <- function(adjace.list.ini) {
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


.BrokenReferences <- function(adjace.list.ini, tab.ini) {
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


.ValueFromTable <- function(tab.ini) {
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


.AbsValue <- function(uid,
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
                .AbsValue(neigh,
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


.RebuildAbsList <- function(tab.ini,  adj.ini) {
    raw.value.lst <- .ValueFromTable(tab.ini)
    pool.abs.val <- list()
    res <- list()
    for (section.key in tab.ini$uid) {
        section <- gsub("\\:.+$", "", section.key)
        key <- gsub("^.+\\:", "", section.key)
        #print(paste0(section, " --- ", key))
        res[[section]][[key]] <-
            .AbsValue(section.key, raw.value.lst, adj.ini, pool.abs.val)
        
    }
    return(res)
}

.RebuildList <- function(tab.ini) {
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


.WriteDownLst <- function(obj.ini, file.path) {
    obj.file <- file(file.path, open = "w")
    
    for (sect in names(obj.ini)) {
        writeLines(paste0("[", sect, "]"), obj.file)
        section <- obj.ini[[sect]]
        for (key in names(section)) {
            writeLines(paste0(key, " = " , section[[key]]), obj.file)
        }
        writeLines("", obj.file)
    }
    
    close(obj.file)
}


# interface to outside --------------------------------------------------

.ExpandRaw <- function(raw) {
    res <- list()
    res$raw <- raw
    res$tab <- .TableFromRawList(raw)
    res$adj.lst <- .AdjaceListFromTable(res$tab)
    class(res) <- c("iniparser", "list")
    return(res)
}

ReadINI <- function(path.ini) {
    raw <- .LoadINIRaw(path.ini)
    res <- .ExpandRaw(raw)
    return(res)
}

.makesure <- function(ini) {
    if (is.character(ini)) {
        obj.ini <- ReadINI(ini)
    } else {
        stopifnot(any(class(ini) == "iniparser"))
        obj.ini <- ini
    }
    return(obj.ini)
}

HasCycleRef <- function(ini) {
    obj.ini <- .makesure(ini)
    .DetectCircleInAdjacentList(obj.ini$adj.lst)
}

HasNullRef <- function(ini, verbose = F) {
    obj.ini <- .makesure(ini)
    ref.null <- .BrokenReferences(obj.ini$adj.lst, obj.ini$tab)
    if (verbose) {
        return(ref.null)
    } else {
        return(length(ref.null) > 0)
    }
}

Export <- function(ini.object,
                    use.abs = F,
                    tabulate = F) {
    this.obj <- .ExpandRaw(ini.object$raw)
    
    if (HasCycleRef(this.obj)) {
        stop("has cycle ref in your ini object")
    }
    
    if (HasNullRef(this.obj)) {
        stop("has NULL ref in your ini object")
    }
    
    if (use.abs) {
        lst.final <- .RebuildAbsList(ini.object$tab, ini.object$adj.lst)
    } else {
        lst.final <- ini.object$raw
    }
    
    if (tabulate) {
        return(.TableFromRawList(lst.final))
        
    } else {
        return(lst.final)
    }
}

WriteINI <- function(obj.ini, path.to.export, use.abs = F) {
    .WriteDownLst(Export(obj.ini, use.abs, tabulate = F),
                  path.to.export)
}

WriteCSV <- function(obj.ini, path.to.export, use.abs = F) {
    write.table(
        Export(obj.ini, use.abs, tabulate = T),
        path.to.export,
        quote = F,
        sep = "," ,
        col.names = T,
        row.names = F
    )
}
