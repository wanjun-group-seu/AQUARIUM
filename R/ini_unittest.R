# test --------------------------------------------------------------------
base::rm(list = base::ls())

library(magrittr)

# "../test/minimal.cfg" %>% Load_INI_RAW() %>% lst2table() %>% make_adjace_list_from_table() %>% DetectCircleInAdjacentList()
#
#
# # path.ini <- "./test/minimal.functional.cfg"
# path.ini <- "../test/minimal.cfg"
# res.raw <- Load_INI_RAW(path.ini)
# tab.raw <- lst2table(res.raw)
#
# # sapply(tab.raw$value[1], helper_extract_infer)
# al.ini <- make_adjace_list_from_table(tab.raw)



# source("ini.R")
# al.ini2 <-
#     "../test/minimal_missing_value.cfg" %>% LoadINIRaw() %>% TableFromRawList() %>% AdjaceListFromTable()
# stopifnot(!DetectCircleInAdjacentList(al.ini2))
# 
# 
# al.ini <-
#     "../test/minimal_cycle.cfg" %>% LoadINIRaw() %>% TableFromRawList() %>% AdjaceListFromTable() #%>% DetectCircleInAdjacentList()
# stopifnot(DetectCircleInAdjacentList(al.ini))
# 
# 
# 
# 
# base::rm(list = base::ls())
# source('ini.R')
# obj.ini <- "../test/minimal.functional.cfg" %>% LoadINIRaw()
# tab.ini <- obj.ini %>% TableFromRawList()
# al.ini <- tab.ini %>% AdjaceListFromTable()
# 
# uid <- "a:3"
# val.raw <- ValueFromTable(tab.ini)
# adjact.list <- al.ini
# 
# # ExtractAbsValue <- function(uid, abs.val.raw, adjact.list){}
# 
# 
# 
# abs.val <- list()
# 
# stopifnot(AbsValue("a:1", val.raw, adjact.list, abs.val) == "a1")
# 
# stopifnot(AbsValue("a:3", val.raw, adjact.list, abs.val) == "b2")
# 
# stopifnot(AbsValue("b:1", val.raw, adjact.list, abs.val) == "a2 + ")
# 
# 
# stopifnot("b2" == AbsValue("c:y", val.raw, adjact.list, abs.val))
# 
# 
# obj.ini <- "../test/minimal_missing_value.cfg" %>% LoadINIRaw()
# tab.ini <-  obj.ini %>% TableFromRawList() 
# al.ini <- tab.ini %>% AdjaceListFromTable()
# 
# stopifnot(BrokenReferences(al.ini, tab.ini) == "c:1")
# 
# #WriteDownLst(obj.ini, "../test/rewrite.cfg")
# 
# AssginAbsValueWhole(LoadINIRaw("../test/minimal.functional.cfg"))
# 
# read.ini('../test/minimal.functional.cfg')




# test phase 2 ------------------------------------------------------------
source("./ini.R")


stopifnot(!is.null(ReadINI("../test/minimal_cycle.cfg")))

stopifnot(HasCycleRef("../test/minimal_cycle.cfg"))

stopifnot(HasNullRef("../test/minimal_missing_value.cfg"))


WriteINI(ReadINI("../test/minimal.functional.cfg"), "../test/copy.minimal.functional.cfg")

WriteINI(ReadINI("../test/minimal.functional.cfg"), "../test/abs.copy.minimal.functional.cfg", use.abs = T)


WriteCSV(ReadINI("../test/minimal.functional.cfg"), "../test/abs.copy.minimal.functional.tab", use.abs = T)


