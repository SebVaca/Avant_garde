library("devtools")
library(roxygen2)
setwd("C:/Users/svaca/Documents/Code_saved/GeneticAlgorithm/20190118_AvantGardeForPublication/r_package/")
create("AvantGardeDIA")
setwd("./AvantGardeDIA")
# copy r file
document()
devtools::use_package("dplyr")
devtools::use_package("tidyr")
devtools::use_package("stringr")
devtools::use_package("data.table")
devtools::use_package("sqldf")
devtools::use_package("foreach")
devtools::use_package("doSNOW")
devtools::use_package("GA")

setwd("..")
install("AvantGardeDIA")

# A<-sessionInfo()
# B<-rbindlist(lapply(A$otherPkgs,function(L){data.frame(Pckg=L$Package,version=L$Version)})) %>% mutate(version2=paste0("(>= ",version,")")) %>% mutate(version3=paste0(Pckg," ",version2))
# B
# paste(B$version3,collapse=",")


# setwd("C:/Users/svaca/Documents/Code_saved/GeneticAlgorithm/20180228_AvantGardeDIA_v1.15/TestPackage/")
# library("devtools")
# library(roxygen2)
# install("AvantGardeDIATest")

# library(AvantGardeDIATest8)
# #LatestVersion="C:/Users/svaca/Documents/Code_saved/GeneticAlgorithm/20180228_AvantGardeDIA_v1.15/20180312_Avant-gardeDIA_v18_Unified_Annotated_foreach.R"
# ParamsFile="C:/Users/svaca/Documents/Code_saved/GeneticAlgorithm/20180228_AvantGardeDIA_v1.15/20180308_ParamsFile.R"
# source(ParamsFile)
# #source(LatestVersion)
# AvantGardeDIA_DB(D.file.name = D.file.name,RefinementWorkflow = "GlobalRefinement",Quantitative.Report.Input,ParamsFile)


#### Create TAR.GZ file
library("devtools")
library(roxygen2)
setwd("C:/Users/svaca/Documents/Code_saved/GeneticAlgorithm/20190118_AvantGardeForPublication/r_package/AvantGardeDIA/")
build()


#### Remove package
remove.packages("AvantGardeDIATest7", lib="~/R/win-library/3.4")