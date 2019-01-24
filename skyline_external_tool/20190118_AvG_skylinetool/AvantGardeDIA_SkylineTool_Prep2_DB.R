rm(list=ls())
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(sqldf)
library(foreach)
library(doSNOW)
library(GA, quietly = T)
library(parallel)

args <- commandArgs()
print(args)

## Install AvantGarde Package if it is not already installed
a<-installed.packages()
packages<-a[,1]
if (!is.element("AvantGardeDIA",packages) || packageVersion("AvantGardeDIA") < "0.0.0.0" ){
  install.packages(Avg_RPackage_path,repos = NULL)
}


setwd(file.path(args[5]))
getwd()
setwd("AvantGardeDIA/")

ParamsFile="AvG_Params.R"
source(ParamsFile)

D.file.name=args[4]


## Folders
Folder_1="DataFormatting"
dir.create(file.path(getwd(),Folder_1),showWarnings = F)
Output.file=file.path(getwd(),Folder_1)
dir.create(file.path(Output.file,"MultiFile"),showWarnings = F)
MultiFile.path=file.path(Output.file,"MultiFile")
Folder_2="TempFiles"
dir.create(file.path(getwd(),Folder_1,Folder_2),showWarnings = F)
dir.output=file.path(getwd(),Folder_1,Folder_2)
Multifile_path=file.path(Output.file,"MultiFile")

library(AvantGardeDIA)
AvantGardeDIA::AvantGardeDIA_InChunks_DB(D.file.name = D.file.name,RefinementWorkflow = "GlobalRefinement",ClusterWorkflow = TRUE)


rm(list=ls())
