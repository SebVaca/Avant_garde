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
setwd(file.path(args[4]))
setwd("AvantGardeDIA/")
getwd()


#Folders 
Folder_1="DataFormatting"
dir.create(file.path(getwd(),Folder_1),showWarnings = F)
Output.file=file.path(getwd(),Folder_1)
dir.create(file.path(Output.file,"MultiFile"),showWarnings = F)
MultiFile.path=file.path(Output.file,"MultiFile")
Folder_2="TempFiles"
dir.create(file.path(getwd(),Folder_1,Folder_2),showWarnings = F)
dir.output=file.path(getwd(),Folder_1,Folder_2)
Multifile_path=file.path(Output.file,"MultiFile")


ParamsFile="AvG_Params.R"
source(ParamsFile)
RefinementWorkflow="OnlyScoring"
write.table(paste0("RefinementWorkflow='",RefinementWorkflow,"'"),
            file = file.path(getwd(),"AvG_Params.R"),
            append = T,quote = F,row.names = F,col.names = F)


library(AvantGardeDIA,quietly = T)
cat(paste(paste0("AvantGardeDIA_",RefinementWorkflow," is running."),"The analisis is done using N-1 cores where N is the number of cores detected in this device. This might take a few minutes.",sep="\n"))


LaunchParallelTasks_ReadFromDB(ParamsFile,RefinementWorkflow)
Read_AndFormatResults()
print("AvantGardeDIA_NO: Done! Have a nice day!")

rm(list=ls())
