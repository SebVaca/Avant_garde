### Formatting and Indexing

#' AvantGardeDIA_CreateParamsFile
#'
#' This function is AvantGardeDIA_CreateParamsFile
#' @param 
#' @keywords AvantGardeDIA
#' @export
#' @examples
AvantGardeDIA_CreateParamsFile<-function(){
  dir.create(file.path(getwd(),"AvantGardeDIA"),showWarnings = F)
  
  # Create AvantGardeDIA Parameters file
  A<-"Name_Tag='Name_Tag' #Short Name that will be added to all output files. (Must be in quotes)
  
  # I- DataFormating Parameters
  MinimalInitialNumberOfTransitions_Filter=5 # Peptides with less than N transitions will be removed and not included in the SQLite Database
  RemoveTransitionBelowOrdinal_Filter=3  # Transition ordinal below or equal to this value will be removed (for example if Max_TransitionIonOrdinal=3 transitions b1,b2,b3,y1,y2,y3 are removed) 
  RemoveSharedTransitionsBetweenLightAndHeavy_Filter=TRUE # Should be TRUE for DIA,FALSE for PRM.
  ReadNumberOfLines=20000 # The CSV file is read in chunks. This parameter sets the number of lines in each chunk.
  
  # II- Optimization Parameters
  NonZeroBaselineChromatogram=FALSE # TRUE for data having high noise baseline level (e.g. SWATH data on Q-TOF), FALSE (noise-reduced data, e.g demultiplexed Overlap DIA Data acquired on Q-Exactive series)
  MinimalNumberOfTransitionsAfterOptimization=4 # The optimization will find the best solution that will have at least this number of transitions.
  KeepPeptidesWithLowerNumberOfTransitions=FALSE # Keep (or not) a peptide if it has a lower number of transitions than the limit above (recomended FALSE for DIA). The transition refinement refinement will not be done on these peptides but the peak boundary refienement and scoring will be performed.
  TopN_RankedbyIntensity=6 # Number of most intense transitions in teh spectral library (N) among which at least n transitions need to be present.
  MinimalNumberOfTransitionsAmongTopN=2 # Number of transitions (n) that need to be present among the N most intense transitions in the spectral library.
  alpha=0.005
  Beta=0.05
  SpectralLibraryDotProduct_limit=0.7 # 
  MassError_Tolerance=10 # Mass error tolerance (in PPM)
  MassError_CutOff=20 # Mass error cut-off (in PPM)
  MinimalNumberOfConsecutivePoints=3 # Minimal number of consecutive points above the limit of noise to consider a signal as a potential peptide.
  MinimalIntensityPercentagePerTransition=3 # Minimal percentage of the maximum intensity of each transtion below which any point is no longer considered.
  UseHeavyPeakBoundariesForLight=FALSE # Use the heavy peptide peak boundaries for the light peptides.
  
  ## Folders
  Folder_1='DataFormatting'
  Output.file=file.path(getwd(),Folder_1)
  MultiFile.path=file.path(Output.file,'MultiFile')
  Folder_2='TempFiles'
  dir.output=file.path(getwd(),Folder_1,Folder_2)
  Multifile_path=file.path(Output.file,'MultiFile')
  RefinementWorkflow='GlobalRefinement'"
  
  write.table(A,file = file.path(getwd(),"AvantGardeDIA","AvG_Params.R"),row.names = F,col.names = F,quote = F)
  print("Create AvantGardeDIA parameters file: Done!")
}

#' Indexing
#'
#' This function is Indexing
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Indexing()
Indexing<-function(D){
  ### Create indices
  D$ID_FragmentIon_charge<- D %>% group_indices(FragmentIon,ProductCharge)
  #D$ID_Analyte<- D %>% group_indices(PeptideModifiedSequence,PrecursorCharge,IsDecoy)
  D$ID_Analyte<- D %>% group_indices(ProteinName,PeptideModifiedSequence,PrecursorCharge,IsDecoy)
  D$ID_Analyte <- with(rle(as.character(D$ID_Analyte)), rep(seq_along(values), lengths))
  D$ID_Rep<-D %>% group_indices(FileName)
  
  
  
  ## MetaData
  M<-D %>% select(-which(str_detect(names(D),pattern = "Interpolated"))) %>%
    select(FileName,ProteinName,PeptideModifiedSequence,ModifiedSequence,PrecursorMz,
           PrecursorCharge,ProductMz,ProductCharge,FragmentIon,IsotopeLabelType,TransitionLocator,Quantitative,IsDecoy,
           ID_FragmentIon_charge,ID_Analyte,ID_Rep,PrecursorResultLocator) %>% distinct()
  M_Analyte<-M %>%select(ID_Analyte,ProteinName,PeptideModifiedSequence,PrecursorCharge,IsDecoy) %>% distinct()
  M_Replicate<-M %>% select(ID_Rep,FileName) %>% distinct()
  M_TransitionsType<-M %>% select(ID_FragmentIon_charge,ProductCharge,FragmentIon) %>% distinct() %>% arrange(ID_FragmentIon_charge)
  M_Transitions<-M %>% select(ID_Analyte,ID_FragmentIon_charge,ProteinName,PeptideModifiedSequence,ModifiedSequence,PrecursorMz,IsotopeLabelType,IsDecoy,
                              PrecursorCharge,ProductMz,ProductCharge,FragmentIon,TransitionLocator,Quantitative) %>% distinct()
  M_PrecursorResult<-M %>%select(PrecursorResultLocator,ProteinName,PeptideModifiedSequence,IsotopeLabelType,PrecursorCharge,IsDecoy,FileName) %>% distinct()
  
  return(MetaData=list(M_Analyte=M_Analyte,
                       M_Replicate=M_Replicate,
                       M_TransitionsType=M_TransitionsType,
                       M_Transitions=M_Transitions,
                       M_PrecursorResult=M_PrecursorResult))
}

#' Filter_Na_Shared_Or_LowMassTransitions
#'
#' This function is Filter_Na_Shared_Or_LowMassTransitions
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Filter_Na_Shared_Or_LowMassTransitions()
Filter_Na_Shared_Or_LowMassTransitions<-function(D1){
  # Remove NA replictaes (Non integrated peaks) and remove transitions y1-3 and b1-3
  D1<-as.data.frame(D1) %>% filter(!is.na(as.numeric(D1$MinStartTime))) %>%
    mutate(Transition_Filter=paste0(FragmentIon," "))%>%
    mutate(Transition_Filter=substr(Transition_Filter,start=1,stop = str_locate(Transition_Filter,pattern = " ")-1)) %>%
    filter(!grepl(paste(c("precursor",paste0("^y", 0:RemoveTransitionBelowOrdinal_Filter,"$"),paste0("^b",0:RemoveTransitionBelowOrdinal_Filter,"$")), collapse="|"), Transition_Filter)) %>%
    select(-Transition_Filter)
  
  ## Remove transition that are not in light AND in the heavy version
  A0<-D1 %>% select(ID_Analyte,ID_FragmentIon_charge,IsotopeLabelType) %>% 
    distinct() %>% group_by(ID_Analyte) %>% summarise(TotalNumLabels=length(unique(IsotopeLabelType)))
  
  A1<-D1 %>% select(ID_Analyte,ID_FragmentIon_charge,IsotopeLabelType) %>% 
    distinct() %>%
    group_by(ID_Analyte,ID_FragmentIon_charge) %>% 
    summarise(NumLabels=length(unique(IsotopeLabelType))) %>%
    left_join(A0, by = "ID_Analyte") %>%
    mutate(diff=TotalNumLabels-NumLabels) %>%
    filter(diff!=0) %>% select(ID_Analyte,ID_FragmentIon_charge) %>%
    ungroup %>%
    mutate(ID_Analyte=as.integer(ID_Analyte),ID_FragmentIon_charge=as.integer(ID_FragmentIon_charge),Remove="DoNotKeep")
  
  if (dim(A1)[1]>=1) {
    D1<-D1 %>% left_join(A1, by = c("ID_Analyte", "ID_FragmentIon_charge")) %>% filter(is.na(Remove)) %>% select(-Remove)
  }
  
  if(RemoveSharedTransitionsBetweenLightAndHeavy_Filter==TRUE){
    ## Remove shared transitions between light and heavy (b ions if the labeling is in C-ter)
    RemoveSharedTransitionsLightAndHeavy<-D1 %>% select(ProteinName,PeptideModifiedSequence,PrecursorCharge,ProductMz,ProductCharge,FragmentIon,IsotopeLabelType,IsDecoy) %>%
      distinct() %>%
      group_by(ProteinName,PeptideModifiedSequence,PrecursorCharge,ProductMz,IsDecoy) %>%
      tally %>%
      filter(n>1) %>% ungroup %>%
      select(ProteinName,PeptideModifiedSequence,PrecursorCharge,ProductMz,IsDecoy) %>%
      mutate(Remove="Remove")
    
    ## Filter less than MinimalInitialNumberOfTransitions_Filter transitions
    FilterLessThanNTrans<-D1 %>%
      select(ProteinName,PeptideModifiedSequence,PrecursorCharge,IsotopeLabelType,FragmentIon,ProductCharge,IsDecoy) %>%
      distinct()  %>%
      group_by(ProteinName,PeptideModifiedSequence,PrecursorCharge,IsotopeLabelType,IsDecoy) %>%
      tally() %>%
      mutate(Keep=ifelse(n>=MinimalInitialNumberOfTransitions_Filter,"keep","DoNOTKeep")) %>%
      select(ProteinName,PeptideModifiedSequence,PrecursorCharge,IsotopeLabelType,IsDecoy,Keep)
    
    
    D1<- D1 %>% left_join(RemoveSharedTransitionsLightAndHeavy,by=c("ProteinName","PeptideModifiedSequence","PrecursorCharge","ProductMz","IsDecoy")) %>%
      filter(is.na(Remove)) %>%
      select(-Remove) %>%
      left_join(data.frame(FilterLessThanNTrans),by = c("ProteinName","PeptideModifiedSequence", "PrecursorCharge", "IsotopeLabelType","IsDecoy")) %>%
      filter(Keep=="keep") %>% select(-Keep)%>%
      #select(ID_Analyte,IsotopeLabelType,ID_FragmentIon_charge,ID_Rep,InterpolatedTimes,
      #       InterpolatedIntensities,InterpolatedMassErrors,Area,LibraryIntensity,MinStartTime,MaxEndTime)%>%
      arrange(ID_Analyte,ID_FragmentIon_charge,ID_Rep)
  } else {
    
    ## Filter less than MinimalInitialNumberOfTransitions_Filter transitions
    FilterLessThanNTrans<-D1 %>%
      select(ProteinName,PeptideModifiedSequence,PrecursorCharge,IsotopeLabelType,FragmentIon,ProductCharge,IsDecoy) %>%
      distinct()  %>%
      group_by(ProteinName,PeptideModifiedSequence,PrecursorCharge,IsotopeLabelType,IsDecoy) %>%
      tally() %>%
      mutate(Keep=ifelse(n>=MinimalInitialNumberOfTransitions_Filter,"keep","DoNOTKeep")) %>%
      select(ProteinName,PeptideModifiedSequence,PrecursorCharge,IsotopeLabelType,IsDecoy,Keep)
    
    D1<- D1 %>%
      left_join(data.frame(FilterLessThanNTrans),by = c("ProteinName","PeptideModifiedSequence", "PrecursorCharge", "IsotopeLabelType","IsDecoy")) %>%
      filter(Keep=="keep") %>% select(-Keep)%>%
      #select(ID_Analyte,IsotopeLabelType,ID_FragmentIon_charge,ID_Rep,InterpolatedTimes,
      #       InterpolatedIntensities,InterpolatedMassErrors,Area,LibraryIntensity,MinStartTime,MaxEndTime)%>%
      arrange(ID_Analyte,ID_FragmentIon_charge,ID_Rep)
  }
  
  return(D1)}

#' MetaDataLibraryConcatenator
#'
#' This function is MetaDataLibraryConcatenator
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' MetaDataLibraryConcatenator()
MetaDataLibraryConcatenator<-function(MetaData1,MetaData2){
  M_Analyte<-full_join(MetaData1$M_Analyte,MetaData2$M_Analyte,by = c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge","IsDecoy"))%>%
    select(-ID_Analyte.y) %>% rename(ID_Analyte=ID_Analyte.x) %>% arrange(ID_Analyte)
  M_Analyte2<-M_Analyte %>% filter(is.na(ID_Analyte)) %>% distinct()
  if(dim(M_Analyte2)[1]>0){ M_Analyte2$ID_Analyte<-(max(unique(M_Analyte$ID_Analyte),na.rm = T)+1):  (max(unique(M_Analyte$ID_Analyte),na.rm = T)+length(M_Analyte2$ID_Analyte))
  M_Analyte <- M_Analyte %>% left_join(M_Analyte2, by = c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "IsDecoy")) %>%
    mutate(ID_Analyte.x=ifelse(is.na(ID_Analyte.x),ID_Analyte.y,ID_Analyte.x)) %>% 
    select(-ID_Analyte.y) %>% 
    rename(ID_Analyte=ID_Analyte.x)}
  
  
  M_Replicate<-full_join(MetaData1$M_Replicate,MetaData2$M_Replicate, by = c("FileName"))%>%
    select(-ID_Rep.y) %>% rename(ID_Rep=ID_Rep.x) %>% arrange(ID_Rep)
  M_Replicate2<- M_Replicate %>% filter(is.na(ID_Rep)) %>% 
    mutate(ID_Rep=(1:length(M_Replicate$ID_Rep))[-which( (1:length(M_Replicate$ID_Rep))  %in% c(M_Replicate$ID_Rep[!is.na(M_Replicate$ID_Rep)]))])
  M_Replicate<- M_Replicate %>% filter(!is.na(ID_Rep))
  M_Replicate<-rbind(M_Replicate,M_Replicate2)
  
  M_TransitionsType<-full_join(MetaData1$M_TransitionsType,MetaData2$M_TransitionsType, by = c("ProductCharge", "FragmentIon"))%>%
    select(-ID_FragmentIon_charge.y) %>% rename(ID_FragmentIon_charge=ID_FragmentIon_charge.x) %>% arrange(ID_FragmentIon_charge)
  M_TransitionsType2<- M_TransitionsType %>% filter(is.na(ID_FragmentIon_charge)) %>%
    mutate(ID_FragmentIon_charge=(1:length(M_TransitionsType$ID_FragmentIon_charge))[-which( (1:length(M_TransitionsType$ID_FragmentIon_charge))  %in% c(M_TransitionsType$ID_FragmentIon_charge[!is.na(M_TransitionsType$ID_FragmentIon_charge)]))])
  M_TransitionsType<- M_TransitionsType %>% filter(!is.na(ID_FragmentIon_charge))
  M_TransitionsType<-rbind(M_TransitionsType,M_TransitionsType2)
  
  M_Transitions<-full_join(MetaData1$M_Transitions,MetaData2$M_Transitions, by = c("ProteinName", "PeptideModifiedSequence", "ModifiedSequence", "PrecursorMz", "IsotopeLabelType", "IsDecoy", "PrecursorCharge", "ProductMz", "ProductCharge", "FragmentIon", "Quantitative","TransitionLocator")) %>%
    select(-ID_Analyte.x,-ID_FragmentIon_charge.x,-ID_Analyte.y,-ID_FragmentIon_charge.y) %>%
    left_join(M_TransitionsType, by = c("ProductCharge", "FragmentIon")) %>%
    left_join(M_Analyte,by = c("ProteinName", "PeptideModifiedSequence", "IsDecoy", "PrecursorCharge"))
  
  M_PrecursorResult<-rbind(MetaData1$M_PrecursorResult,MetaData2$M_PrecursorResult) %>% distinct()
  
  return(list(M_Analyte=M_Analyte,
              M_Replicate=M_Replicate,
              M_TransitionsType=M_TransitionsType,
              M_Transitions=M_Transitions,
              M_PrecursorResult=M_PrecursorResult))}

#' ReadFileInChunks_Format_And_Filter
#'
#' This function is ReadFileInChunks_Format_And_Filter
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' ReadFileInChunks_Format_And_Filter()
ReadFileInChunks_Format_And_Filter<-function(D.file.name,ColumnNames,PreviousMetaData,Skip,is.FirstIteration){
  # D<-data.frame(fread(file = D.file.name,header = F,stringsAsFactors = FALSE, na.strings = c("#N/A","NA","NaN"),skip=Skip,nrows = ReadNumberOfLines,col.names = ColumnNames,sep = ","))
  
  D<-data.frame(read_skyline_report_chunk(csv_file_path = D.file.name,
                                          skip_rows = Skip,
                                          n_rows = ReadNumberOfLines,
                                          column_names = ColumnNames))
  ChunkLastRowNumber=max(as.numeric(row.names(D)))
  Stop.Condition=!ChunkLastRowNumber==ReadNumberOfLines
  
  ChunkMetaData<-Indexing(D)
  #Correct MetaData
  if(is.FirstIteration==TRUE){
    NewMetaData<-ChunkMetaData} else{
      NewMetaData<-MetaDataLibraryConcatenator(PreviousMetaData,ChunkMetaData)}
  
  # Add MetaData to Chunk
  D<-D %>% left_join(NewMetaData$M_Analyte,by = c("ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "IsDecoy")) %>%
    left_join(NewMetaData$M_Replicate, by = "FileName") %>%
    left_join(NewMetaData$M_TransitionsType, by = c("ProductCharge", "FragmentIon"))
  
  # Find last Analyte, remove it and save line to include it in the next Chunk
  if(Stop.Condition){D<-D} else{
    Skip2<-which(D$ID_Analyte==D$ID_Analyte[length(D$ID_Analyte)])[1]-1
    D<-D[1:(Skip2),]
    Skip=Skip+Skip2}
  
  # Filter Data
  D<-Filter_Na_Shared_Or_LowMassTransitions(D)
  
  # Precursor Breaks
  PrecursorBreaks<-which(!duplicated(D$ID_Analyte))
  PrecursorBreaks[length(PrecursorBreaks)+1]<-length(D$ID_Analyte)+1
  
  return(list(Data=D,Skip=Skip,MetaData=NewMetaData,PrecursorBreaks=PrecursorBreaks,ChunkLastRowNumber=ChunkLastRowNumber))
}


#' AvG_writeParamsUsed
#'
#' This function is AvG_writeParamsUsed
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' AvG_writeParamsUsed()
AvG_writeParamsUsed<-function(SubSetOfParams){
  
  RefinementTAG<-if(RefinementWorkflow=="GlobalRefinement") {"GR"} else{
    if(RefinementWorkflow=="TransitionRefinement") {"TR"} else{
      if(RefinementWorkflow=="PeakBoundariesRefinement") {"PB"} else{
        if(RefinementWorkflow=="OnlyScoring") {"NO"}
      }}}
  
  
  
  if(SubSetOfParams=="DB"){
    write.csv(data.frame(Name_Tag,
                         RefinementTAG,
                         MinimalInitialNumberOfTransitions_Filter,
                         RemoveTransitionBelowOrdinal_Filter,
                         RemoveSharedTransitionsBetweenLightAndHeavy_Filter,
                         ReadNumberOfLines,
                         stringsAsFactors = F),
              file = paste0(Output.file,"/",gsub(gsub(paste("ParamsUsed",Name_Tag,"DB",Sys.time(), sep="_"),pattern = " ",replacement = "_"),pattern = ":",replacement = ""),".csv"),
              quote = F,row.names = F)}
  
  if(SubSetOfParams=="Optimization"){
    write.csv(data.frame(Name_Tag,
                         RefinementTAG,
                         NonZeroBaselineChromatogram,
                         MinimalNumberOfTransitionsAfterOptimization,
                         KeepPeptidesWithLowerNumberOfTransitions,
                         TopN_RankedbyIntensity,
                         MinimalNumberOfTransitionsAmongTopN,
                         alpha,
                         Beta,
                         SpectralLibraryDotProduct_limit,
                         MassError_Tolerance,
                         MassError_CutOff,
                         MinimalNumberOfConsecutivePoints,
                         MinimalIntensityPercentagePerTransition,
                         UseHeavyPeakBoundariesForLight,stringsAsFactors = F),
              file = paste0(Output.file,"/",gsub(gsub(paste("ParamsUsed",Name_Tag,RefinementTAG,Sys.time(), sep="_"),pattern = " ",replacement = "_"),pattern = ":",replacement = ""),".csv"),
              quote = F,row.names = F)}
}


### Transition refinement functions

#' dot.p
#'
#' This function is dot.p
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' dot.p()
dot.p<-function(u,v){
  c<-crossprod(as.numeric(u),as.numeric(v))/(sqrt(crossprod(as.numeric(u),as.numeric(u)))*sqrt(crossprod(as.numeric(v),as.numeric(v))))
  return(c)    
}

#' LinearEquationByRange
#'
#' This function is LinearEquationByRange
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' LinearEquationByRange()
LinearEquationByRange<-function(x,x_range,y_range){
  A <- matrix(data=c(x_range[1],1,x_range[2],1), nrow=2, ncol=2, byrow=TRUE)    
  b <- matrix(data=c(y_range[1],y_range[2]), nrow=2, ncol=1, byrow=FALSE)
  eq<-round(solve(A, b),digits = 3)
  x<-eq[1]*x+eq[2]
  return(x)
}

#' ScoreTransform.dotP
#'
#' This function is ScoreTransform.dotP
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' ScoreTransform.dotP()
ScoreTransform.dotP<-function(x){
  x=ifelse(as.numeric(x)<=SpectralLibraryDotProduct_limit*0.95,LinearEquationByRange(x,x_range = c(0,SpectralLibraryDotProduct_limit*0.95),y_range = c(0,0.1)),
           ifelse(as.numeric(x)<=SpectralLibraryDotProduct_limit,LinearEquationByRange(x,x_range = c(SpectralLibraryDotProduct_limit*0.95,SpectralLibraryDotProduct_limit),y_range = c(0.1,0.85)),
                  ifelse(as.numeric(x)<=(1+SpectralLibraryDotProduct_limit)/2,LinearEquationByRange(x,x_range = c(SpectralLibraryDotProduct_limit,(1+SpectralLibraryDotProduct_limit)/2),y_range = c(0.85,1)),
                         LinearEquationByRange(x,x_range = c((1+SpectralLibraryDotProduct_limit)/2,1),y_range = c(1,1)))))
  return(x)}

#' Rank.filter
#'
#' This function is Rank.filter
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Rank.filter()
Rank.filter=function(y,TopN_RankedbyIntensity,MinimalNumberOfTransitionsAmongTopN,Transition.Rank){
  Eval.Rank=Transition.Rank[which(y==1)]
  length(Eval.Rank[Eval.Rank<=TopN_RankedbyIntensity])>=MinimalNumberOfTransitionsAmongTopN}

#' FindTranswithConsecutivePoints
#'
#' This function is FindTranswithConsecutivePoints
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' FindTranswithConsecutivePoints()
FindTranswithConsecutivePoints<-function(Norm.Chrom,MinimalIntensityPercentagePerTransition,MinimalNumberOfConsecutivePoints){
  L0=apply(Norm.Chrom,2,function(x){
    x<-ifelse(x>MinimalIntensityPercentagePerTransition,1,0)})
  
  L1=apply(Norm.Chrom,2,function(x){ #Counts how many consecutive points are about 4% of the max and then tells if the transitions has more than 5 points, i.e. 1 if ok 0 if not
    x<-ifelse(is.na(x),0,x)
    x<-ifelse(x>MinimalIntensityPercentagePerTransition,1,0)
    y<-rle(as.numeric(x))
    z=ifelse(max(y$lengths[which(y$values==1)])>=MinimalNumberOfConsecutivePoints,1,0)})
  M0=sum(L1) # Num of transitions with more than 5 consecutive points with more than 4% of the max Normalized intensity
  if(M0<=1){return(0)}
  
  QQ<-combn(colnames(Norm.Chrom[,which(as.numeric(L1)>=1)]),2)
  L2=apply(QQ,2,function(x){
    l=L0[,x]
    l=l[,1]*l[,2]})
  
  colnames(L2)=apply(QQ,2,function(x){
    paste0(x,collapse = "_")})
  L2
  
  L3=apply(L2,2,function(x){
    y<-rle(as.numeric(x))
    z=ifelse(max(y$lengths[which(y$values==1)])>=MinimalNumberOfConsecutivePoints,1,0)})
  M=sum(L3,na.rm = T)
  return(ifelse(M>=1,1,0))}

#' Similarity.Score_informsMPRA
#'
#' This function is Similarity.Score_informsMPRA
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Similarity.Score_informsMPRA()
Similarity.Score_informsMPRA<- function(Norm.Chrom, y){
  
  list2 <- lapply(Norm.Chrom, function(m){
    l<-FindTranswithConsecutivePoints(Norm.Chrom = m,MinimalIntensityPercentagePerTransition,MinimalNumberOfConsecutivePoints)
    m[, which(y==1)]*l
  })
  A=sapply(list2, function(m){
    m=data.matrix(m)
    m.mean <- apply(m, 1, mean,na.rm=T)
    m.mean<-m.mean/max(m.mean,na.rm = T)*100
    m=cbind(m.mean, m)
    m=apply(m, 2, function(x2) {dot.p(x2, m[,1])})
    m[is.na(m)]<-0
    m=m[-c(1)]
    return(m)
  })
  Score=mean(A)*mean(apply(A,2,FUN = function(x){ x=ifelse(is.na(as.numeric(x)),0,as.numeric(x));
  ifelse(max(x,na.rm = T)==0,x,x[!x==max(x,na.rm = T)])}))
  
  # Best_reps=substr(names(sort(apply(A,2,mean),decreasing = T)[1:floor((dim(A)[2])*0.25)]),
  #                  start = 1,
  #                  stop = str_locate(names(sort(apply(A,2,mean),decreasing = T)[1:floor((dim(A)[2])*0.25)]),pattern = " ")-1)
  # 
  Best_reps=names(sort(apply(A,2,mean),decreasing = T)[1:floor((dim(A)[2])*0.25)])
  return(list(Score=Score,Best_reps=Best_reps))}

#' Similarity.Score.Report_informsMPRA
#'
#' This function is Similarity.Score.Report_informsMPRA
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Similarity.Score.Report_informsMPRA()
Similarity.Score.Report_informsMPRA<- function(Norm.Chrom, y){
  list2 <- lapply(Norm.Chrom, function(m){
    l<-FindTranswithConsecutivePoints(Norm.Chrom = m,MinimalIntensityPercentagePerTransition,MinimalNumberOfConsecutivePoints)
    m[, which(y==1)]*l
  })
  A=sapply(list2, function(m){
    m=data.matrix(m)
    m.mean <- apply(m, 1, mean,na.rm=T)
    m.mean<-m.mean/max(m.mean,na.rm = T)*100
    m=cbind(m.mean, m)
    m=apply(m, 2, function(x2) {dot.p(x2, m[,1])})
    m[is.na(m)]<-0
    m=m[-c(1)]
    return(m)
  })
  Report=apply(A,2,FUN = function(x){mean(x[!x==max(x,na.rm = T)])*x})
  Report<-apply(Report,2,function(x){ifelse(is.nan(x),0,x)})
  # Best_reps=substr(names(sort(apply(A,2,mean),decreasing = T)[1:floor((dim(A)[2])*0.25)]),
  #                  start = 1,
  #                  stop = str_locate(names(sort(apply(A,2,mean),decreasing = T)[1:floor((dim(A)[2])*0.25)]),pattern = " ")-1)
  # 
  Best_reps=names(sort(apply(A,2,mean),decreasing = T)[1:floor((dim(A)[2])*0.25)])
  return(list(Report=Report,Best_reps=Best_reps))}

#' Similarity.Score.Report.AllTransitions
#'
#' This function is Similarity.Score.Report.AllTransitions
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Similarity.Score.Report.AllTransitions()
Similarity.Score.Report.AllTransitions<- function(Norm.ChromX, y){
  list2 <- lapply(Norm.ChromX, function(m){
    m[, which(y==1)]
  })
  A=sapply(list2, function(m){
    m=data.matrix(m)
    m.mean <- apply(m, 1, mean,na.rm=T)
    m.mean<-m.mean/max(m.mean,na.rm = T)*100
    #m=cbind(m.mean, m)
    #m=apply(m, 2, function(x2) {dot.p(x2, m[,1])})
    #m[is.na(m)]<-0
    #m=m[-c(1)]
    return(m.mean)
  })
  
  B=Map(cbind, A,Norm.ChromX)
  B=sapply(B, function(m){
    #m=data.matrix(m)
    #m.mean <- apply(m, 1, mean,na.rm=T)
    #m.mean<-m.mean/max(m.mean,na.rm = T)*100
    #m=cbind(m.mean, m)
    m=apply(m, 2, function(x2) {dot.p(x2, m[,1])})
    m[is.na(m)]<-0
    m=m[-c(1)]
    return(m)
  })
  
  return(apply(B,2,FUN = function(x){mean(x[!x==max(x,na.rm = T)])*x}))}

#' MPRA.Score.dotp
#'
#' This function is MPRA.Score.dotp
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' MPRA.Score.dotp()
MPRA.Score.dotp <- function(list){
  mean(
    sapply(list, function(m){
      m=apply(m, 1, function(x2) {dot.p(x2, m[1,])})
      mean(m[-c(1,2)],na.rm = T)
    }),
    na.rm = T)
}

#' Library.dotp
#'
#' This function is Library.dotp
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Library.dotp()
Library.dotp <- function(list){
  mean(
    sapply(list, function(m){
      m=apply(m, 1, function(x2) {dot.p(x2, m[2,])})
      mean(m[-c(1,2)],na.rm = T)
    }),
    na.rm = T)
}

#' MPRA_y_SpectLib.Fitness
#'
#' This function is MPRA_y_SpectLib.Fitness
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' MPRA_y_SpectLib.Fitness()
MPRA_y_SpectLib.Fitness<- function(Transition.Area, y){
  list2<- lapply(Transition.Area, function(m){
    m=m[,-dim(m)[2]]
    m=m[, which(y==1)]
    m= m/apply(m,1,function(x){sum(x,na.rm = T)})*100
  })
  MPRA.score=MPRA.Score.dotp(list2)
  Library.dotp=Library.dotp(list2)
  return(data.frame(MPRA.score=MPRA.score,Library.dotp=Library.dotp))
}

#' MPRA_y_SpectLib.Report
#'
#' This function is MPRA_y_SpectLib.Report
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' MPRA_y_SpectLib.Report()
MPRA_y_SpectLib.Report<- function(Transition.Area, y){
  list2<- lapply(Transition.Area, function(m){
    m=m[,-dim(m)[2]]
    m=m[, which(y==1)]
    m= m/apply(m,1,sum)*100
  })
  
  MPRA.Score.dotp.report<- function(list2){
    sapply(list2, function(m){
      m=apply(m, 1, function(x2) {dot.p(x2, m[1,])})
      m[-c(1,2)]
    })
  }
  
  Library.dotp.report<- function(list2){
    sapply(list2, function(m){
      m=apply(m, 1, function(x2) {dot.p(x2, m[2,])})
      m[-c(1,2)]
    })
  }
  
  return(list(MPRA.score=MPRA.Score.dotp.report(list2),Library.dotp=Library.dotp.report(list2)))
}

#' MPRA_y_SpectLib.Fitness_informed
#'
#' This function is MPRA_y_SpectLib.Fitness_informed
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' MPRA_y_SpectLib.Fitness_informed()
MPRA_y_SpectLib.Fitness_informed<- function(Transition.Area, y,Best_reps,Best_reps_MassErrors){
  list2<- lapply(Transition.Area, function(m){
    m=m[,-dim(m)[2]]
    m=m[, which(y==1)]
    m= m/apply(m,1,function(x){sum(x,na.rm = T)})*100
  })
  MPRA.score=MPRA.Score.dotp(list2)
  
  Best_reps2=intersect(Best_reps,Best_reps_MassErrors)
  if(length(Best_reps2)==0){Best_reps2=unique(c(Best_reps[1:floor(length(Best_reps)/2)],Best_reps_MassErrors[1:floor(length(Best_reps_MassErrors)/2)]))}
  
  Best.Areas=lapply(Transition.Area, function(m){
    m=m[-c(1,2),-dim(m)[2]]
  })
  
  Best.Areas<-apply(data.frame(do.call(rbind,Best.Areas)) %>%
                      mutate(Rep=row.names(.)) %>%
                      mutate(Rep=paste0(Rep,".")) %>% 
                      mutate(Rep=gsub(Rep,pattern=" ",replacement=".")) %>% 
                      mutate(Rep2=substr(Rep,start = str_locate(Rep,pattern = "Rep_"),stop = nchar(Rep)),
                             Analyte=substr(Rep,start = str_locate(Rep,pattern = "Analyte_"),stop = nchar(Rep)),
                             IsotopeLabelType=substr(Rep,start = str_locate(Rep,pattern = "IsotopeLabelType_"),stop = nchar(Rep))) %>%
                      mutate(Rep2=substr(Rep2,stop = str_locate(Rep2,pattern = "\\.")-1,start = 1),
                             Analyte=substr(Analyte,stop = str_locate(Analyte,pattern = "\\.")-1,start = 1),
                             IsotopeLabelType=substr(IsotopeLabelType,stop = str_locate(IsotopeLabelType,pattern = "\\.")-1,start = 1)) %>%
                      mutate(Reps=paste(Rep2,Analyte,IsotopeLabelType,sep=" "))%>% select(-Rep,-Rep2,-Analyte,-IsotopeLabelType) %>% filter(Reps %in% Best_reps2) %>% select(-Reps)
                    ,2,function(x){sum(as.numeric(x))})
  Best.Areas[length(Best.Areas)+1]<-0
  
  Transition.Area2=lapply(Transition.Area,function(L){
    L[1,]=Best.Areas
    L
  })
  
  list2<- lapply(Transition.Area2, function(m){
    m=m[,-dim(m)[2]]
    m=m[, which(y==1)]
    m= m/apply(m,1,function(x){sum(x,na.rm = T)})*100
  })
  MPRA.score_informed=MPRA.Score.dotp(list2)
  Library.dotp=Library.dotp(list2)
  return(data.frame(MPRA.score=MPRA.score*MPRA.score_informed,Library.dotp=Library.dotp))
}

#' MPRA_y_SpectLib.Report_informed
#'
#' This function is MPRA_y_SpectLib.Report_informed
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' MPRA_y_SpectLib.Report_informed()
MPRA_y_SpectLib.Report_informed<- function(Transition.Area, y,Best_reps,Best_reps_MassErrors){
  list<- lapply(Transition.Area, function(m){
    m=m[,-dim(m)[2]]
    m=m[, which(y==1)]
    m= m/apply(m,1,sum)*100
  })
  
  MPRA.Score.dotp.report<- function(list){
    sapply(list, function(m){
      m=apply(m, 1, function(x2) {dot.p(x2, m[1,])})
      m[-c(1,2)]
    })
  }
  
  
  
  Best_reps2=intersect(Best_reps,Best_reps_MassErrors)
  if(length(Best_reps2)==0){Best_reps2=unique(c(Best_reps[1:floor(length(Best_reps)/2)],Best_reps_MassErrors[1:floor(length(Best_reps_MassErrors)/2)]))}
  
  Best.Areas=lapply(Transition.Area, function(m){
    m=m[-c(1,2),-dim(m)[2]]
  })
  
  
  Best.Areas=lapply(Transition.Area, function(m){
    m=m[-c(1,2),-dim(m)[2]]
  })
  
  Best.Areas<-apply(data.frame(do.call(rbind,Best.Areas)) %>%
                      mutate(Rep=row.names(.)) %>%
                      mutate(Rep=paste0(Rep,".")) %>% 
                      mutate(Rep=gsub(Rep,pattern=" ",replacement=".")) %>% 
                      mutate(Rep2=substr(Rep,start = str_locate(Rep,pattern = "Rep_"),stop = nchar(Rep)),
                             Analyte=substr(Rep,start = str_locate(Rep,pattern = "Analyte_"),stop = nchar(Rep)),
                             IsotopeLabelType=substr(Rep,start = str_locate(Rep,pattern = "IsotopeLabelType_"),stop = nchar(Rep))) %>%
                      mutate(Rep2=substr(Rep2,stop = str_locate(Rep2,pattern = "\\.")-1,start = 1),
                             Analyte=substr(Analyte,stop = str_locate(Analyte,pattern = "\\.")-1,start = 1),
                             IsotopeLabelType=substr(IsotopeLabelType,stop = str_locate(IsotopeLabelType,pattern = "\\.")-1,start = 1)) %>%
                      mutate(Reps=paste(Rep2,Analyte,IsotopeLabelType,sep=" "))%>% select(-Rep,-Rep2,-Analyte,-IsotopeLabelType) %>% 
                      filter(Reps %in% Best_reps2) %>%
                      select(-Reps)
                    ,2,function(x){sum(as.numeric(x))})
  Best.Areas[length(Best.Areas)+1]<-0
  
  Transition.Area2=lapply(Transition.Area,function(L){
    L[1,]=Best.Areas
    L
  })
  
  list2<- lapply(Transition.Area2, function(m){
    m=m[,-dim(m)[2]]
    m=m[, which(y==1)]
    m= m/apply(m,1,function(x){sum(x,na.rm = T)})*100
  })
  MPRA.Score.dotp.report_informed<- function(list2){
    sapply(list2, function(m){
      m=apply(m, 1, function(x2) {dot.p(x2, m[1,])})
      m[-c(1,2)]
    })
  }
  
  Library.dotp.report<- function(list2){
    sapply(list2, function(m){
      m=apply(m, 1, function(x2) {dot.p(x2, m[2,])})
      m[-c(1,2)]
    })
  }
  
  
  
  
  return(list(MPRA.score=MPRA.Score.dotp.report(list)*MPRA.Score.dotp.report_informed(list2),
              Library.dotp=Library.dotp.report(list2)))
}

#' Intensity.Fitness
#'
#' This function is Intensity.Fitness
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Intensity.Fitness()
Intensity.Fitness<- function(Transition.Area, y){
  Intensity.Score<-lapply(Transition.Area, function(m){
    Total.To.Max.Ratio=m[-c(1,2),-dim(m)[2]]
    Total.To.Max.Ratio.Initial=m[-c(1,2),dim(m)[2]]
    Total.To.Max.Ratio=Total.To.Max.Ratio[, which(y==1)]
    MAX.L=apply(Total.To.Max.Ratio,1,max,na.rm =T)
    SUM.L=apply(Total.To.Max.Ratio,1,sum,na.rm =T)
    #Total.To.Max.Ratio=cbind(Total.To.Max.Ratio,Total.To.Max.Ratio=SUM.L/MAX.L)
    Total.To.Max.Ratio=Total.To.Max.Ratio=SUM.L/MAX.L
    #cbind(Total.To.Max.Ratio,Total.To.Max.Ratio.Initial)
    Intensity.Score=Total.To.Max.Ratio/Total.To.Max.Ratio.Initial
    mean(Intensity.Score,na.rm=T)
  })
  Intensity.Score<-mean(unlist(Intensity.Score))
  
  return(data.frame(Intensity.Score=Intensity.Score))
}

#' Intensity.Fitness.Report
#'
#' This function is Intensity.Fitness.Report
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Intensity.Fitness.Report()
Intensity.Fitness.Report<- function(Transition.Area, y){
  Intensity.Score<-lapply(Transition.Area, function(m){
    Total.To.Max.Ratio=m[-c(1,2),-dim(m)[2]]
    Total.To.Max.Ratio.Initial=m[-c(1,2),dim(m)[2]]
    Total.To.Max.Ratio=Total.To.Max.Ratio[, which(y==1)]
    MAX.L=apply(Total.To.Max.Ratio,1,max,na.rm =T)
    SUM.L=apply(Total.To.Max.Ratio,1,sum,na.rm =T)
    #Total.To.Max.Ratio=cbind(Total.To.Max.Ratio,Total.To.Max.Ratio=SUM.L/MAX.L)
    Total.To.Max.Ratio=Total.To.Max.Ratio=SUM.L/MAX.L
    #cbind(Total.To.Max.Ratio,Total.To.Max.Ratio.Initial)
    Intensity.Score=Total.To.Max.Ratio/Total.To.Max.Ratio.Initial
    Intensity.Score
  })
  Intensity.Score<-unlist(Intensity.Score)
  
  return(data.frame(Intensity.Score=Intensity.Score))
}

#' MassError.Score.Fitness_informsMPRA
#'
#' This function is MassError.Score.Fitness_informsMPRA
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' MassError.Score.Fitness_informsMPRA()
MassError.Score.Fitness_informsMPRA<- function(MassErrors, y){
  list2 <- lapply(MassErrors, function(m){
    m[, which(y==1)]
  })
  A=lapply(list2,function(L){
    #if(any(as.numeric(L["SkorMassErrors.IntegratedZone",])==0)) {return(0)} else{
    mean(as.numeric(L["SkorMassErrors.IntegratedZone",]),na.rm = T)
    #}
  })
  A=do.call(rbind,A)
  Score=mean(A,na.rm = T)
  
  Best_reps= (data.frame(A) %>% mutate(Reps=row.names(.)) %>%
                top_n(A,n = floor(dim(.)[1]/4)) %>% select(Reps))[,1]
  return(list(Score=Score,Best_reps_MassErrors=Best_reps))}

#' MassError.Score.Report
#'
#' This function is MassError.Score.Report
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' MassError.Score.Report()
MassError.Score.Report<- function(MassErrors, y){
  list2 <- lapply(MassErrors, function(m){
    m[, which(y==1)]
  })
  A=lapply(list2,function(L){
    mean(as.numeric(L["SkorMassErrors.IntegratedZone",]),na.rm = T)
  })
  A=do.call(rbind,A)
  
  Best_reps= (data.frame(A) %>% mutate(Reps=row.names(.)) %>%
                top_n(A,n = floor(dim(.)[1]/4)) %>% select(Reps))[,1]
  return(list(Report=A,Best_reps_MassErrors=Best_reps))}

#' Transition.Remover
#'
#' This function is Transition.Remover
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Transition.Remover()
Transition.Remover<- function(Matrix,y){
  Matrix2<- lapply(Matrix, function(m){
    m[, which(y==1)]
  })
  return(Matrix2)
}

#' Transition.Classifier
#'
#' This function is Transition.Classifier
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Transition.Classifier()
Transition.Classifier<-function(Norm.Chrom,Transition.Area,num_trans,Transition.Rank){
  U<-Similarity.Score.Report_informsMPRA(Norm.Chrom = Norm.Chrom ,y = rep(1,num_trans))
  U<-U$Report
  U2<-data.frame(U) %>% 
    mutate(ID_FragmentIon_charge=row.names(U)) %>%
    gather(key = "Rep",value = "Similarity.Score",1:(dim(U)[2])) %>% 
    ungroup %>% 
    group_by(ID_FragmentIon_charge) %>% 
    summarise(median.Initial.SimilarityScore=median(Similarity.Score,na.rm = T ),
              The25quartile.Initial.SimilarityScore=quantile(Similarity.Score,0.25,na.rm = T)
    ) %>% 
    mutate(Rank.median.Initial.SimilarityScore=rank(-median.Initial.SimilarityScore,ties.method = "min"),
           Rank.The25quartile.Initial.SimilarityScore=rank(-The25quartile.Initial.SimilarityScore,ties.method = "min")) %>%
    select(-median.Initial.SimilarityScore,-The25quartile.Initial.SimilarityScore )
  
  Q<-do.call(Transition.Area,what = rbind)
  Q<-Q[-c(1,2),]
  Q<-data.frame(Q) %>% select(-Total.To.Max.Ratio) %>%
    mutate(Names=row.names(Q)) %>%
    gather(key = "ID_FragmentIon_charge",value = "Area",1:(dim(Q)[2]-1)) %>%
    ungroup %>%
    group_by(ID_FragmentIon_charge) %>% 
    summarise(median.Initial.Area=median(Area,na.rm = T ), The25quartile.Initial.Area=quantile(Area,0.25,na.rm = T)) %>% 
    mutate(Rank.median.Initial.Area=rank(-median.Initial.Area,ties.method = "min"),
           Rank.The25quartile.Initial.Area=rank(-The25quartile.Initial.Area,ties.method = "min")) %>%
    select(-median.Initial.Area,-The25quartile.Initial.Area) %>% mutate(ID_FragmentIon_charge=gsub(ID_FragmentIon_charge,pattern="X",replacement = ""))
  
  U2<-U2 %>% left_join(Q,by = "ID_FragmentIon_charge")
  U4<-data.frame(t(Transition.Rank)) %>% mutate(ID_FragmentIon_charge=row.names(.))
  U2<-U2 %>% left_join(U4,by = "ID_FragmentIon_charge")
  return(U2)
}

#' Transition.Classifier2
#'
#' This function is Transition.Classifier2
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Transition.Classifier2()
Transition.Classifier2<-function(Norm.Chrom,Transition.Area,num_trans,y,Transition.Rank){
  U<-Similarity.Score.Report.AllTransitions(Norm.Chrom = Norm.Chrom ,y = y)
  U2<-data.frame(U) %>% 
    mutate(ID_FragmentIon_charge=row.names(U)) %>%
    gather(key = "Rep",value = "Similarity.Score",1:(dim(U)[2])) %>% 
    ungroup %>% 
    group_by(ID_FragmentIon_charge) %>% 
    summarise(median.Initial.SimilarityScore=median(Similarity.Score,na.rm = T ),
              The25quartile.Initial.SimilarityScore=quantile(Similarity.Score,0.25,na.rm = T)
    ) %>% 
    mutate(Rank.median.Initial.SimilarityScore=rank(-median.Initial.SimilarityScore,ties.method = "min"),
           Rank.The25quartile.Initial.SimilarityScore=rank(-The25quartile.Initial.SimilarityScore,ties.method = "min")) %>%
    select(-median.Initial.SimilarityScore,-The25quartile.Initial.SimilarityScore )
  
  Q<-do.call(Transition.Area,what = rbind)
  Q<-Q[-c(1,2),]
  Q<-data.frame(Q) %>% select(-Total.To.Max.Ratio) %>%
    mutate(Names=row.names(Q)) %>%
    gather(key = "ID_FragmentIon_charge",value = "Area",1:(dim(Q)[2]-1)) %>%
    ungroup %>%
    group_by(ID_FragmentIon_charge) %>% 
    summarise(median.Initial.Area=median(Area,na.rm = T ), The25quartile.Initial.Area=quantile(Area,0.25,na.rm = T)) %>% 
    mutate(Rank.median.Initial.Area=rank(-median.Initial.Area,ties.method = "min"),
           Rank.The25quartile.Initial.Area=rank(-The25quartile.Initial.Area,ties.method = "min")) %>%
    select(-median.Initial.Area,-The25quartile.Initial.Area) %>% mutate(ID_FragmentIon_charge=gsub(ID_FragmentIon_charge,pattern="X",replacement = ""))
  
  U2<-U2 %>% left_join(Q,by = "ID_FragmentIon_charge")
  U4<-data.frame(t(Transition.Rank)) %>% mutate(ID_FragmentIon_charge=row.names(.))
  U2<-U2 %>% left_join(U4,by = "ID_FragmentIon_charge")
  return(U2)
}

#' Sug.Matrix.FUN
#'
#' This function is Sug.Matrix.FUN
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Sug.Matrix.FUN()
Sug.Matrix.FUN<-function(Trans.Classified,Trans.VectorX){
  K=list()
  K[[1]]=data.frame(Trans.Classified %>% top_n(Rank.median.Initial.SimilarityScore,n = -6) %>% arrange(Rank.median.Initial.SimilarityScore))[,1]
  K[[2]]=K[[1]][1:5]
  K[[3]]=K[[1]][1:4]
  
  K[[4]]=data.frame(Trans.Classified %>% top_n(Rank.The25quartile.Initial.SimilarityScore,n = -6) %>% arrange(Rank.The25quartile.Initial.SimilarityScore))[,1]
  K[[5]]=K[[4]][1:5]
  K[[6]]=K[[4]][1:4]
  
  K[[7]]=data.frame(Trans.Classified %>% top_n(Rank.median.Initial.SimilarityScore,n = -5) %>% top_n(Rank.median.Initial.Area,n = -4))[,1]
  K[[8]]=data.frame(Trans.Classified %>% top_n(Rank.median.Initial.SimilarityScore,n = -6) %>% top_n(Rank.median.Initial.Area,n = -5))[,1]
  K[[9]]=data.frame(Trans.Classified %>% top_n(Rank.median.Initial.SimilarityScore,n = -6) %>% top_n(Rank.median.Initial.Area,n = -4))[,1]
  K[[10]]=data.frame(Trans.Classified %>% top_n(Rank.median.Initial.SimilarityScore,n = -7) %>% top_n(Rank.median.Initial.Area,n = -5))[,1]
  K[[11]]=data.frame(Trans.Classified %>% top_n(Rank.median.Initial.SimilarityScore,n = -7) %>% top_n(Rank.median.Initial.Area,n = -4))[,1]
  
  K[[12]]=data.frame(Trans.Classified %>% top_n(Rank.The25quartile.Initial.SimilarityScore,n = -5) %>% top_n(Rank.median.Initial.Area,n = -4))[,1]
  K[[13]]=data.frame(Trans.Classified %>% top_n(Rank.The25quartile.Initial.SimilarityScore,n = -6) %>% top_n(Rank.median.Initial.Area,n = -5))[,1]
  K[[14]]=data.frame(Trans.Classified %>% top_n(Rank.The25quartile.Initial.SimilarityScore,n = -6) %>% top_n(Rank.median.Initial.Area,n = -4))[,1]
  K[[15]]=data.frame(Trans.Classified %>% top_n(Rank.The25quartile.Initial.SimilarityScore,n = -7) %>% top_n(Rank.median.Initial.Area,n = -5))[,1]
  K[[16]]=data.frame(Trans.Classified %>% top_n(Rank.The25quartile.Initial.SimilarityScore,n = -7) %>% top_n(Rank.median.Initial.Area,n = -4))[,1]
  
  K[[17]]=data.frame(Trans.Classified %>% top_n(Rank.median.Initial.Area,n = -4))[,1]
  K[[18]]=Trans.VectorX
  
  sug.matrix=list()
  for(k in 1:18){
    sug.matrix[[k]]=sapply(Trans.VectorX,function(x) ifelse(x %in% K[[k]],1,0))
  }
  sug.matrix<-unique(do.call("rbind", sug.matrix))
}

#' GA.Fitness
#'
#' This function is GA.Fitness
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' GA.Fitness()
GA.Fitness<-function(y,Norm.Chrom,Transition.Area,Transition.Rank,MassErrors){
  if(sum(y) < MinimalNumberOfTransitionsAfterOptimization) return(0)
  if(Rank.filter(y,TopN_RankedbyIntensity,MinimalNumberOfTransitionsAmongTopN,Transition.Rank)==FALSE) return(0)
  
  SimS<-Similarity.Score_informsMPRA(Norm.Chrom,y)
  MassS<-MassError.Score.Fitness_informsMPRA(MassErrors,y)
  A=MPRA_y_SpectLib.Fitness_informed(Transition.Area,y,Best_reps = SimS$Best_reps,Best_reps_MassErrors = MassS$Best_reps_MassErrors)
  if(is.na(A[,"Library.dotp"])) return(0)
  #if(A[,"Library.dotp"]<SpectralLibraryDotProduct_limit) return(0)
  
  Scores=data.frame(MPRA.score=A[,"MPRA.score"],
                    Sim.score=SimS$Score,
                    Intensity.Score=Intensity.Fitness(Transition.Area,y),
                    MassError.score=MassS$Score,
                    Dotp.TransformedScore=ScoreTransform.dotP(A[,"Library.dotp"])) %>% 
    mutate(Score=(1-alpha-Beta)*(MPRA.score+Sim.score)+alpha*Intensity.Score+Beta*MassError.score+0.05*Dotp.TransformedScore)
  
  return(as.numeric(Scores$Score))
}

#' num.trans.To.Remove.Fun
#'
#' This function is num.trans.To.Remove.Fun
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' num.trans.To.Remove.Fun()
num.trans.To.Remove.Fun<-function(num_transX){
  if(num_transX<=5) {N=num_transX} else
    if(num_transX<=12) {N=ifelse((num_transX-1)<5,5,(num_transX-1))} else
      if(num_transX<20) {N=ifelse(round(num_transX*0.75,digits = 0)<10,5,round(num_transX*0.75,digits = 0))} else 
      {N=ifelse(round(num_transX*0.6,digits = 0)<5,5,round(num_transX*0.6,digits = 0))}
  return(num_transX-N)
}

#' Report.Replicate_informedMPRA
#'
#' This function is Report.Replicate_informedMPRA
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Report.Replicate_informedMPRA()
Report.Replicate_informedMPRA<-function(Norm.ChromX,Transition.AreaX,MassErrorsX,y,Comment){
  #Similarity Scores  
  H1<-Similarity.Score.Report_informsMPRA(Norm.ChromX,y)
  H=H1$Report
  H<-data.frame(H)%>%
    mutate(ID_FragmentIon_charge=row.names(H)) %>%
    gather(key = "Rep",value = "Similarity.Score",1:(dim(H)[2])) %>%
    group_by(Rep) %>% summarise(Similarity.Score=mean(Similarity.Score))
  H<-data.frame(H) %>% mutate(Rep=paste0(Rep,".")) %>% 
    mutate(Rep2=substr(Rep,start = str_locate(Rep,pattern = "Rep_"),stop = nchar(Rep)),
           Analyte=substr(Rep,start = str_locate(Rep,pattern = "Analyte_"),stop = nchar(Rep)),
           IsotopeLabelType=substr(Rep,start = str_locate(Rep,pattern = "IsotopeLabelType_"),stop = nchar(Rep))) %>%
    mutate(Rep2=substr(Rep2,stop = str_locate(Rep2,pattern = "\\.")-1,start = 5),
           Analyte=substr(Analyte,stop = str_locate(Analyte,pattern = "\\.")-1,start = 9),
           IsotopeLabelType=substr(IsotopeLabelType,stop = str_locate(IsotopeLabelType,pattern = "\\.")-1,start = 18)) %>%
    rename(ID_Rep=Rep2, ID_Analyte=Analyte) %>% select(-Rep)
  
  K0<-MassError.Score.Report(MassErrors = MassErrorsX,y = y)
  K1=as.character(row.names(K0$Report))
  K<- data.frame(K0$Report) %>% mutate(Rep=K1) %>% 
    mutate(Rep=paste0(Rep,".")) %>% 
    mutate(Rep=gsub(Rep,pattern=" ",replacement=".")) %>% 
    mutate(Rep2=substr(Rep,start = str_locate(Rep,pattern = "Rep_"),stop = nchar(Rep)),
           Analyte=substr(Rep,start = str_locate(Rep,pattern = "Analyte_"),stop = nchar(Rep)),
           IsotopeLabelType=substr(Rep,start = str_locate(Rep,pattern = "IsotopeLabelType_"),stop = nchar(Rep))) %>%
    mutate(Rep2=substr(Rep2,stop = str_locate(Rep2,pattern = "\\.")-1,start = 5),
           Analyte=substr(Analyte,stop = str_locate(Analyte,pattern = "\\.")-1,start = 9),
           IsotopeLabelType=substr(IsotopeLabelType,stop = str_locate(IsotopeLabelType,pattern = "\\.")-1,start = 18)) %>%
    rename(ID_Rep=Rep2, ID_Analyte=Analyte) %>% select(-Rep) %>% rename(Score.MassError=K0.Report)
  
  
  
  
  I<-MPRA_y_SpectLib.Report_informed(Transition.Area = Transition.AreaX,y = y,Best_reps = H1$Best_reps,K0$Best_reps_MassErrors)
  
  MPRA<-data.frame(Rep=row.names(I$MPRA.score),I$MPRA.score)
  MPRA<-MPRA%>% gather(key = "Rep2",value = "MPRA.Score",2:dim(MPRA)[2]) %>%
    mutate(Rep=paste0(Rep,".",Rep2)) %>% select(-Rep2)
  
  Lib.dotp<-data.frame(Rep=row.names(I$Library.dotp),I$Library.dotp)
  Lib.dotp<-Lib.dotp%>% gather(key = "Rep2",value = "Library.dotp",2:dim(Lib.dotp)[2]) %>%
    mutate(Rep=paste0(Rep,".",Rep2)) %>% select(-Rep2)
  
  I<-left_join(MPRA,Lib.dotp,by="Rep") %>% mutate(Rep=paste0(Rep,".")) %>% mutate(Rep=substr(Rep,start = str_locate(Rep,pattern = "Rep_"),stop = nchar(Rep)),
                                                                                  Analyte=substr(Rep,start = str_locate(Rep,pattern = "Analyte_"),stop = nchar(Rep)),
                                                                                  IsotopeLabelType=substr(Rep,start = str_locate(Rep,pattern = "IsotopeLabelType_"),stop = nchar(Rep))
  ) %>%
    mutate(Rep=substr(Rep,stop = str_locate(Rep,pattern = "\\.")-1,start = 5),
           Analyte=substr(Analyte,stop = str_locate(Analyte,pattern = "\\.")-1,start = 9),
           IsotopeLabelType=substr(IsotopeLabelType,stop = str_locate(IsotopeLabelType,pattern = "\\.")-1,start = 18)) %>%
    rename(ID_Rep=Rep, ID_Analyte=Analyte)
  
  
  
  J<-Intensity.Fitness.Report(Transition.Area = Transition.AreaX,y=y)
  J<- data.frame(J) %>% mutate(Rep=row.names(J)) %>% 
    mutate(Rep=paste0(Rep,".")) %>% 
    mutate(Rep=gsub(Rep,pattern=" ",replacement=".")) %>% 
    mutate(Rep2=substr(Rep,start = str_locate(Rep,pattern = "Rep_"),stop = nchar(Rep)),
           Analyte=substr(Rep,start = str_locate(Rep,pattern = "Analyte_"),stop = nchar(Rep)),
           IsotopeLabelType=substr(Rep,start = str_locate(Rep,pattern = "IsotopeLabelType_"),stop = nchar(Rep))) %>%
    mutate(Rep2=substr(Rep2,stop = str_locate(Rep2,pattern = "\\.")-1,start = 5),
           Analyte=substr(Analyte,stop = str_locate(Analyte,pattern = "\\.")-1,start = 9),
           IsotopeLabelType=substr(IsotopeLabelType,stop = str_locate(IsotopeLabelType,pattern = "\\.")-1,start = 18)) %>%
    rename(ID_Rep=Rep2, ID_Analyte=Analyte) %>% select(-Rep)
  
  
  
  
  Report<-H %>% left_join(I,by = c("ID_Rep","ID_Analyte","IsotopeLabelType")) %>%
    left_join(J,by = c("ID_Rep","ID_Analyte","IsotopeLabelType")) %>%
    left_join(K,by = c("ID_Rep","ID_Analyte","IsotopeLabelType")) %>%
    select(ID_Analyte,IsotopeLabelType,ID_Rep,Similarity.Score,MPRA.Score,Library.dotp,Intensity.Score,Score.MassError) %>% mutate(Comment=Comment)
  
  
  return(Report)
}

#' Report.Transition
#'
#' This function is Report.Transition
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Report.Transition()
Report.Transition<-function(Chrom.AnalyteX,Norm.ChromX,y){
  #Similarity Scores  
  SimS<-Similarity.Score.Report_informsMPRA(Norm.ChromX,y)
  H<-SimS$Report
  H<-data.frame(H)%>%
    mutate(ID_FragmentIon_charge=row.names(H)) %>%
    gather(key = "Rep",value = "Similarity.Score",1:(dim(H)[2]))
  H<- data.frame(H)%>% mutate(Rep=paste0(Rep,".")) %>% 
    mutate(Rep=substr(Rep,start = str_locate(Rep,pattern = "Rep_"),stop = nchar(Rep)),
           Analyte=substr(Rep,start = str_locate(Rep,pattern = "Analyte_"),stop = nchar(Rep)),
           IsotopeLabelType=substr(Rep,start = str_locate(Rep,pattern = "IsotopeLabelType_"),stop = nchar(Rep))) %>%
    mutate(Rep=substr(Rep,stop = str_locate(Rep,pattern = "\\.")-1,start = 5),
           Analyte=substr(Analyte,stop = str_locate(Analyte,pattern = "\\.")-1,start = 9),
           IsotopeLabelType=substr(IsotopeLabelType,stop = str_locate(IsotopeLabelType,pattern = "\\.")-1,start = 18)) %>%
    rename(ID_Rep=Rep, ID_Analyte=Analyte) %>% mutate(ID_Rep=as.numeric(ID_Rep),ID_Analyte=as.numeric(ID_Analyte))
  Y<-Chrom.AnalyteX %>% select(ID_Rep,ID_Analyte,IsotopeLabelType) %>%
    distinct()
  H<-H %>%left_join(Y, by = c("ID_Rep","ID_Analyte","IsotopeLabelType"))
  return(H)
}

#' Report.Transition.All
#'
#' This function is Report.Transition.All
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Report.Transition.All()
Report.Transition.All<-function(Norm.ChromX,y){
  y2=t(matrix(ifelse(Trans.Vector %in% Trans.Vector3[which(y==1)],1,0)))
  colnames(y2)=Trans.Vector
  
  H<-Similarity.Score.Report.AllTransitions(Norm.Chrom,y2)
  H<-data.frame(H)%>%
    mutate(ID_FragmentIon_charge=row.names(H)) %>%
    gather(key = "Rep",value = "Similarity.Score",1:(dim(H)[2]))
  H<- data.frame(H)%>% mutate(Rep=paste0(Rep,".")) %>%
    mutate(Rep=substr(Rep,start = str_locate(Rep,pattern = "Rep_"),stop = nchar(Rep)),
           Analyte=substr(Rep,start = str_locate(Rep,pattern = "Analyte_"),stop = nchar(Rep)),
           IsotopeLabelType=substr(Rep,start = str_locate(Rep,pattern = "IsotopeLabelType_"),stop = nchar(Rep))) %>%
    mutate(Rep=substr(Rep,stop = str_locate(Rep,pattern = "\\.")-1,start = 5),
           Analyte=substr(Analyte,stop = str_locate(Analyte,pattern = "\\.")-1,start = 9),
           IsotopeLabelType=substr(IsotopeLabelType,stop = str_locate(IsotopeLabelType,pattern = "\\.")-1,start = 18)) %>%
    rename(ID_Rep=Rep, ID_Analyte=Analyte) %>% mutate(ID_Rep=as.numeric(ID_Rep),ID_Analyte=as.numeric(ID_Analyte))
  Y<-Chrom.Analyte %>% select(ID_Rep,ID_Analyte,IsotopeLabelType) %>%
    distinct()
  H<-H %>%left_join(Y, by = c("ID_Rep","ID_Analyte","IsotopeLabelType"))
  return(H)
} ## Does it do something?

### Peak Boundaries functions:

#' Chromatographic.MPRA
#'
#' This function is Chromatographic.MPRA
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Chromatographic.MPRA()
Chromatographic.MPRA<- function(list,Chrom.IsPotentialPeak){
  list<-Map(cbind,list,Chrom.IsPotentialPeak)
  lapply(list, function(m){
    #m1=m[,-c(dim(m)[2])]
    #Is.Potential=m[,c(dim(m)[2])]
    max.l=apply(m,1,function(x){
      max.l=if(x[length(x)]==0){rep(0,length(x))} else {c(x[-length(x)]/max(x[-length(x)]),x[length(x)])}})
    max.l=t(max.l)
    o=apply(max.l, 1, function(x) {
      if(x[length(x)]==0){0} else{ 
        dot.p(x[-length(x)], m[1,-c(length(x))])}
    })
    p=matrix(o)
    row.names(p)=names(o)
    colnames(p)="MPRA"
    p
  })
}

#' Chromatographic.DotP
#'
#' This function is Chromatographic.DotP
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Chromatographic.DotP()
Chromatographic.DotP<- function(list,Chrom.IsPotentialPeak){
  list<-Map(cbind,list,Chrom.IsPotentialPeak)
  lapply(list, function(m){
    #m1=m[,-c(dim(m)[2])]
    #Is.Potential=m[,c(dim(m)[2])]
    max.l=apply(m,1,function(x){
      max.l=if(x[length(x)]==0){rep(0,length(x))} else {c(x[-length(x)]/max(x[-length(x)]),x[length(x)])}})
    max.l=t(max.l)
    o=apply(max.l, 1, function(x) {
      if(x[length(x)]==0){0} else{ 
        dot.p(x[-length(x)], m[2,-c(length(x))])}
    })
    p=matrix(o)
    row.names(p)=names(o)
    colnames(p)="Dotp"
    p
  })
}

#' Chromatographic.IsPotentialPeak
#'
#' This function is Chromatographic.IsPotentialPeak
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Chromatographic.IsPotentialPeak()
Chromatographic.IsPotentialPeak<-function(list,MinimalNumberOfConsecutivePoints){
  p=lapply(list, function(L){
    namesL=row.names(L)
    L=L[-c(1,2),]
    L0=apply(L,2,function(x){
      x<-ifelse(x>0,1,0)})
    
    L1=apply(L,2,function(x){ #Counts how many consecutive points are have with non-zero intensity and then tells if the transitions have more than 5 points, i.e. 1 if ok 0 if not
      x<-ifelse(is.na(x),0,x)
      x<-ifelse(x>0,1,0)
      y<-rle(as.numeric(x))
      z=ifelse(max(y$lengths[which(y$values==1)])>=MinimalNumberOfConsecutivePoints,1,0)})
    M0=sum(L1) # Num of transitions with more than N consecutive points with non-zero intensity
    if(M0<=2){p=rep(0,(dim(L)[1]+2))
    p=matrix(p)
    row.names(p)=namesL
    colnames(p)="isPotentialPeak"
    return(p)}
    
    #QQ<-combn(colnames(L[,which(as.numeric(L1)>=1)]),2) ## Multiply all combinations of pairwise columns 
    QQ<-combn(colnames(L[,which(as.numeric(L1)>=1)]),3) ## Multiply all combinations of a number of columns  
    
    L2=apply(QQ,2,function(x){
      l=L0[,x]
      l=apply(l,1,prod)
    })
    
    colnames(L2)=apply(QQ,2,function(x){
      paste0(x,collapse = "_")})
    #L2
    
    L22<-apply(L2,2,function(o){ ## Find regions where the intensity of N transitions is above 0 for 5 consecutive points
      M0<-which(o==0)
      if(!length(M0)==0){
        M<-data.frame(PointsZero=M0,NextPointZero=c(M0[-c(1)],length(o))) %>% mutate(Peak=NextPointZero-PointsZero-1) %>% filter(Peak>=MinimalNumberOfConsecutivePoints) 
        M<-apply(M,1,function(x){paste(x[1]:x[2], collapse = "_")})
        M<-paste(M, collapse = "_")
        M<-strsplit(as.character(M), "_")
        M<-c(as.numeric(M[[1]]))
        M<-M[!M %in% M0]
        
        p<-rep(0,length(o))
        p[M]<-1} else{
          p<-rep(1,length(o))}
      
      p=matrix(p)
      #row.names(p)=row.names(m)
      colnames(p)="isPotentialPeak"
      p=c(1,1,p)})
    
    p=ifelse(apply(L22,1,sum)>=1,1,0)
    p=matrix(p)
    row.names(p)=namesL
    colnames(p)="isPotentialPeak"
    p
    return(p)})
  return(p)
}

#' Chromatographic.IntScore
#'
#' This function is Chromatographic.IntScore
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Chromatographic.IntScore()
Chromatographic.IntScore<- function(list){
  lapply(list, function(m){
    sum.m=apply(m,1,function(x){sum(as.numeric(x),na.rm = T)})
    sum.m=sum.m/max(sum.m[-c(1,2)])
    sum.m[c(1,2)]=0
    p=matrix(sum.m)
    row.names(p)=names(sum.m)
    colnames(p)="IntScore"
    p
  })
}

#' Chromatographic.IntProductScore
#'
#' This function is Chromatographic.IntProductScore
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Chromatographic.IntProductScore()
Chromatographic.IntProductScore<- function(list){
  lapply(list, function(m){
    m=m+1
    prod.m=apply(m,1,function(x){log(prod(as.numeric(x),na.rm = T),base = 2)})
    prod.m=prod.m/max(prod.m[-c(1,2)])
    prod.m[c(1,2)]=0
    p=matrix(prod.m)
    row.names(p)=names(prod.m)
    colnames(p)="IntProductScore"
    p
  })
}

#' smooth.Criminal
#'
#' This function is smooth.Criminal
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' smooth.Criminal()
smooth.Criminal<-function(w){
  w=as.numeric(w)
  u=w
  u[1]=median(c(w[1],w[1+1]),na.rm = T)
  u[length(w)]=median(c(w[length(w)],w[length(w)-1]),na.rm = T)
  for(i in 2:length(w)-1){
    u[i]=median(c(w[i-1],w[i],w[i+1]),na.rm = T)}
  u=as.numeric(u)
  return(u)}

#' New_Boundaries
#'
#' This function is New_Boundaries
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' New_Boundaries()
New_Boundaries<-function(L,RowNum_at_max_Skor,Intensity_at_max_Skor){
  v=L$Skor
  w=L$IntScore
  h=RowNum_at_max_Skor
  j=RowNum_at_max_Skor
  
  while(h>1 & w[h]/w[RowNum_at_max_Skor]>=0.05 & v[[h]]>0.25 & h<length(w)-1){
    h=h-1
  }
  while(j>1 & j<(length(w)-1) & w[[j]]/w[[RowNum_at_max_Skor]]>=0.05 & v[j]>0.25){
    j=j+1
  }
  if((j-h+1)==1) {#avoid any point that only has one integration point
    if(h==1){j=j+1} else {h=h-1} # prevent having h <1, this would produce and error when calling for the hth row of the dataframe.
  } 
  while((j-h+1)<8 & # we need at least 8 points
        h>=1 & #the left side integration boundary is at least the first chromatographic point
        h<=(length(w)-1) & # the left integration point is not the last chromatographic point
        j>1 & # the right side integration boundary is not the first chromatographic point
        j<=length(w)# the right side integration boundary is lower than the last chromatographic point
  ){
    if(w[[j]]>=w[[h]]){
      if(h==1){j=j+1} else {h=h-1} # prevent having h <1, this would produce and error when calling for the hth row of the dataframe.
      
    } else {
      if(j==length(w)){h=h-1} else {j=j+1} } # prevent having j> length(w), this would produce and error when calling a rownumber larger than the dimension of the dataframe.
  } ## At least 8 points chosing the ones with the lowest intensity
  
  New_Boundaries=data.frame(left=L[h,]$Times,right=L[j,]$Times,numPoints=(j-h+1))
  return(New_Boundaries)
}

#' Data.Loader_DB
#'
#' This function loads the data for a given analyte (defined by its index) from the SQlite database. It creates several lists containing each different information.
#' For each list, each element of the list corresponds to an MS run.
#' 
#' @param index int: analyte for which the data is extracted.
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Data.Loader_DB(index = 1)
Data.Loader_DB<-function(index,D){
  # Loading the data
  
  db <- dbConnect(SQLite(), dbname=paste0("DB_",Name_Tag,".sqlite"))
  
  Chrom.Analyte =  dbGetQuery(db,paste0("select * from MainTable where ID_Analyte = ", index))
  
  dbDisconnect(db)
  
  
  ### Boundaries
  Boundaries<-tapply(paste(Chrom.Analyte$ID_FragmentIon_charge,Chrom.Analyte$MinStartTime,Chrom.Analyte$MaxEndTime,Chrom.Analyte$InterpolatedTimes,sep=','), paste0("Rep_",Chrom.Analyte$ID_Rep," Analyte_",Chrom.Analyte$ID_Analyte," IsotopeLabelType_",Chrom.Analyte$IsotopeLabelType), function(x){
    m=strsplit(x, ',') %>% unlist() %>% gsub(pattern=' *', replacement='') %>% matrix(nrow=length(x), byrow=T)
    m=t(m)
    m=m[-1,1]
    n=ifelse(round(as.numeric(m),5)>=round(as.numeric(m[1]),5) & round(as.numeric(m),5)<=round(as.numeric(m[2]),5),1,0)
    m=cbind(m,n)
    m=m[-c(1,2),]
    
    # Find replicate where only a single data point was integrated.
    o= as.numeric(m[,dim(m)[2]])  
    if(sum(o)==1) {
      # ^ If there is only one data point integarted then another one is added to 
      # avoid the autamatic conversion of single-row dataframes into vectors.
      r = ifelse(which(o==1)==1, which(o==1)+1, which(o==1)-1)
      # ^ Thishandles the case where the single data point is the first row of the dataframe,
      # then it adds the second row. If not, it adds the previous row.
      m[r, dim(m)[2]] = 1
    }
    
    colnames(m)=c("Times","IntegrationZone")
    m
  })
  
  ##### Chromatograms
  
  Chrom_Full<-tapply(paste(Chrom.Analyte$ID_FragmentIon_charge,Chrom.Analyte$InterpolatedIntensities,sep=','), paste0("Rep_",Chrom.Analyte$ID_Rep," Analyte_",Chrom.Analyte$ID_Analyte," IsotopeLabelType_",Chrom.Analyte$IsotopeLabelType), function(x){
    m=strsplit(x, ',') %>% unlist() %>% gsub(pattern=' *', replacement='') %>% matrix(nrow=length(x), byrow=T)
    m=t(m)
    
    colnames(m) <- m[1,]
    m=m[-1,]
    m= apply(m,2,as.numeric) 
    
    row.names(m) <- paste('Point', 1:nrow(m), sep='.')
    m
  })
  
  Chrom_Full<-rapply( Chrom_Full, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  Chrom_Full<-rapply( Chrom_Full, f=function(x) ifelse((x)<0,0,x), how="replace" )
  
  Chrom<-Map(cbind, Chrom_Full,Boundaries)
  Chrom<-lapply(X = Chrom,FUN = function(W){
    W=W[(W[,which(colnames(W)=="IntegrationZone")]==1), ,drop=FALSE]
    W=W[,1:(dim(W)[2]-2), drop=FALSE]
  })
  
  ##### Normalized Chromatograms
  
  Norm.Chrom<-lapply(Chrom,FUN = function(m){
    p=row.names(m)
    m= apply(m,2,as.numeric)
    m= apply(m, 2, function(x2)x2/max(x2, na.rm=T))*100
    row.names(m) <- p
    m
  })
  
  ##### Mass.Errors
  
  MassErrors_Full<-tapply(paste(Chrom.Analyte$ID_FragmentIon_charge,Chrom.Analyte$InterpolatedMassErrors,sep=','), paste0("Rep_",Chrom.Analyte$ID_Rep," Analyte_",Chrom.Analyte$ID_Analyte," IsotopeLabelType_",Chrom.Analyte$IsotopeLabelType), function(x){
    m=strsplit(x, ',') %>% unlist() %>% gsub(pattern=' *', replacement='') %>% matrix(nrow=length(x), byrow=T)
    m=t(m)
    colnames(m) <- m[1,]
    m=m[-1,]
    m= apply(m,2,as.numeric) 
    row.names(m) <- paste('Point', 1:nrow(m), sep='.')
    m
  })
  
  MassErrors_Full<-rapply( MassErrors_Full, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  MassErrors_Full<-Map(cbind, MassErrors_Full,Boundaries)
  ### Do not change the order of these operations!!!
  MassErrors<-lapply(X = MassErrors_Full,FUN = function(W){
    W=W[(W[,which(colnames(W)=="IntegrationZone")]==1),]
    W=W[,1:(dim(W)[2]-2)]
  })
  
  MassErrors<-Map(rbind,MassErrors,Norm.Chrom)
  colnames(MassErrors[[1]])<-ifelse(duplicated(colnames(MassErrors[[1]])),paste0("W_",colnames(MassErrors[[1]])),colnames(MassErrors[[1]]))
  
  
  MassErrors<-lapply(MassErrors,FUN =function(L){
    a=dim(L)[1]/2
    b=dim(L)[1]
    meanMassErrors.IntegratedZone<-apply(L,2,function(x){
      y=abs(as.numeric(x[1:a]))
      z=as.numeric(x[(a+1):b])
      weighted.mean(y,z,na.rm = T)
    })
    SkorMassErrors.IntegratedZone<-ifelse(as.numeric(meanMassErrors.IntegratedZone)<=MassError_Tolerance,1,ifelse(as.numeric(meanMassErrors.IntegratedZone)<=MassError_CutOff,MassError_CutOff/(MassError_CutOff-MassError_Tolerance)+as.numeric(meanMassErrors.IntegratedZone)*(1/(MassError_Tolerance-MassError_CutOff)),0))
    t(data.frame(meanMassErrors.IntegratedZone,SkorMassErrors.IntegratedZone))
  })  
  
  
  ## Library Intensities
  SpctLib <- Chrom.Analyte %>% 
    select(Area,LibraryIntensity, ID_Rep,ID_Analyte,IsotopeLabelType,ID_FragmentIon_charge) %>% 
    filter(ID_Rep== unique(ID_Rep)[1]) %>% filter(IsotopeLabelType== unique(IsotopeLabelType)[1]) %>% 
    mutate(Rank=rank(-LibraryIntensity,ties.method= "min")) %>% distinct()
  Transition.Lib.Intensity=as.matrix(t(SpctLib[,"LibraryIntensity"]))
  colnames(Transition.Lib.Intensity)=SpctLib[,"ID_FragmentIon_charge"]
  row.names(Transition.Lib.Intensity)="Lib.Intensity"
  
  ## Rank transitions by  Spectral Library Intensity
  Transition.Rank=as.matrix(t(SpctLib[,"Rank"]))
  colnames(Transition.Rank)=SpctLib[,"ID_FragmentIon_charge"]
  row.names(Transition.Rank)="Rank"
  
  ## DIA Area  per transition
  
  Transition.Area<- Chrom.Analyte %>%
    select(ID_FragmentIon_charge,Area,ID_Rep,ID_Analyte,IsotopeLabelType) %>% 
    mutate(Area = as.numeric(Area)) %>%
    distinct %>%
    group_by(ID_Rep,ID_Analyte,IsotopeLabelType) %>%
    spread(key = ID_FragmentIon_charge,value = Area)  %>% ungroup### Normalized DIA.Areas
  
  Transition.Area<- split(as.data.frame(Transition.Area,stringsAsFactors = F), f =paste0("Analyte_",Transition.Area$ID_Analyte," IsotopeLabelType_",Transition.Area$IsotopeLabelType), drop=FALSE) 
  
  Transition.Area<- lapply(Transition.Area,FUN = function(L) {
    L=L %>%
      select(-ID_Analyte,-IsotopeLabelType) %>%
      ungroup
    
    row.names(L)=paste0("Rep_",L$ID_Rep)
    
    L=L %>%select(-ID_Rep)
    return(L)})
  
  MPRA.MeanArea<-lapply(data.matrix(Transition.Area), FUN = function(L) {
    P=data.frame(t(apply(L,2, mean, na.rm=T)),stringsAsFactors = F)
    names(P)=gsub(names(P),pattern = "X",replacement = "")
    row.names(P)=c("Mean.Area")
    return(P)})
  
  #### Replicate-wide data
  
  for(i in 1:length(names(Transition.Area))){
    Transition.Area[[i]]=rbind(MPRA.MeanArea[[i]],Transition.Lib.Intensity,Transition.Area[[i]])
  }
  
  
  Transition.Area= lapply(Transition.Area,FUN = function(L){ MAX.L=apply(L,1,max,na.rm =T)
  SUM.L=apply(L,1,sum,na.rm =T)
  L=cbind(L,Total.To.Max.Ratio=SUM.L/MAX.L)})
  
  
  ### Chromatogram_Score
  for(i in 1:length(names(Chrom_Full))){
    Chrom_Full[[i]]=rbind(MPRA.MeanArea[[1]],Transition.Lib.Intensity,Chrom_Full[[i]])
  }
  
  return(list(Chrom.Analyte=Chrom.Analyte,
              Boundaries=Boundaries,
              Chrom_Full=Chrom_Full,
              Chrom=Chrom,
              Norm.Chrom=Norm.Chrom,
              MassErrors_Full=MassErrors_Full,
              MassErrors=MassErrors,
              Transition.Area=Transition.Area,
              Transition.Rank=Transition.Rank))
}

### Tools

#' Run_Transition_Refinment_Tool
#'
#' This function is Run_Transition_Refinment_Tool
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Run_Transition_Refinment_Tool()
Run_Transition_Refinment_Tool<-function(Chrom.Analyte,Norm.Chrom,MassErrors,Transition.Area,Transition.Rank){
  
  # %%%%%%%%%%%%% Parameters
  
  num_trans=dim(Norm.Chrom[[1]])[2]
  
  if (num_trans<MinimalNumberOfTransitionsAfterOptimization) {
    if (KeepPeptidesWithLowerNumberOfTransitions==TRUE) {
      Report.Replicate.Values<-Report.Replicate_informedMPRA(Norm.ChromX = Norm.Chrom,Transition.Area=Transition.Area,MassErrorsX = MassErrors,y = rep(1,num_trans), Comment ="aa-NonOptimized_LowerNumOfMinTrans")
      
      Report.Transition.Values<-Report.Transition(Chrom.AnalyteX = Chrom.Analyte,Norm.ChromX = Norm.Chrom,y = rep(1,num_trans))
      
      return(list(Report.Transition.Values=Report.Transition.Values,
                  Report.Replicate.Values=Report.Replicate.Values))
    } else{return(0)}
  } else{
    
    num_rep=dim(Transition.Area[[1]])[1]-2
    sug=rep(1,num_trans)
    Trans.Vector=colnames(Norm.Chrom[[1]])
    
    #     Functions
    
    Trans.Classified<-Transition.Classifier(Norm.Chrom = Norm.Chrom,Transition.Area = Transition.Area,num_trans = num_trans,Transition.Rank = Transition.Rank)
    sug.matrix<-Sug.Matrix.FUN(Trans.Classified,Trans.Vector = Trans.Vector)
    
    
    
    Transitions.To.Remove= data.frame(
      Trans.Classified %>% 
        select(ID_FragmentIon_charge,Rank.The25quartile.Initial.SimilarityScore,Rank.The25quartile.Initial.Area,Rank) %>%
        arrange(-Rank.The25quartile.Initial.SimilarityScore,-Rank) %>%
        select(ID_FragmentIon_charge))[c(0:num.trans.To.Remove.Fun(num_trans)),]
    ### Remove transition if the median mass error is higher than the mass error cutoff
    Transitions_To_Remove_Mass_Error<-data.frame(rbindlist(lapply(MassErrors,FUN = function(L){data.frame(Transitions=colnames(L),SkorMassErrors_IntegratedZone=L[1,])})))
    
    Transitions_To_Remove_Mass_Error<-Transitions_To_Remove_Mass_Error %>% group_by(Transitions) %>%
      summarise(median_error=median(SkorMassErrors_IntegratedZone,na.rm = T)) %>%
      arrange(-median_error) %>%
      filter(median_error>=MassError_CutOff) %>%
      select(Transitions) %>%
      mutate(Transitions=as.character(Transitions))
    
    Transitions.To.Remove<-unique(c(Transitions.To.Remove,Transitions_To_Remove_Mass_Error$Transitions))
    #
    
    
    Transitions.To.Remove=sapply(Trans.Vector,function(x) ifelse(x %in% Transitions.To.Remove,0,1))
    
    
    
    
    Trans.Classified_filtered<-Trans.Classified[Trans.Classified$ID_FragmentIon_charge %in% names(Transitions.To.Remove[Transitions.To.Remove==1]),]
    Trans.Vector2= names(Transitions.To.Remove[Transitions.To.Remove==1])
    num_trans2=length(Trans.Vector2)
    
    
    #########  New dataset without the worst transitions
    
    Norm.Chrom2<-Transition.Remover(Norm.Chrom,Transitions.To.Remove)
    
    Transition.Area2<-Transition.Remover(Transition.Area,Transitions.To.Remove)
    Transition.Area2=lapply(Transition.Area2,FUN = function(L){ MAX.L=apply(L,1,max,na.rm =T);SUM.L=apply(L,1,sum,na.rm =T); L=cbind(L,Total.To.Max.Ratio=SUM.L/MAX.L)})
    
    Transition.Rank2=t(matrix(Transition.Rank[which(Transitions.To.Remove==1)]))
    colnames(Transition.Rank2)=colnames(Transition.Rank)[which(Transitions.To.Remove==1)]
    row.names(Transition.Rank2)="Rank"
    
    MassErrors2=Transition.Remover(MassErrors,Transitions.To.Remove)
    
    sug.matrix2<-Sug.Matrix.FUN(Trans.Classified_filtered,Trans.Vector = Trans.Vector2)
    
    
    GA_TransitionOpt<-ga(type ="binary",
                         fitness = GA.Fitness,
                         Norm.Chrom=Norm.Chrom2,
                         Transition.Area=Transition.Area2,
                         Transition.Rank=Transition.Rank2,
                         MassErrors=MassErrors2,
                         nBits = num_trans2, 
                         suggestions = sug.matrix2, 
                         maxiter=10,run=3, 
                         popSize=ifelse(dim(sug.matrix2)[1]<=20,20,dim(sug.matrix2)[1]),
                         pcrossover=0.9,
                         pmutation=0.1,
                         elitism=0.1, 
                         keepBest = T,
                         seed = 112358)
    
    
    GA_1stRunSolution=Trans.Vector2[which(GA_TransitionOpt@solution[1,]==1)]
    
    ## {{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}
    
    Trans.Classified2<-Transition.Classifier2(Norm.Chrom = Norm.Chrom2,Transition.Area = Transition.Area2,num_trans = num_trans2,y = GA_TransitionOpt@solution[1,],Transition.Rank = Transition.Rank)
    
    Transitions.To.Remove2= data.frame(
      Trans.Classified2 %>% 
        select(ID_FragmentIon_charge,Rank.The25quartile.Initial.SimilarityScore,Rank.The25quartile.Initial.Area,Rank) %>%
        arrange(-Rank.The25quartile.Initial.SimilarityScore,-Rank) %>%
        select(ID_FragmentIon_charge))[c(0:num.trans.To.Remove.Fun(num_trans2)),]
    
    Transitions.To.Remove2=sapply(Trans.Vector2,function(x) ifelse(x %in% Transitions.To.Remove2,0,1))
    
    Trans.Classified_filtered2<-Trans.Classified2[Trans.Classified2$ID_FragmentIon_charge %in% names(Transitions.To.Remove2[Transitions.To.Remove2==1]),]
    
    Trans.Vector3= names(Transitions.To.Remove2[Transitions.To.Remove2==1])
    num_trans3=length(Trans.Vector3)
    
    
    #########  New dataset without the worst transitions
    
    Norm.Chrom3<-Transition.Remover(Norm.Chrom2,Transitions.To.Remove2)
    
    Transition.Area3<-Transition.Remover(Transition.Area2,Transitions.To.Remove2)
    Transition.Area3=lapply(Transition.Area3,FUN = function(L){ MAX.L=apply(L,1,max,na.rm =T);SUM.L=apply(L,1,sum,na.rm =T); L=cbind(L,Total.To.Max.Ratio=SUM.L/MAX.L)})
    
    Transition.Rank3=t(matrix(Transition.Rank2[which(Transitions.To.Remove2==1)]))
    colnames(Transition.Rank3)=colnames(Transition.Rank2)[which(Transitions.To.Remove2==1)]
    row.names(Transition.Rank3)="Rank"
    
    MassErrors3=Transition.Remover(MassErrors2,Transitions.To.Remove2)
    
    sug.matrix3<-Sug.Matrix.FUN(Trans.Classified_filtered2,Trans.Vector = Trans.Vector3)
    
    sug.matrix3<-rbind(sug.matrix3,ifelse(Trans.Vector3 %in% GA_1stRunSolution,1,0))
    
    GA_TransitionOpt_2<-ga(type ="binary",
                           fitness = GA.Fitness,
                           Norm.Chrom=Norm.Chrom3,
                           Transition.Area=Transition.Area3,
                           Transition.Rank=Transition.Rank3,
                           MassErrors=MassErrors3,
                           nBits = num_trans3, 
                           suggestions = sug.matrix3, 
                           maxiter=10,run=5, 
                           popSize=ifelse(dim(sug.matrix2)[1]<=20,20,dim(sug.matrix2)[1]),
                           pcrossover=0.9,
                           pmutation=0.1,
                           elitism=0.1, 
                           keepBest = T,
                           seed = 112358)
    
    
    Report.Replicate.Values<-rbind(Report.Replicate_informedMPRA(Norm.ChromX = Norm.Chrom,Transition.Area=Transition.Area,MassErrorsX = MassErrors,y = rep(1,num_trans), Comment ="a-Initial"),
                                   Report.Replicate_informedMPRA(Norm.ChromX = Norm.Chrom2,Transition.Area=Transition.Area2,MassErrorsX = MassErrors2,y = rep(1,num_trans2), Comment ="b-After 1st similarity score filter"),
                                   Report.Replicate_informedMPRA(Norm.ChromX = Norm.Chrom2,Transition.Area=Transition.Area2,MassErrorsX = MassErrors2,y = GA_TransitionOpt@solution[1,], Comment ="c-After 1st GA optimization"),
                                   Report.Replicate_informedMPRA(Norm.ChromX = Norm.Chrom3,Transition.Area=Transition.Area3,MassErrorsX = MassErrors3,y = rep(1,num_trans3), Comment ="d-After 2nd similarity score filter"),
                                   Report.Replicate_informedMPRA(Norm.ChromX = Norm.Chrom3,Transition.Area=Transition.Area3,MassErrorsX = MassErrors3,y = GA_TransitionOpt_2@solution[1,],Comment = "e-After 2nd GA optimization"))
    
    Report.Transition.Values<-Report.Transition(Chrom.AnalyteX = Chrom.Analyte,Norm.ChromX = Norm.Chrom3,y = GA_TransitionOpt_2@solution[1,])
    
    return(list(Report.Transition.Values=Report.Transition.Values,
                Report.Replicate.Values=Report.Replicate.Values))}
}




#' Run_PeakBoundaries_tool
#'
#' This function is Run_PeakBoundaries_tool
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Run_PeakBoundaries_tool()
Run_PeakBoundaries_tool<-function(Chrom_Full,MassErrors_Full,TransitionRefinementSolution,Boundaries){
  ### TransitionRemover
  Chrom_Full_2<-lapply(Chrom_Full,function(L){
    L=L[,which(colnames(L) %in% TransitionRefinementSolution)]})
  MassErrors_Full_2<-lapply(MassErrors_Full,function(L){
    L=L[,which(colnames(L) %in% TransitionRefinementSolution)]})
  
  MassErrors_Full_2<-lapply(MassErrors_Full_2,FUN = function(L){
    Mean.massErrors<-apply(L,1,function(x){mean(abs(as.numeric(x)))})
    #N<-ifelse(as.numeric(M)<=5,1,ifelse(as.numeric(M)<=15,1-as.numeric(M)/15,0))
    #N<-ifelse(as.numeric(M)<=5,1,ifelse(as.numeric(M)<=15,3/2+as.numeric(M)*(1-3/2)/5,0))
    Skor.MassErrors<-ifelse(as.numeric(Mean.massErrors)<=MassError_Tolerance,1,ifelse(as.numeric(Mean.massErrors)<=MassError_CutOff,MassError_CutOff/(MassError_CutOff-MassError_Tolerance)+as.numeric(Mean.massErrors)*(1/(MassError_Tolerance-MassError_CutOff)),0))
    #N=N^3
    data.frame(Mean.massErrors,Skor.MassErrors)
    #L=mean(L)
  })
  
  for(i in 1:length(names(MassErrors_Full_2))){
    MassErrors_Full_2[[i]]=rbind(rep(0,2),rep(0,2),MassErrors_Full_2[[i]])
  }
  
  # 
  ## finding a better potential peak for QTof data
  if(NonZeroBaselineChromatogram==TRUE){
    Chrom_Full_NoiseLimit<-lapply(Chrom_Full,FUN = function(L){
      min_int_AtEachPoint<-apply(L,1,function(x){ifelse(is.na(min(x,na.rm = T)),0,min(x,na.rm = T))})
      min_int_AtEachPoint<-min_int_AtEachPoint[-c(1,2)]
      #sd_int_AtEachPoint<-apply(L,1,function(x){ifelse(is.na(sd(x,na.rm = T)),0,sd(x,na.rm = T))})
      NoiseLimit<-median(min_int_AtEachPoint)+2*sd(min_int_AtEachPoint)
      L[3:dim(L)[1],]=t(apply(L[3:dim(L)[1],],1,function(x){ifelse(x-NoiseLimit<0,0,x-NoiseLimit)}))
      L
    })
    NoiseLimits<-lapply(Chrom_Full,FUN = function(L){
      min_int_AtEachPoint<-apply(L,1,function(x){ifelse(is.na(min(x,na.rm = T)),0,min(x,na.rm = T))})
      min_int_AtEachPoint<-min_int_AtEachPoint[-c(1,2)]
      #sd_int_AtEachPoint<-apply(L,1,function(x){ifelse(is.na(sd(x,na.rm = T)),0,sd(x,na.rm = T))})
      NoiseLimit<-median(min_int_AtEachPoint)+2*sd(min_int_AtEachPoint)
    })
    Chrom_Full_2_NoiseLimit<-lapply(Chrom_Full_NoiseLimit,function(L){
      L=L[,which(colnames(L) %in% TransitionRefinementSolution)]})
    Chrom.IsPotentialPeak<-Chromatographic.IsPotentialPeak(Chrom_Full_2_NoiseLimit,MinimalNumberOfConsecutivePoints)
  } else {
    Chrom.IsPotentialPeak<-Chromatographic.IsPotentialPeak(Chrom_Full_2,MinimalNumberOfConsecutivePoints)}
  
  
  
  Chrom.DotP<-Chromatographic.DotP(Chrom_Full_2,Chrom.IsPotentialPeak)
  Chrom.MPRA<-Chromatographic.MPRA(Chrom_Full_2,Chrom.IsPotentialPeak)
  Chrom.Int<-Chromatographic.IntScore(Chrom_Full_2)
  Chrom.Int.Prod<-Chromatographic.IntProductScore(Chrom_Full_2)
  
  ChromatoScores<-Map(f = cbind,Chrom.DotP,Chrom.MPRA,Chrom.Int,Chrom.Int.Prod,MassErrors_Full_2,Chrom.IsPotentialPeak)
  ChromatoScores<- lapply(ChromatoScores,FUN = function(L){
    L=data.frame(L[-c(1,2),])%>%mutate(MassError_smooth=smooth.Criminal(Skor.MassErrors),
                                       Skor=Dotp^3*MPRA^3*isPotentialPeak,
                                       Skor_smooth=smooth.Criminal(Dotp^3*MPRA^3)*IntProductScore^3*isPotentialPeak,
                                       Skor_smooth_withMassErrors=smooth.Criminal(Dotp^3*MPRA^3)*IntProductScore^3*smooth.Criminal(Skor.MassErrors)*isPotentialPeak) %>%
      mutate(smoothed_Skor=smooth.Criminal(Skor_smooth_withMassErrors))})
  ChromatoScores<-Map(cbind,ChromatoScores,Boundaries)
  ChromatoScores<-lapply(ChromatoScores,function(L){
    
    # cat(L, '\n')
    # L=ChromatoScores[[L]]
    #max.Skor=max(L$Skor_smooth_withMassErrors,na.rm = T)
    max.Skor=as.numeric(data.frame(L) %>% select(Skor_smooth_withMassErrors,IntScore) %>%
                          top_n(n = 1,Skor_smooth_withMassErrors) %>%
                          top_n(n = 1,IntScore) %>%
                          select(Skor_smooth_withMassErrors) #%>%
                        #mutate(Skor_smooth_withMassErrors=as.numeric(as.character(Skor_smooth_withMassErrors)))
    )
    
    Time_at_max_Skor=as.numeric(as.character((data.frame(L) %>% filter(Skor_smooth_withMassErrors==max.Skor) %>% mutate(Times=as.character(Times)) %>% select(Times))))[1]
    
    if(is.na(max.Skor)) return(data.frame(left=NA,right=NA,numPoints=NA,max.Skor=NA,RT_maxSkor=NA))
    if(max.Skor==0) return(data.frame(left=NA,right=NA,numPoints=NA,max.Skor=NA,RT_maxSkor=NA))
    
    #RowNum_at_max_Skor=as.numeric(row.names(L[L$Skor_smooth_withMassErrors==max.Skor,]))
    L$Row.Num<-1:dim(L)[1]
    RowNum_at_max_Skor=as.numeric((data.frame(L) %>% filter(Skor_smooth_withMassErrors==max.Skor) %>% select(Row.Num)))[1]
    
    Intensity_at_max_Skor=as.numeric((data.frame(L) %>% filter(Skor_smooth_withMassErrors==max.Skor) %>% select(IntScore)))[1]
    
    O<-New_Boundaries(L,RowNum_at_max_Skor,Intensity_at_max_Skor)
    O<-data.frame(O,max.Skor=max.Skor,RT_maxSkor=Time_at_max_Skor)
  })
  
  New_PeakBoundaries<-do.call(rbind,ChromatoScores)
  New_PeakBoundaries<-data.frame(New_PeakBoundaries) %>% mutate(Rep=row.names(.)) %>% mutate(Rep=gsub(Rep,pattern=" ",replacement=".")) %>% mutate(Rep=paste0(Rep,".")) %>%
    mutate(Rep2=substr(Rep,start = str_locate(Rep,pattern = "Rep_"),stop = nchar(Rep)),
           Analyte=substr(Rep,start = str_locate(Rep,pattern = "Analyte_"),stop = nchar(Rep)),
           IsotopeLabelType=substr(Rep,start = str_locate(Rep,pattern = "IsotopeLabelType_"),stop = nchar(Rep))) %>%
    mutate(Rep2=substr(Rep2,stop = str_locate(Rep2,pattern = "\\.")-1,start = 5),
           Analyte=substr(Analyte,stop = str_locate(Analyte,pattern = "\\.")-1,start = 9),
           IsotopeLabelType=substr(IsotopeLabelType,stop = str_locate(IsotopeLabelType,pattern = "\\.")-1,start = 18)) %>%
    rename(ID_Rep=Rep2, ID_Analyte=Analyte) %>% select(-Rep)
  
  return(Results_PeakBoundariesTool=list(New_PeakBoundaries=New_PeakBoundaries))}

#' ChromatogramsSkorPlots
#'
#' This function is ChromatogramsSkorPlots
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' ChromatogramsSkorPlots()
ChromatogramsSkorPlots<-function(Chrom_Full,MassErrors_Full,TransitionRefinementSolution,Boundaries,Num_rep){
  # Chrom_Full=Chrom_Full[which(str_detect(names(Chrom_Full),pattern = paste("^","Rep_",Num_rep," ",sep = "")))]
  # MassErrors_Full=MassErrors_Full[which(str_detect(names(MassErrors_Full),pattern = paste("^","Rep_",Num_rep," ",sep = "")))]
  # Boundaries=Boundaries[which(str_detect(names(Boundaries),pattern = paste("^","Rep_",Num_rep," ",sep = "")))]
  #
  
  ### TransitionRemover
  Chrom_Full_2<-lapply(Chrom_Full,function(L){
    L=L[,which(colnames(L) %in% TransitionRefinementSolution)]})
  MassErrors_Full_2<-lapply(MassErrors_Full,function(L){
    L=L[,which(colnames(L) %in% TransitionRefinementSolution)]})
  
  MassErrors_Full_2<-lapply(MassErrors_Full_2,FUN = function(L){
    Mean.massErrors<-apply(L,1,function(x){mean(abs(as.numeric(x)))})
    #N<-ifelse(as.numeric(M)<=5,1,ifelse(as.numeric(M)<=15,1-as.numeric(M)/15,0))
    #N<-ifelse(as.numeric(M)<=5,1,ifelse(as.numeric(M)<=15,3/2+as.numeric(M)*(1-3/2)/5,0))
    Skor.MassErrors<-ifelse(as.numeric(Mean.massErrors)<=MassError_Tolerance,1,ifelse(as.numeric(Mean.massErrors)<=MassError_CutOff,MassError_CutOff/(MassError_CutOff-MassError_Tolerance)+as.numeric(Mean.massErrors)*(1/(MassError_Tolerance-MassError_CutOff)),0))
    #N=N^3
    data.frame(Mean.massErrors,Skor.MassErrors)
    #L=mean(L)
  })
  
  for(i in 1:length(names(MassErrors_Full_2))){
    MassErrors_Full_2[[i]]=rbind(rep(0,2),rep(0,2),MassErrors_Full_2[[i]])
  }
  
  
  ## finding a better potential peak for QTof data
  if(NonZeroBaselineChromatogram==TRUE){
    Chrom_Full_NoiseLimit<-lapply(Chrom_Full,FUN = function(L){
      min_int_AtEachPoint<-apply(L,1,function(x){ifelse(is.na(min(x,na.rm = T)),0,min(x,na.rm = T))})
      min_int_AtEachPoint<-min_int_AtEachPoint[-c(1,2)]
      #sd_int_AtEachPoint<-apply(L,1,function(x){ifelse(is.na(sd(x,na.rm = T)),0,sd(x,na.rm = T))})
      NoiseLimit<-median(min_int_AtEachPoint)+2*sd(min_int_AtEachPoint)
      L[3:dim(L)[1],]=t(apply(L[3:dim(L)[1],],1,function(x){ifelse(x-NoiseLimit<0,0,x-NoiseLimit)}))
      L
    })
    NoiseLimits<-lapply(Chrom_Full,FUN = function(L){
      min_int_AtEachPoint<-apply(L,1,function(x){ifelse(is.na(min(x,na.rm = T)),0,min(x,na.rm = T))})
      min_int_AtEachPoint<-min_int_AtEachPoint[-c(1,2)]
      #sd_int_AtEachPoint<-apply(L,1,function(x){ifelse(is.na(sd(x,na.rm = T)),0,sd(x,na.rm = T))})
      NoiseLimit<-median(min_int_AtEachPoint)+2*sd(min_int_AtEachPoint)
    })
    Chrom_Full_2_NoiseLimit<-lapply(Chrom_Full_NoiseLimit,function(L){
      L=L[,which(colnames(L) %in% TransitionRefinementSolution)]})
    Chrom.IsPotentialPeak<-Chromatographic.IsPotentialPeak(Chrom_Full_2_NoiseLimit,MinimalNumberOfConsecutivePoints)
  } else {
    Chrom.IsPotentialPeak<-Chromatographic.IsPotentialPeak(Chrom_Full_2,MinimalNumberOfConsecutivePoints)}
  
  Chrom.DotP<-Chromatographic.DotP(Chrom_Full_2,Chrom.IsPotentialPeak)
  Chrom.MPRA<-Chromatographic.MPRA(Chrom_Full_2,Chrom.IsPotentialPeak)
  Chrom.Int<-Chromatographic.IntScore(Chrom_Full_2)
  Chrom.Int.Prod<-Chromatographic.IntProductScore(Chrom_Full_2)
  
  ChromatoScores<-Map(f = cbind,Chrom.DotP,Chrom.MPRA,Chrom.Int,Chrom.Int.Prod,MassErrors_Full_2,Chrom.IsPotentialPeak)
  ChromatoScores<- lapply(ChromatoScores,FUN = function(L){
    L=data.frame(L[-c(1,2),])%>%mutate(MassError_smooth=smooth.Criminal(Skor.MassErrors),
                                       Skor=Dotp^3*MPRA^3*isPotentialPeak,
                                       Skor_smooth=smooth.Criminal(Dotp^3*MPRA^3)*IntProductScore^3*isPotentialPeak,
                                       Skor_smooth_withMassErrors=smooth.Criminal(Dotp^3*MPRA^3)*IntProductScore^3*smooth.Criminal(Skor.MassErrors)*isPotentialPeak) %>%
      mutate(smoothed_Skor=smooth.Criminal(Skor_smooth_withMassErrors))})
  ChromatoScores<-Map(cbind,ChromatoScores,Boundaries)
  
  
  ChromatoScores2<-lapply(ChromatoScores,function(L){
    #max.Skor=max(L$Skor_smooth_withMassErrors,na.rm = T)
    max.Skor=as.numeric(data.frame(L) %>% select(Skor_smooth_withMassErrors,IntScore) %>%
                          top_n(n = 1,Skor_smooth_withMassErrors) %>%
                          top_n(n = 1,IntScore) %>%
                          select(Skor_smooth_withMassErrors) #%>%
                        #mutate(Skor_smooth_withMassErrors=as.numeric(as.character(Skor_smooth_withMassErrors)))
    )
    
    Time_at_max_Skor=as.numeric(as.character((data.frame(L) %>% filter(Skor_smooth_withMassErrors==max.Skor) %>% mutate(Times=as.character(Times)) %>% select(Times))))[1]
    
    if(is.na(max.Skor)) return(data.frame(left=NA,right=NA,numPoints=NA,max.Skor=NA,RT_maxSkor=NA))
    if(max.Skor==0) return(data.frame(left=NA,right=NA,numPoints=NA,max.Skor=NA,RT_maxSkor=NA))
    
    #RowNum_at_max_Skor=as.numeric(row.names(L[L$Skor_smooth_withMassErrors==max.Skor,]))
    L$Row.Num<-1:dim(L)[1]
    RowNum_at_max_Skor=as.numeric((data.frame(L) %>% filter(Skor_smooth_withMassErrors==max.Skor) %>% select(Row.Num)))[1]
    Intensity_at_max_Skor=as.numeric((data.frame(L) %>% filter(Skor_smooth_withMassErrors==max.Skor) %>% select(IntScore)))[1]
    
    O<-New_Boundaries(L,RowNum_at_max_Skor,Intensity_at_max_Skor)
    O<-data.frame(O,max.Skor=max.Skor,RT_maxSkor=Time_at_max_Skor)
  })
  
  New_PeakBoundaries<-do.call(rbind,ChromatoScores2)
  New_PeakBoundaries<-data.frame(New_PeakBoundaries) %>% mutate(Rep=row.names(.)) %>% mutate(Rep=gsub(Rep,pattern=" ",replacement=".")) %>% mutate(Rep=paste0(Rep,".")) %>%
    mutate(Rep2=substr(Rep,start = str_locate(Rep,pattern = "Rep_"),stop = nchar(Rep)),
           Analyte=substr(Rep,start = str_locate(Rep,pattern = "Analyte_"),stop = nchar(Rep)),
           IsotopeLabelType=substr(Rep,start = str_locate(Rep,pattern = "IsotopeLabelType_"),stop = nchar(Rep))) %>%
    mutate(Rep2=substr(Rep2,stop = str_locate(Rep2,pattern = "\\.")-1,start = 5),
           Analyte=substr(Analyte,stop = str_locate(Analyte,pattern = "\\.")-1,start = 9),
           IsotopeLabelType=substr(IsotopeLabelType,stop = str_locate(IsotopeLabelType,pattern = "\\.")-1,start = 18)) %>%
    rename(ID_Rep=Rep2, ID_Analyte=Analyte) %>% select(-Rep)
  
  U=ChromatoScores[[Num_rep]]
  #names(ChromatoScores[Num_rep])
  #min(as.numeric(as.character(U$Times)))
  #max(as.numeric(as.character(U$Times)))
  
  
  #U$Point=1:dim(U)[1]
  U<-U %>% select(-Mean.massErrors) %>%
    select(Times,isPotentialPeak,Dotp,MPRA,IntScore,IntProductScore,Skor.MassErrors,MassError_smooth,Skor,Skor_smooth,Skor_smooth_withMassErrors,Times,IntegrationZone,smoothed_Skor) %>%
    gather(key = "score",value = "Value",2:dim(.)[2])
  
  V<-U %>% filter(!score=="Skor.MassErrors",!score=="Dotp",!score=="MPRA",!score=="MassError_smooth")
  C=cbind(Chrom_Full[[Num_rep]][-c(1,2),],Boundaries[[Num_rep]]) %>% select(-IntegrationZone) %>%
    gather(key="Transition","Intensity", 1:(dim(.)[2]-1))
  
  U<-U%>%filter(score%in%c("isPotentialPeak","Dotp","MPRA","IntScore","IntProductScore","MassError_smooth","Skor_smooth_withMassErrors","IntegrationZone"))%>%
    mutate(score=gsub(score,pattern="IntegrationZone",replacement="a.IntegrationZone"))%>%
    mutate(score=gsub(score,pattern="isPotentialPeak",replacement="b.isPotentialPeak"))%>%
    mutate(score=gsub(score,pattern="Dotp",replacement="c.Dotp"))%>%
    mutate(score=gsub(score,pattern="MPRA",replacement="d.MPRA"))%>%
    mutate(score=gsub(score,pattern="IntScore",replacement="e.IntScore"))%>%
    mutate(score=gsub(score,pattern="IntProductScore",replacement="f.IntProductScore"))%>%
    mutate(score=gsub(score,pattern="MassError_smooth",replacement="g.MassError_smooth"))%>%
    mutate(score=gsub(score,pattern="Skor_smooth_withMassErrors",replacement="h.Skor_smooth_withMassErrors"))
  
  
  
  library(ggplot2)
  library(gridExtra)
  
  Graph1=ggplot(U,aes(x=as.numeric(as.character(Times)),y=as.numeric(Value),color=score))+geom_line()+facet_grid(score~.)+
    theme_bw()+
    scale_y_continuous(breaks = seq(0,1,0.5))+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$left)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$right)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$RT_maxSkor)),linetype=2,color="red")+theme(legend.position="top")
  
  Graph2=ggplot(U,aes(x=as.numeric(as.character(Times)),y=as.numeric(Value),color=score))+geom_line()+facet_grid(score~.)+theme_classic()+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$left)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$right)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$RT_maxSkor)),linetype=2,color="red")+theme(legend.position="top")+
    scale_x_continuous(limits = c(as.numeric(as.character(ChromatoScores2[[Num_rep]]$left)),as.numeric(as.character(ChromatoScores2[[Num_rep]]$right))))
  
  b=ggplot(C[which( C$Transition %in% TransitionRefinementSolution),],aes(x=as.numeric(as.character(Times)),y=as.numeric(Intensity),color=Transition))+
    geom_line()+theme_classic()+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$left)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$right)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$RT_maxSkor)),linetype=2,color="red")+theme(legend.position="top")+
    geom_hline(yintercept = 0,linetype=1,color="grey")+
    geom_hline(yintercept = NoiseLimits[[Num_rep]],linetype=2,color="black")
  
  
  # c=ggplot(V[V$score=="Skor_smooth_withMassErrors",],aes(x=as.numeric(as.character(Times)),y=as.numeric(Value)))+
  #   geom_line()+theme_classic()+
  #   geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$left)),linetype=2)+
  #   geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$right)),linetype=2)+
  #   geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$RT_maxSkor)),linetype=2,color="red")+theme(legend.position="top")
  
  PP=V[V$score=="Skor_smooth_withMassErrors",]
  #PP=V[V$score=="IntProductScore",]
  PP_median<-median(as.numeric(PP$Value),na.rm = T)
  PP_sd<-sd(as.numeric(PP$Value),na.rm = T)
  
  c=ggplot(PP,aes(x=as.numeric(as.character(Times)),y=as.numeric(Value)))+
    geom_line()+theme_classic()+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$left)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$right)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$RT_maxSkor)),linetype=2,color="red")+
    theme(legend.position="top")+
    geom_hline(yintercept = PP_median,linetype=2)+
    geom_hline(yintercept = PP_median+3*PP_sd,linetype=2,color="red")
  Graph3=grid.arrange(c,b)
  
  
  
  
  b=ggplot(C[which( C$Transition %in% TransitionRefinementSolution),],aes(x=as.numeric(as.character(Times)),y=as.numeric(Intensity),color=Transition))+
    geom_line()+theme_classic()+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$left)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$right)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$RT_maxSkor)),linetype=2,color="red")+theme(legend.position="top")+
    scale_x_continuous(limits = c(as.numeric(as.character(ChromatoScores2[[Num_rep]]$left)),as.numeric(as.character(ChromatoScores2[[Num_rep]]$right))))
  c=ggplot(V[V$score=="Skor_smooth_withMassErrors",],aes(x=as.numeric(as.character(Times)),y=as.numeric(Value)))+
    geom_line()+theme_classic()+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$left)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$right)),linetype=2)+
    geom_vline(xintercept = as.numeric(as.character(ChromatoScores2[[Num_rep]]$RT_maxSkor)),linetype=2,color="red")+theme(legend.position="top")+
    scale_x_continuous(limits = c(as.numeric(as.character(ChromatoScores2[[Num_rep]]$left)),as.numeric(as.character(ChromatoScores2[[Num_rep]]$right))))
  
  Graph4=grid.arrange(c,b)
  
  return(list(Graph1=Graph1,Graph2=Graph2,Graph3=Graph3,Graph4=Graph4))}

#' Run_Rescoring_Tool
#'
#' This function is Run_Rescoring_Tool
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Run_Rescoring_Tool()
Run_Rescoring_Tool<-function(Chrom.AnalyteX,TransitionRefinementSolution,Results_PeakBoundaries_tool,Comment){
  
  # Change integration Boundaries in Chrom_Analyte
  if(UseHeavyPeakBoundariesForLight==TRUE){
    Filtered_New_PeakBoundaries<-Results_PeakBoundaries_tool$New_PeakBoundaries %>%
      filter(IsotopeLabelType=="heavy") %>%
      select(-IsotopeLabelType)
    Chrom.Analyte_PB_Changed<-Chrom.AnalyteX %>% 
      mutate(ID_Rep=as.character(ID_Rep),ID_Analyte=as.character(ID_Analyte) )%>%
      left_join(Filtered_New_PeakBoundaries,by = c("ID_Rep", "ID_Analyte")) %>%
      mutate(MinStartTime=left,
             MaxEndTime=right) %>%
      select(-left,-right) %>%
      filter(ID_FragmentIon_charge %in% TransitionRefinementSolution)} 
  else {
    
    Chrom.Analyte_PB_Changed<-Chrom.AnalyteX %>% 
      mutate(ID_Rep=as.character(ID_Rep),ID_Analyte=as.character(ID_Analyte) )%>%
      left_join(Results_PeakBoundaries_tool$New_PeakBoundaries,by = c("ID_Rep", "ID_Analyte", "IsotopeLabelType")) %>%
      mutate(MinStartTime=left,
             MaxEndTime=right) %>%
      select(-left,-right) %>%
      filter(ID_FragmentIon_charge %in% TransitionRefinementSolution)}
  
  
  
  
  ### Boundaries
  Boundaries<-tapply(paste(Chrom.Analyte_PB_Changed$ID_FragmentIon_charge,Chrom.Analyte_PB_Changed$MinStartTime,Chrom.Analyte_PB_Changed$MaxEndTime,Chrom.Analyte_PB_Changed$InterpolatedTimes,sep=','), paste0("Rep_",Chrom.Analyte_PB_Changed$ID_Rep," Analyte_",Chrom.Analyte_PB_Changed$ID_Analyte," IsotopeLabelType_",Chrom.Analyte_PB_Changed$IsotopeLabelType), function(x){
    m=strsplit(x, ',') %>% unlist() %>% gsub(pattern=' *', replacement='') %>% matrix(nrow=length(x), byrow=T)
    m=t(m)
    m=m[-1,1]
    n=ifelse(round(as.numeric(m),5)>=round(as.numeric(m[1]),5) & round(as.numeric(m),5)<=round(as.numeric(m[2]),5),1,0)
    m=cbind(m,n)
    m=m[-c(1,2),]
    colnames(m)=c("Times","IntegrationZone")
    m
  })
  
  
  
  
  ##### Chromatograms
  
  Chrom_Full<-tapply(paste(Chrom.Analyte_PB_Changed$ID_FragmentIon_charge,Chrom.Analyte_PB_Changed$InterpolatedIntensities,sep=','), paste0("Rep_",Chrom.Analyte_PB_Changed$ID_Rep," Analyte_",Chrom.Analyte_PB_Changed$ID_Analyte," IsotopeLabelType_",Chrom.Analyte_PB_Changed$IsotopeLabelType), function(x){
    m=strsplit(x, ',') %>% unlist() %>% gsub(pattern=' *', replacement='') %>% matrix(nrow=length(x), byrow=T)
    m=t(m)
    #colnames(m) <- paste('Transition', 1:ncol(m), sep='.')
    
    colnames(m) <- m[1,]
    m=m[-1,]
    m= apply(m,2,as.numeric) 
    #m=apply(m, 2, function(x2)x2/max(x2, na.rm=T))*100
    row.names(m) <- paste('Point', 1:nrow(m), sep='.')
    m
  })
  
  Chrom_Full<-rapply( Chrom_Full, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  Chrom_Full<-rapply( Chrom_Full, f=function(x) ifelse((x)<0,0,x), how="replace" )
  Chrom<-Map(cbind, Chrom_Full,Boundaries)
  Chrom<-lapply(X = Chrom,FUN = function(W){
    W=W[(W[,which(colnames(W)=="IntegrationZone")]==1),,drop=F]
    W=W[,1:(dim(W)[2]-2)]
  })
  
  
  Areas_Recalc<-lapply(Chrom,FUN = function(L){
    if(is.null(dim(L)[1])){L} else
    {L=apply(X = L,MARGIN = 2,function(x){
      x=as.numeric(x)
      sum(x)})}
  })
  
  Areas_Recalc<-do.call(rbind,Areas_Recalc) 
  Areas_Recalc2 <-data.frame(Areas_Recalc) %>%
    mutate(Rep=row.names(.)) %>% mutate(Rep=gsub(Rep,pattern=" ",replacement=".")) %>% mutate(Rep=paste0(Rep,".")) %>% 
    mutate(Rep2=substr(Rep,start = str_locate(Rep,pattern = "Rep_"),stop = nchar(Rep)),
           Analyte=substr(Rep,start = str_locate(Rep,pattern = "Analyte_"),stop = nchar(Rep)),
           IsotopeLabelType=substr(Rep,start = str_locate(Rep,pattern = "IsotopeLabelType_"),stop = nchar(Rep))) %>%
    mutate(Rep2=substr(Rep2,stop = str_locate(Rep2,pattern = "\\.")-1,start = 5),
           Analyte=substr(Analyte,stop = str_locate(Analyte,pattern = "\\.")-1,start = 9),
           IsotopeLabelType=substr(IsotopeLabelType,stop = str_locate(IsotopeLabelType,pattern = "\\.")-1,start = 18)) %>%
    rename(ID_Rep=Rep2, ID_Analyte=Analyte) %>% select(-Rep) %>%
    select(ID_Rep,ID_Analyte,IsotopeLabelType,1:(dim(.)[2]-3)) %>%
    gather(key = "ID_FragmentIon_charge",value = "Area.Recalc",4:dim(.)[2]) %>%
    mutate(ID_FragmentIon_charge=gsub(ID_FragmentIon_charge,pattern="X",replacement=""))
  
  Chrom.Analyte_PB_Changed<-Chrom.Analyte_PB_Changed %>%
    mutate(ID_FragmentIon_charge=as.character(ID_FragmentIon_charge))%>%
    left_join(Areas_Recalc2, by = c("ID_Rep", "ID_FragmentIon_charge", "ID_Analyte", "IsotopeLabelType")) %>%
    mutate(Area=Area.Recalc) %>% select(-Area.Recalc)
  
  #Normalized Chromatogram
  Norm.Chrom<-lapply(Chrom,FUN = function(m){
    m= apply(m,2,as.numeric)
    m= apply(m, 2, function(x2)x2/max(x2, na.rm=T))*100
    row.names(m) <- paste('Point', 1:nrow(m), sep='.')
    m
  })
  
  ##### Mass.Errors
  
  MassErrors_Full<-tapply(paste(Chrom.Analyte_PB_Changed$ID_FragmentIon_charge,Chrom.Analyte_PB_Changed$InterpolatedMassErrors,sep=','), paste0("Rep_",Chrom.Analyte_PB_Changed$ID_Rep," Analyte_",Chrom.Analyte_PB_Changed$ID_Analyte," IsotopeLabelType_",Chrom.Analyte_PB_Changed$IsotopeLabelType), function(x){
    m=strsplit(x, ',') %>% unlist() %>% gsub(pattern=' *', replacement='') %>% matrix(nrow=length(x), byrow=T)
    m=t(m)
    colnames(m) <- m[1,]
    m=m[-1,]
    m= apply(m,2,as.numeric) 
    row.names(m) <- paste('Point', 1:nrow(m), sep='.')
    #m<-apply(m,2,mean)
    m
  })
  
  MassErrors_Full<-rapply( MassErrors_Full, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  MassErrors_Full<-Map(cbind, MassErrors_Full,Boundaries)
  ### Do not change the order of these operations!!!
  MassErrors<-lapply(X = MassErrors_Full,FUN = function(W){
    W=W[(W[,which(colnames(W)=="IntegrationZone")]==1),]
    W=W[,1:(dim(W)[2]-2)]
  })
  
  MassErrors<-Map(rbind,MassErrors,Norm.Chrom)
  colnames(MassErrors[[1]])<-ifelse(duplicated(colnames(MassErrors[[1]])),paste0("W_",colnames(MassErrors[[1]])),colnames(MassErrors[[1]]))
  
  
  MassErrors<-lapply(MassErrors,FUN =function(L){
    a=dim(L)[1]/2
    b=dim(L)[1]
    meanMassErrors.IntegratedZone<-apply(L,2,function(x){
      y=abs(as.numeric(x[1:a]))
      z=as.numeric(x[(a+1):b])
      weighted.mean(y,z,na.rm = T)
    })
    SkorMassErrors.IntegratedZone<-ifelse(as.numeric(meanMassErrors.IntegratedZone)<=MassError_Tolerance,1,ifelse(as.numeric(meanMassErrors.IntegratedZone)<=MassError_CutOff,MassError_CutOff/(MassError_CutOff-MassError_Tolerance)+as.numeric(meanMassErrors.IntegratedZone)*(1/(MassError_Tolerance-MassError_CutOff)),0))
    t(data.frame(meanMassErrors.IntegratedZone,SkorMassErrors.IntegratedZone))
    
    
  })  
  
  ## Library Intensities
  SpctLib <- Chrom.Analyte_PB_Changed %>% 
    select(Area,LibraryIntensity, ID_Rep,ID_Analyte,IsotopeLabelType,ID_FragmentIon_charge) %>% 
    filter(ID_Rep== unique(ID_Rep)[1]) %>% filter(IsotopeLabelType== unique(IsotopeLabelType)[1]) %>% 
    mutate(Rank=rank(-LibraryIntensity,ties.method= "min"))
  Transition.Lib.Intensity=as.matrix(t(SpctLib[,"LibraryIntensity"]))
  colnames(Transition.Lib.Intensity)=SpctLib[,"ID_FragmentIon_charge"]
  row.names(Transition.Lib.Intensity)="Lib.Intensity"
  #SummedTotalArea.Lib=apply(Transition.Lib.Intensity,1,sum)
  #Transition.Lib.Intensity=cbind(Transition.Lib.Intensity,SummedTotalArea=SummedTotalArea.Lib)
  
  ## Rank transitions by  Spectral Library Intensity
  Transition.Rank=as.matrix(t(SpctLib[,"Rank"]))
  colnames(Transition.Rank)=SpctLib[,"ID_FragmentIon_charge"]
  row.names(Transition.Rank)="Rank"
  
  
  ## DIA Area  per transition
  
  ### Change Transition Area!!!!
  Transition.Area<- Chrom.Analyte_PB_Changed %>% select(ID_FragmentIon_charge,Area,ID_Rep,ID_Analyte,IsotopeLabelType) %>% 
    distinct %>%
    group_by(ID_Rep,ID_Analyte,IsotopeLabelType) %>%
    spread(key = ID_FragmentIon_charge,value = Area)  %>% ungroup### Normalized DIA.Areas
  Transition.Area<- split(as.data.frame(Transition.Area,stringsAsFactors = F), f =paste0("Analyte_",Transition.Area$ID_Analyte," IsotopeLabelType_",Transition.Area$IsotopeLabelType), drop=FALSE) 
  Transition.Area<- lapply(Transition.Area,FUN = function(L) {L=L %>% select(-ID_Analyte,-IsotopeLabelType) %>% ungroup
  row.names(L)=paste0("Rep_",L$ID_Rep)
  L=L %>% select(-ID_Rep)
  return(L)})
  
  #MPRA.SummedTotalArea<-lapply(data.matrix(Transition.Area), FUN = function(L) apply(L,1, sum, na.rm=T))
  MPRA.MeanArea<-lapply(data.matrix(Transition.Area), FUN = function(L) {
    P=data.frame(t(apply(L,2, mean, na.rm=T)),stringsAsFactors = F)
    names(P)=gsub(names(P),pattern = "X",replacement = "")
    #SummedTotalArea=apply(P,1,sum)
    #P=cbind(P,SummedTotalArea=SummedTotalArea)
    row.names(P)=c("Mean.Area")
    return(P)})
  
  #### Replicate-wide data
  
  for(i in 1:length(names(Transition.Area))){
    Transition.Area[[i]]=rbind(MPRA.MeanArea[[i]],Transition.Lib.Intensity,Transition.Area[[i]])
  }
  
  
  Transition.Area= lapply(Transition.Area,FUN = function(L){ MAX.L=apply(L,1,max,na.rm =T)
  SUM.L=apply(L,1,sum,na.rm =T)
  L=cbind(L,Total.To.Max.Ratio=SUM.L/MAX.L)})
  
  
  num_trans=dim(Norm.Chrom[[1]])[2]
  num_rep=dim(Transition.Area[[1]])[1]-2
  sug=rep(1,num_trans)
  Trans.Vector=colnames(Norm.Chrom[[1]])
  
  
  Results_ReScore<-Report.Replicate_informedMPRA(Norm.ChromX = Norm.Chrom,
                                                 Transition.AreaX = Transition.Area,
                                                 MassErrorsX = MassErrors,
                                                 y = rep(1,num_trans),Comment = Comment)
  return(Results_ReScore)
}


### Wrappers
#' AvantGardeDIA_GlobalRefinement
#'
#' This function is AvantGardeDIA_GlobalRefinement
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' AvantGardeDIA_GlobalRefinement()
AvantGardeDIA_GlobalRefinement<-function(index,D){
  AG<-Data.Loader_DB(index,D = D)
  Chrom.Analyte=AG$Chrom.Analyte
  Boundaries=AG$Boundaries
  Chrom_Full=AG$Chrom_Full
  Chrom=AG$Chrom
  Norm.Chrom=AG$Norm.Chrom
  MassErrors_Full=AG$MassErrors_Full
  MassErrors=AG$MassErrors
  Transition.Area=AG$Transition.Area
  Transition.Rank=AG$Transition.Rank
  rm(AG)
  
  
  ########### Running the tools
  Results_TransitionRefinementTool<-Run_Transition_Refinment_Tool(Chrom.Analyte = Chrom.Analyte,
                                                                  Norm.Chrom = Norm.Chrom,
                                                                  MassErrors = MassErrors,
                                                                  Transition.Area = Transition.Area,
                                                                  Transition.Rank = Transition.Rank)
  TransitionRefinementSolution<-unique(Results_TransitionRefinementTool$Report.Transition.Values$ID_FragmentIon_charge)
  
  Results_PeakBoundaries_tool<-Run_PeakBoundaries_tool(Chrom_Full,MassErrors_Full,TransitionRefinementSolution,Boundaries)
  
  Results_ReScore<-Run_Rescoring_Tool(Chrom.AnalyteX = Chrom.Analyte,
                                      TransitionRefinementSolution =TransitionRefinementSolution,
                                      Results_PeakBoundaries_tool=Results_PeakBoundaries_tool,Comment = "f.GlobalRefinement")
  ColNames_Report=NULL
  ColNames_Report[1]=paste0(names(Results_TransitionRefinementTool$Report.Transition.Values),collapse = ";")
  ColNames_Report[2]=paste0(names(Results_TransitionRefinementTool$Report.Replicate.Values),collapse = ";")
  ColNames_Report[3]=paste0(names(Results_PeakBoundaries_tool$New_PeakBoundaries),collapse = ";")
  return(list(Results_TransitionRefinementTool=Results_TransitionRefinementTool,
              Results_PeakBoundaries_tool=Results_PeakBoundaries_tool,
              Results_ReScore=Results_ReScore,
              ColNames_Report=ColNames_Report))
  
}

#' AvantGardeDIA_TransitionRefinementAndReScore
#'
#' This function is AvantGardeDIA_TransitionRefinementAndReScore
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' AvantGardeDIA_TransitionRefinementAndReScore()
AvantGardeDIA_TransitionRefinementAndReScore<-function(index,D){
  
  AG<-Data.Loader_DB(index,D = D)
  
  Chrom.Analyte=AG$Chrom.Analyte
  Boundaries=AG$Boundaries
  Chrom_Full=AG$Chrom_Full
  Chrom=AG$Chrom
  Norm.Chrom=AG$Norm.Chrom
  MassErrors_Full=AG$MassErrors_Full
  MassErrors=AG$MassErrors
  Transition.Area=AG$Transition.Area
  Transition.Rank=AG$Transition.Rank
  
  ########### Running the tools
  
  Starting.PeakBoundaries<-list()
  Starting.PeakBoundaries$New_PeakBoundaries<-Chrom.Analyte %>% 
    select(left=MinStartTime,right=MaxEndTime,ID_Rep,ID_Analyte,IsotopeLabelType) %>% 
    distinct %>%
    mutate(numPoints=NA,max.Skor=NA,RT_maxSkor=NA) %>%
    select(left,right,numPoints,max.Skor,RT_maxSkor,ID_Rep,ID_Analyte,IsotopeLabelType) %>%
    mutate(ID_Rep=as.character(ID_Rep),
           ID_Analyte=as.character(ID_Analyte))
  
  Results_TransitionRefinementTool<-Run_Transition_Refinment_Tool(Chrom.Analyte = Chrom.Analyte,
                                                                  Norm.Chrom = Norm.Chrom,
                                                                  MassErrors = MassErrors,
                                                                  Transition.Area = Transition.Area,
                                                                  Transition.Rank = Transition.Rank)
  TransitionRefinementSolution<-unique(Results_TransitionRefinementTool$Report.Transition.Values$ID_FragmentIon_charge)
  
  Results_ReScore<-Run_Rescoring_Tool(Chrom.AnalyteX = Chrom.Analyte,
                                      TransitionRefinementSolution =TransitionRefinementSolution,
                                      Results_PeakBoundaries_tool = Starting.PeakBoundaries,
                                      Comment = "f.TransitionRefinementAndReScore")
  ColNames_Report=NULL
  ColNames_Report[1]=paste0(names(Results_TransitionRefinementTool$Report.Transition.Values),collapse = ";")
  ColNames_Report[2]=paste0(names(Results_TransitionRefinementTool$Report.Replicate.Values),collapse = ";")
  ColNames_Report[3]=NA
  return(list(Results_TransitionRefinementTool=Results_TransitionRefinementTool,
              Results_ReScore=Results_ReScore,
              ColNames_Report=ColNames_Report))
}

#' AvantGardeDIA_PeakBoundariesAndReScore
#'
#' This function is AvantGardeDIA_PeakBoundariesAndReScore
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' AvantGardeDIA_PeakBoundariesAndReScore()
AvantGardeDIA_PeakBoundariesAndReScore<-function(index,D){
  
  AG<-Data.Loader_DB(index,D = D)
  
  Chrom.Analyte=AG$Chrom.Analyte
  Boundaries=AG$Boundaries
  Chrom_Full=AG$Chrom_Full
  Chrom=AG$Chrom
  Norm.Chrom=AG$Norm.Chrom
  MassErrors_Full=AG$MassErrors_Full
  MassErrors=AG$MassErrors
  Transition.Area=AG$Transition.Area
  Transition.Rank=AG$Transition.Rank
  
  ########### Running the tools
  All.Trans<-unique(Chrom.Analyte$ID_FragmentIon_charge)
  
  Results_PeakBoundaries_tool<-Run_PeakBoundaries_tool(Chrom_Full,
                                                       MassErrors_Full,
                                                       TransitionRefinementSolution=All.Trans,
                                                       Boundaries)
  
  Results_ReScore<-Run_Rescoring_Tool(Chrom.AnalyteX = Chrom.Analyte,
                                      TransitionRefinementSolution =All.Trans,
                                      Results_PeakBoundaries_tool=Results_PeakBoundaries_tool,Comment = "f.PeakBoundariesAndReScore")
  ColNames_Report=NULL
  ColNames_Report[1]=NA
  ColNames_Report[2]=paste0(names(Results_ReScore),collapse = ";")
  ColNames_Report[3]=paste0(names(Results_PeakBoundaries_tool$New_PeakBoundaries),collapse = ";")
  return(list(Results_PeakBoundaries_tool=Results_PeakBoundaries_tool,
              Results_ReScore=Results_ReScore,
              ColNames_Report=ColNames_Report))
  
}

#' AvantGardeDIA_ReScore
#'
#' This function is AvantGardeDIA_ReScore
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' AvantGardeDIA_ReScore()
AvantGardeDIA_ReScore<-function(index,D){
  AG<-Data.Loader_DB(index,D = D)
  
  Chrom.Analyte=AG$Chrom.Analyte
  Boundaries=AG$Boundaries
  Chrom_Full=AG$Chrom_Full
  Chrom=AG$Chrom
  Norm.Chrom=AG$Norm.Chrom
  MassErrors_Full=AG$MassErrors_Full
  MassErrors=AG$MassErrors
  Transition.Area=AG$Transition.Area
  Transition.Rank=AG$Transition.Rank
  
  ########### Running the tools
  All.Trans<-unique(Chrom.Analyte$ID_FragmentIon_charge)
  Starting.PeakBoundaries<-list()
  Starting.PeakBoundaries$New_PeakBoundaries<-Chrom.Analyte %>% 
    select(left=MinStartTime,right=MaxEndTime,ID_Rep,ID_Analyte,IsotopeLabelType) %>% 
    distinct %>%
    mutate(numPoints=NA,max.Skor=NA,RT_maxSkor=NA) %>%
    select(left,right,numPoints,max.Skor,RT_maxSkor,ID_Rep,ID_Analyte,IsotopeLabelType) %>%
    mutate(ID_Rep=as.character(ID_Rep),
           ID_Analyte=as.character(ID_Analyte))
  
  ScoreNonOptimized<-Run_Rescoring_Tool(Chrom.Analyte = Chrom.Analyte,
                                        TransitionRefinementSolution =All.Trans,
                                        Results_PeakBoundaries_tool=Starting.PeakBoundaries,Comment = "NonOptimizedData")
  ColNames_Report=NULL
  ColNames_Report[1]=1
  ColNames_Report[2]=paste0(names(ScoreNonOptimized),collapse = ";")
  ColNames_Report[3]=3
  return(list(ScoreNonOptimized=ScoreNonOptimized,
              ColNames_Report=ColNames_Report))
  
}


#' AvantGardeDIA_SkorPlots
#'
#' This function is AvantGardeDIA_SkorPlots
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' AvantGardeDIA_SkorPlots()
AvantGardeDIA_SkorPlots<-function(index,D){
  AG<-Data.Loader_DB(index,D = D)
  Chrom.Analyte=AG$Chrom.Analyte
  Boundaries=AG$Boundaries
  Chrom_Full=AG$Chrom_Full
  Chrom=AG$Chrom
  Norm.Chrom=AG$Norm.Chrom
  MassErrors_Full=AG$MassErrors_Full
  MassErrors=AG$MassErrors
  Transition.Area=AG$Transition.Area
  Transition.Rank=AG$Transition.Rank
  rm(AG)
  
  
  ########### Running the tools
  Results_TransitionRefinementTool<-Run_Transition_Refinment_Tool(Chrom.Analyte = Chrom.Analyte,
                                                                  Norm.Chrom = Norm.Chrom,
                                                                  MassErrors = MassErrors,
                                                                  Transition.Area = Transition.Area,
                                                                  Transition.Rank = Transition.Rank)
  TransitionRefinementSolution<-unique(Results_TransitionRefinementTool$Report.Transition.Values$ID_FragmentIon_charge)
  
  #TransitionRefinementSolution<-unique(Chrom.Analyte$ID_FragmentIon_charge)
  A<-ChromatogramsSkorPlots(Chrom_Full,MassErrors_Full,TransitionRefinementSolution,Boundaries,6)
  p1=grid.arrange(A$Graph1,A$Graph3,ncol=2)
  B<-ChromatogramsSkorPlots(Chrom_Full,MassErrors_Full,TransitionRefinementSolution,Boundaries,1)
  p2=grid.arrange(B$Graph1,B$Graph3,ncol=2)
  grid.arrange(p1,p2,ncol=1)
}

#### Launchers
## LaunchAvantGardeInParallel And reading in chunks 
## Formatting

#' read_skyline_report_chunk
#'
#' This function is read_skyline_report_chunk
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' read_skyline_report_chunk()
read_skyline_report_chunk<-function(csv_file_path,n_rows,skip_rows,column_names){
  csv_chunk <- read_csv(file =csv_file_path,
                        n_max = n_rows,
                        skip = skip_rows,
                        col_names = column_names,
                        na= c("#N/A","NA","NaN"),
                        col_types = cols(
                          .default = col_character(),
                          IsDecoy = col_logical(),
                          PrecursorCharge = col_double(),
                          PrecursorMz = col_double(),
                          ProductCharge = col_double(),
                          ProductMz = col_double(),
                          Quantitative = col_logical(),
                          LibraryDotProduct = col_double(),
                          MinStartTime = col_double(),
                          MaxEndTime = col_double(),
                          Area = col_double(),
                          LibraryIntensity = col_double()))
  
  return(csv_chunk)}

#' AvantGardeDIA_InChunks_DB
#'
#' This function is AvantGardeDIA_InChunks_DB
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' AvantGardeDIA_InChunks_DB()
AvantGardeDIA_InChunks_DB<-function(D.file.name,RefinementWorkflow){
  
  AvG_writeParamsUsed("DB")

  # ColumnNames<-fread(file = D.file.name,header = T,stringsAsFactors = FALSE, na.strings = c("#N/A","NA","NaN"),skip=0,nrows = 1, sep = ",")
  
  ColumnNames <- read_csv(file =D.file.name,
                col_names = T,
                na= c("#N/A","NA","NaN"),
                n_max = 10)
  
  ColumnNames<-gsub(names(ColumnNames),pattern = "\\.",replacement = "")
  ColumnNames<-gsub(ColumnNames,pattern = " ",replacement = "")
  
  if(file.exists(paste0("DB_",Name_Tag,".sqlite"))){print(paste0("The SQLite database called: '",paste0("DB_",Name_Tag,".sqlite"),"' already exists in the working directory. Erase it or change the Name_Tag in the parameters file."))} else{
    # db <- dbConnect(SQLite(), dbname=paste0("DB_",Name_Tag,".sqlite"))
    
    ID_Analyte_catalog<-list()
    CurrentSkip<-list()
    u=1
    Stop.Condition=FALSE
    ID_Analyte_Survivors<-list()
    
    ### Read in chunks
    while(Stop.Condition==FALSE){
      if(u==1){
        #db <- dbConnect(SQLite(), dbname=paste0("DB_",Name_Tag,".sqlite"))
        
        D<-ReadFileInChunks_Format_And_Filter(D.file.name,ColumnNames = ColumnNames,
                                              PreviousMetaData = NULL,
                                              Skip = 1,
                                              is.FirstIteration = T)
        
        write.table(x = names(D$Data),file = paste0(MultiFile.path,"/Chrom.Analyte.Names_",Name_Tag,".csv"),quote = F, row.names = F, col.names = T,                  sep = ";")
        #dbDisconnect(db)
        
      } else{
        #db <- dbConnect(SQLite(), dbname=paste0("DB_",Name_Tag,".sqlite"))
        
        D<-ReadFileInChunks_Format_And_Filter(D.file.name,ColumnNames = ColumnNames,
                                              PreviousMetaData = CurrentMetaData,
                                              Skip =CurrentSkip[[u-1]],
                                              is.FirstIteration = F)
        #dbDisconnect(db)
        
      }
      
      Stop.Condition=!D$ChunkLastRowNumber==ReadNumberOfLines
      CurrentMetaData=D$MetaData
      CurrentSkip[[u]]=D$Skip
      ID_Analyte_catalog[[u]]=data.frame(u,
                                         Start_ID_Analyte=unique(D$Data$ID_Analyte)[1],
                                         End_ID_Analyte=unique(D$Data$ID_Analyte)[length(unique(D$Data$ID_Analyte))])
      ID_Analyte_Survivors[[u]]<-unique(D$Data$ID_Analyte)
      db <- dbConnect(SQLite(), dbname=paste0("DB_",Name_Tag,".sqlite"))
      
      dbWriteTable(conn = db, name = "MainTable", value = D$Data, row.names = FALSE,append= TRUE)
      dbDisconnect(db)
      
      
      rm(D)
      print(u)
      u=u+1
    }
    
    ID_Analyte_Survivors<-unique(do.call(what = c,ID_Analyte_Survivors))
    CurrentMetaData$M_Analyte$isSurvivor<-ifelse(CurrentMetaData$M_Analyte$ID_Analyte %in% ID_Analyte_Survivors,1,0)
    ID_Analyte_catalog=data.frame(do.call(rbind,ID_Analyte_catalog))
    ID_Analyte_catalog[dim(ID_Analyte_catalog)[1],dim(ID_Analyte_catalog)[2]]<-max(CurrentMetaData$M_Analyte$ID_Analyte)
    #write.table(x = ID_Analyte_catalog,file = paste0(MultiFile.path,"/ID_Analyte_catalog_",Name_Tag,".csv"),quote = F, row.names = F, col.names = T,sep = ";")
    #saveRDS(CurrentMetaData,file = paste0(MultiFile.path,"/MetaData_",Name_Tag,".RDS"))
    
    db <- dbConnect(SQLite(), dbname=paste0("DB_",Name_Tag,".sqlite"))
    dbWriteTable(conn = db, name = "MetaData_Analyte", value = CurrentMetaData$M_Analyte, row.names = FALSE,append= TRUE)
    dbWriteTable(conn = db, name = "MetaData_Replicate", value = CurrentMetaData$M_Replicate, row.names = FALSE,append= TRUE)
    dbWriteTable(conn = db, name = "MetaData_TransitionsType", value = CurrentMetaData$M_TransitionsType, row.names = FALSE,append= TRUE)
    dbWriteTable(conn = db, name = "MetaData_Transitions", value = CurrentMetaData$M_Transitions, row.names = FALSE,append= TRUE)
    dbWriteTable(conn = db, name = "MetaData_PrecursorResults", value = CurrentMetaData$M_PrecursorResult, row.names = FALSE,append= TRUE)
    
    print("Indexing Database...")
    dbGetQuery(conn = db,"CREATE INDEX index_Analyte ON MainTable (ID_Analyte)")
    
    dbDisconnect(db)
    print("Database: Done!")
  }
} # Writes in SQLite DB

#' LaunchParallelTasks_ReadFromDB
#'
#' This function is LaunchParallelTasks_ReadFromDB
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' LaunchParallelTasks_ReadFromDB()
LaunchParallelTasks_ReadFromDB<-function(ParamsFile,RefinementWorkflow){
  library(foreach)
  library(doSNOW)
  library(parallel)
  library(RSQLite)
  
  ParamsFile=ParamsFile
  source(ParamsFile)
  AvG_writeParamsUsed("Optimization")
  
  db <- dbConnect(SQLite(), dbname=paste0("DB_",Name_Tag,".sqlite"))
  X<-dbSendQuery(db, "SELECT ID_analyte FROM MetaData_Analyte WHERE isSurvivor=1 ") 
  X<- (dbFetch(X)%>% select(ID_Analyte) %>% distinct())[,1]
  dbDisconnect(db)
  
  cl<-makeCluster(parallel::detectCores() -1)
  registerDoSNOW(cl)
  I<-foreach(i=X,
             .errorhandling = "remove",
             .export =c('Indexing','Filter_Na_Shared_Or_LowMassTransitions',
                        'MetaDataLibraryConcatenator','ReadFileInChunks_Format_And_Filter','read_skyline_report_chunk',
                        'dot.p','LinearEquationByRange','ScoreTransform.dotP','Rank.filter',
                        'FindTranswithConsecutivePoints','Similarity.Score_informsMPRA','Similarity.Score.Report_informsMPRA',
                        'Similarity.Score.Report.AllTransitions','MPRA.Score.dotp','Library.dotp',
                        'MPRA_y_SpectLib.Fitness','MPRA_y_SpectLib.Report','MPRA_y_SpectLib.Fitness_informed',
                        'MPRA_y_SpectLib.Report_informed','Intensity.Fitness','Intensity.Fitness.Report',
                        'MassError.Score.Fitness_informsMPRA','MassError.Score.Report','Transition.Remover',
                        'Transition.Classifier','Transition.Classifier2','Sug.Matrix.FUN','GA.Fitness',
                        'num.trans.To.Remove.Fun','Report.Replicate_informedMPRA','Report.Transition',
                        'Report.Transition.All','Chromatographic.MPRA',
                        'Chromatographic.DotP',
                        'Chromatographic.IsPotentialPeak','Chromatographic.IntScore','Chromatographic.IntProductScore',
                        'smooth.Criminal',
                        'New_Boundaries','Data.Loader_DB','Run_Transition_Refinment_Tool',
                        'Run_PeakBoundaries_tool','ChromatogramsSkorPlots','Run_Rescoring_Tool',
                        'AvantGardeDIA_GlobalRefinement',
                        'AvantGardeDIA_TransitionRefinementAndReScore','AvantGardeDIA_PeakBoundariesAndReScore',
                        'AvantGardeDIA_ReScore','AvantGardeDIA_InChunks_DB',
                        'AvantGardeDIA_AnalyzeSingleAnalyte','Read_AndFormatResults')) %dopar% {
                          #source(LatestVersion)
                          source(ParamsFile)
                          
                          AvantGardeDIA_AnalyzeSingleAnalyte(i,RefinementWorkflow = RefinementWorkflow)
                          # if(Perform.Optimization==TRUE){
                          #   A<-AvantGardeDIA_GlobalRefinement(i,D)
                          #   
                          #   write.table(A$Results_TransitionRefinementTool$Report.Transition.Values,file=paste0(dir.output,"/Report.Transitions_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
                          #   write.table(A$Results_TransitionRefinementTool$Report.Replicate.Values,file=paste0(dir.output,"/Report.Replicate_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
                          #   write.table(A$Results_PeakBoundaries_tool$New_PeakBoundaries,file=paste0(dir.output,"/Report.PeakBoundaries_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
                          #   write.table(A$Results_ReScore,file=paste0(dir.output,"/Report.ReScore_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
                          #   
                          #   if(!file.exists(paste0(dir.output,"/Report.ColNames",Name_Tag,".csv"))){
                          #     write.table(A$ColNames_Report,file=paste0(dir.output,"/Report.ColNames",Name_Tag,".csv"),quote=F,row.names=F,col.names=F,sep=";")}
                          #   
                          #   rm(A)
                          #   return(i)} else {
                          #     A<-AvantGardeDIA_ReScore(i)
                          #     
                          #     write.table(A$ScoreNonOptimized,file=paste0(dir.output,"/ScoreNonOptimized_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
                          #     
                          #     if(!file.exists(paste0(dir.output,"/Report.ColNames",Name_Tag,".csv"))){
                          #       write.table(A$ColNames_Report,file=paste0(dir.output,"/Report.ColNames",Name_Tag,".csv"),quote=F,row.names=F,col.names=F,sep=";")}
                          #     
                          #     rm(A)
                          #     return(i)
                          #   }
                          dbDisconnect(db)
                        }
  dbDisconnect(db)
  
  stopCluster(cl)
  print("Parallel: Done!")
  
} ## For cluster Workflow

#' AvantGardeDIA_AnalyzeSingleAnalyte
#'
#' This function is AvantGardeDIA_AnalyzeSingleAnalyte
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' AvantGardeDIA_AnalyzeSingleAnalyte()
AvantGardeDIA_AnalyzeSingleAnalyte<-function(i,D,RefinementWorkflow){
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(data.table)
  require(sqldf)
  require(foreach)
  require(doSNOW)
  require(GA, quietly = T)
  require(parallel)
  
  if(RefinementWorkflow=="GlobalRefinement"){
    A<-AvantGardeDIA_GlobalRefinement(i,D)
    
    write.table(A$Results_TransitionRefinementTool$Report.Transition.Values,file=paste0(dir.output,"/Report_GR_Transitions_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
    write.table(A$Results_TransitionRefinementTool$Report.Replicate.Values,file=paste0(dir.output,"/Report_GR_Replicate_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
    write.table(A$Results_PeakBoundaries_tool$New_PeakBoundaries,file=paste0(dir.output,"/Report_GR_PeakBoundaries_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
    write.table(A$Results_ReScore,file=paste0(dir.output,"/Report_GR_ReScore_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
    
    if(!file.exists(paste0(dir.output,"/Report_GR_ColNames",Name_Tag,".csv"))){
      write.table(A$ColNames_Report,file=paste0(dir.output,"/Report_GR_ColNames",Name_Tag,".csv"),quote=F,row.names=F,col.names=F,sep=";")}
    
    rm(A)
    return(i)}
  if(RefinementWorkflow=="TransitionRefinement"){
    A<-AvantGardeDIA_TransitionRefinementAndReScore(i,D)
    
    write.table(A$Results_TransitionRefinementTool$Report.Transition.Values,file=paste0(dir.output,"/Report_TR_Transitions_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
    write.table(A$Results_TransitionRefinementTool$Report.Replicate.Values,file=paste0(dir.output,"/Report_TR_Replicate_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
    write.table(A$Results_ReScore,file=paste0(dir.output,"/Report_TR_ReScore_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
    
    if(!file.exists(paste0(dir.output,"/Report_TR_ColNames",Name_Tag,".csv"))){
      write.table(A$ColNames_Report,file=paste0(dir.output,"/Report_TR_ColNames",Name_Tag,".csv"),quote=F,row.names=F,col.names=F,sep=";")}
    
    rm(A)
    return(i)
  }
  if(RefinementWorkflow=="PeakBoundariesRefinement"){
    A<-AvantGardeDIA_PeakBoundariesAndReScore(i,D)
    
    write.table(A$Results_PeakBoundaries_tool$New_PeakBoundaries,file=paste0(dir.output,"/Report_PB_PeakBoundaries_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
    write.table(A$Results_ReScore,file=paste0(dir.output,"/Report_PB_ReScore_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
    
    
    if(!file.exists(paste0(dir.output,"/Report_PB_ColNames",Name_Tag,".csv"))){
      write.table(A$ColNames_Report,file=paste0(dir.output,"/Report_PB_ColNames",Name_Tag,".csv"),quote=F,row.names=F,col.names=F,sep=";")}
    
    rm(A)
    return(i)
  }
  if(RefinementWorkflow=="OnlyScoring"){
    A<-AvantGardeDIA_ReScore(i,D)
    
    write.table(A$ScoreNonOptimized,file=paste0(dir.output,"/Report_NO_ScoreNonOptimized_",Name_Tag,"_",i,".csv"),quote=F,row.names=F,col.names=F,sep=";")
    
    if(!file.exists(paste0(dir.output,"/Report_NO_ColNames",Name_Tag,".csv"))){
      write.table(A$ColNames_Report,file=paste0(dir.output,"/Report_NO_ColNames",Name_Tag,".csv"),quote=F,row.names=F,col.names=F,sep=";")}
    
    rm(A)
    return(i)
  }
  
} ## For UGER

## Read Results

#' Read_AndFormatResults
#'
#' This function reads, concatenates and formats the results.
#' @param Defaults
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' Read_AndFormatResults()
Read_AndFormatResults<-function(){
  Folder_3="ResultsOptimization"
  dir.create(file.path(getwd(),Folder_3),showWarnings = F)
  dir.output.results=file.path(getwd(),Folder_3)
  TempFiles.Location<-dir.output
  
  ##
  RefinementTAG<-if(RefinementWorkflow=="GlobalRefinement") {"GR"} else{
    if(RefinementWorkflow=="TransitionRefinement") {"TR"} else{
      if(RefinementWorkflow=="PeakBoundariesRefinement") {"PB"} else{
        if(RefinementWorkflow=="OnlyScoring") {"NO"}
      }}}
  
  ## Column Names
  ColumnNames<-read.csv(file = paste0(dir.output,"/Report_",RefinementTAG,"_ColNames",Name_Tag,".csv"),header = F,sep=' ', stringsAsFactors = FALSE)
  
  ## MetaData
  db <- dbConnect(SQLite(), dbname=paste0("DB_",Name_Tag,".sqlite"))
  MetaData_Analytes<-dbReadTable(db, "MetaData_Analyte") %>% mutate(ID_Analyte=as.character(ID_Analyte))
  MetaData_Replicate<-dbReadTable(db, "MetaData_Replicate") %>% mutate(ID_Rep=as.character(ID_Rep))
  MetaData_Transitions<-dbReadTable(db, "MetaData_Transitions") %>% mutate(ID_Analyte=as.character(ID_Analyte))
  MetaData_PrecursorResults<-dbReadTable(db, "MetaData_PrecursorResults")
  
  dbDisconnect(db)
  
  ## Peak_BOundaries
  PBlist<-list.files(dir.output,pattern = paste0("Report_",RefinementTAG,"_PeakBoundaries_"))
  ListFiles_PeakBoundaries<-paste0(dir.output,"/",PBlist)
  if(length(PBlist)>=1){
    l <- lapply(ListFiles_PeakBoundaries, fread, header = F,sep=';', stringsAsFactors = FALSE)
    NewPeakBoundaries_All_Results <- rbindlist( l )
    
    colnames(NewPeakBoundaries_All_Results)<-strsplit(ColumnNames[3,],split = ";",fixed = T)[[1]]
    
    NewPeakBoundaries_All_Results<-NewPeakBoundaries_All_Results %>% 
      mutate(left=ifelse(is.na(left),"#N/A",left),
             right=ifelse(is.na(right),"#N/A",right)) %>% 
      mutate(ID_Analyte=as.character(ID_Analyte),
             ID_Rep=as.character(ID_Rep))%>%
      left_join(MetaData_Analytes, by = c("ID_Analyte")) %>%
      left_join(MetaData_Replicate, by = c("ID_Rep"))
    
    write.table(NewPeakBoundaries_All_Results,file=paste0(dir.output.results,"/",RefinementTAG,"_NewPeakBoundaries_",Name_Tag,".csv"),quote=F,row.names=F,col.names=T,sep=",")
    
    Formatted<-NewPeakBoundaries_All_Results %>%
      select(FileName,PeptideModifiedSequence,left,right,PrecursorCharge,IsDecoy,IsotopeLabelType) %>%
      rename(PrecursorIsDecoy=IsDecoy) %>%
      rename(MinStartTime=left,MaxEndTime=right) %>%
      distinct()
    
    if(UseHeavyPeakBoundariesForLight==TRUE) {Formatted<-Formatted %>% filter(IsotopeLabelType=="heavy")} ##Only for P100
    
    write.csv(Formatted,file=paste0(dir.output.results,"/",RefinementTAG,"_NewPeakBoundaries_",Name_Tag,"_Formated.csv"),quote=F,row.names=F)
  }
  
  ## Report Transitions
  Translist<-list.files(dir.output,pattern = paste0("Report_",RefinementTAG,"_Transitions_"))
  ListFiles_Transitions<-paste0(dir.output,"/",Translist)
  if(length(Translist)>=1){
    
    l_trans <- lapply(ListFiles_Transitions, fread, header = F,sep=';', stringsAsFactors = FALSE)
    NewTransitions_All_Results <- rbindlist( l_trans )
    colnames(NewTransitions_All_Results)<-strsplit(ColumnNames[1,],split = ";",fixed = T)[[1]]
    
    NewTransitions_All_Results<-NewTransitions_All_Results %>% 
      mutate(ID_Analyte=as.character(ID_Analyte),
             ID_Rep=as.character(ID_Rep)) %>%
      left_join(MetaData_Transitions, by = c("ID_FragmentIon_charge","ID_Analyte", "IsotopeLabelType")) %>%
      left_join(MetaData_Replicate, by = c("ID_Rep"))
    
    write.table(NewTransitions_All_Results,file=paste0(dir.output.results,"/",RefinementTAG,"_NewTransitions_",Name_Tag,".csv"),quote=F,row.names=F,col.names=T,sep=",")
    
    
    NewTransitions_All_Results_OnlyTransitionsInfo<-data.frame(NewTransitions_All_Results %>% select(-Similarity.Score,-ID_Rep) %>%
                                                                 distinct())
    write.table(NewTransitions_All_Results_OnlyTransitionsInfo,file=paste0(dir.output.results,"/",RefinementTAG,"_NewTransitions_","OnlyTransitionInfo",Name_Tag,".csv"),quote=F,row.names=F,col.names=T,sep=",")
  }   
  
  
  ## Report Replicates
  RepList<-list.files(dir.output,pattern = paste0("Report_",RefinementTAG,"_Replicate_"))
  ListFiles_Replicates<-paste0(dir.output,"/",RepList)
  if(length(RepList)>=1){
    l_rep <- lapply(ListFiles_Replicates, fread, header = F,sep=';', stringsAsFactors = FALSE)
    NewReplicates_All_Results <- rbindlist( l_rep )
    colnames(NewReplicates_All_Results)<-strsplit(ColumnNames[2,],split = ";",fixed = T)[[1]]
    
    NewReplicates_All_Results<-NewReplicates_All_Results %>% 
      mutate(ID_Analyte=as.character(ID_Analyte),
             ID_Rep=as.character(ID_Rep)) %>%
      left_join(MetaData_Analytes, by = c("ID_Analyte")) %>%
      left_join(MetaData_Replicate, by = c("ID_Rep"))
    
    Exponents<-c(9.5,4.5,2.5,0.5)
    NewReplicates_All_Results<-NewReplicates_All_Results %>%
      mutate(Skor=Similarity.Score^Exponents[1]*Library.dotp^Exponents[2]*Score.MassError^Exponents[3]*MPRA.Score^Exponents[4])
    
    write.table(NewReplicates_All_Results,file=paste0(dir.output.results,"/", RefinementTAG,"_BeforeOpt_Replicates_",Name_Tag,".csv"),quote=F,row.names=F,col.names=T,sep=",")
  }
  
  ## Report ReScore
  ReScoreList<-list.files(dir.output,pattern = paste0("Report_",RefinementTAG,"_ReScore_"))
  ListFiles_ReScore<-paste0(dir.output,"/",ReScoreList)
  if(length(ReScoreList)>=1){
    l_reScore <- lapply(ListFiles_ReScore, fread, header = F,sep=';', stringsAsFactors = FALSE)
    NewReScore_All_Results <- rbindlist(l_reScore)
    colnames(NewReScore_All_Results)<-strsplit(ColumnNames[2,],split = ";",fixed = T)[[1]]
    
    NewReScore_All_Results<-NewReScore_All_Results %>% 
      mutate(ID_Analyte=as.character(ID_Analyte),
             ID_Rep=as.character(ID_Rep)) %>%
      left_join(MetaData_Analytes, by = c("ID_Analyte")) %>%
      left_join(MetaData_Replicate, by = c("ID_Rep"))
    
    Exponents<-c(9.5,4.5,2.5,0.5)
    NewReScore_All_Results<- NewReScore_All_Results %>%
      mutate(Skor=Similarity.Score^Exponents[1]*Library.dotp^Exponents[2]*Score.MassError^Exponents[3]*MPRA.Score^Exponents[4])
    
    write.table(NewReScore_All_Results,file=paste0(dir.output.results,"/",RefinementTAG,"_AfterOpt_Replicate_Score_",Name_Tag,".csv"),quote=F,row.names=F,col.names=T,sep=",")
  }
  
  ## Score Annotation
  
  ReScoreList<-list.files(dir.output,pattern = paste0("Report_",RefinementTAG,"_ReScore_"))
  ListFiles_ReScore<-paste0(dir.output,"/",ReScoreList)
  if(length(ReScoreList)>=1){
    l_reScore <- lapply(ListFiles_ReScore, fread, header = F,sep=';', stringsAsFactors = FALSE)
    ScoreAnnotations <- rbindlist(l_reScore)
    colnames(ScoreAnnotations)<-strsplit(ColumnNames[2,],split = ";",fixed = T)[[1]]
    
    Exponents<-c(9.5,4.5,2.5,0.5)
    ScoreAnnotations<-ScoreAnnotations %>% 
      mutate(ID_Analyte=as.character(ID_Analyte),
             ID_Rep=as.character(ID_Rep)) %>%
      left_join(MetaData_Analytes, by = c("ID_Analyte")) %>%
      left_join(MetaData_Replicate, by = c("ID_Rep")) %>%
      select(IsotopeLabelType, ProteinName, PeptideModifiedSequence,PrecursorCharge,IsDecoy, FileName,Similarity.Score, MPRA.Score, Library.dotp,Score.MassError) %>%
      mutate(Skor=Similarity.Score^Exponents[1]*Library.dotp^Exponents[2]*Score.MassError^Exponents[3]*MPRA.Score^Exponents[4])
    
    Annotations_PrecursorResults<-MetaData_PrecursorResults %>% 
      left_join(ScoreAnnotations, by = c("IsotopeLabelType", "ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "IsDecoy", "FileName"))%>%
      select(PrecursorResultLocator,
             Similarity.Score,MPRA.Score,Library.dotp,Score.MassError,Skor)%>%
      mutate(Similarity.Score=ifelse(is.na(Similarity.Score),"#N/A",Similarity.Score),
             MPRA.Score=ifelse(is.na(MPRA.Score),"#N/A",MPRA.Score),
             Library.dotp=ifelse(is.na(Library.dotp),"#N/A",Library.dotp),
             Score.MassError=ifelse(is.na(Score.MassError),"#N/A",Score.MassError),
             Skor=ifelse(is.na(Skor),"#N/A",Skor)) %>%
      rename(ElementLocator=PrecursorResultLocator,
             annotation_AvG_Similarity_Score=Similarity.Score,
             annotation_AvG_MPRA_Score=MPRA.Score,
             annotation_AvG_SpectralLibSim_Score=Library.dotp,
             annotation_AvG_MassError_Score=Score.MassError,
             annotation_AvG_Score=Skor)
    
    fwrite(Annotations_PrecursorResults,file=paste0(dir.output.results,"/",RefinementTAG,"_AnnotationsPrecursorResults_",Name_Tag,".csv"), sep=",", row.names=F)
  }
  
  ## No Result
  if (RefinementTAG %in% c("TR","PB","GR")) {
    Extract_Indices_FromFileNames<-function(ListFiles){
      ListFiles<-gsub(ListFiles,pattern = Name_Tag,replacement = "TAG")
      
      ListFiles<-data.frame(NameFiles=ListFiles) %>%
        mutate(Index=gsub(
          substr(NameFiles,start = str_locate(pattern = "TAG_",NameFiles)[2]+4,stop = nchar(as.character(NameFiles))),
          pattern=".csv",replacement = "")) %>% arrange(as.numeric(Index))
      
      #HighestIndex=max(as.numeric(ListFiles$Index),na.rm = T)
      #setdiff(1:HighestIndex,ListFiles$Index)
      #length(setdiff(1:HighestIndex,ListFiles$Index))
      
      return(ListIndeces<-ListFiles$Index)}
    
    #NoResult<-unique(setdiff(1:max(as.numeric(Extract_Indices_FromFileNames(ListFiles_ReScore))),Extract_Indices_FromFileNames(ListFiles_ReScore)))
    NoResult<-unique(setdiff(1:max(as.numeric(Extract_Indices_FromFileNames(ReScoreList))),Extract_Indices_FromFileNames(ReScoreList)))
    NoResult<-MetaData_Analytes %>% filter(ID_Analyte %in% NoResult) %>% select(PeptideModifiedSequence, IsDecoy,ID_Analyte) %>% distinct()
    write.table(NoResult,file=paste0(dir.output.results,"/",RefinementTAG,"_NoResult_",Name_Tag,".csv"),quote=F,row.names=F,col.names=T,sep=",")
  }
  
  
  ## Import Annotations into Skyline
  if (RefinementTAG %in% c("TR","GR")) {
    NewTransitions_All_Results_OnlyTransitionsInfo_2<-NewTransitions_All_Results_OnlyTransitionsInfo %>% 
      select(TransitionLocator) %>% mutate(Quantitative="TRUE") %>%
      #mutate(PrecursorMz=round(PrecursorMz,4),
      #       ProductMz=round(ProductMz,4)) %>%
      distinct()
    
    A<-MetaData_Transitions %>% 
      select(TransitionLocator) %>%
      distinct() %>%
      #mutate(PrecursorMz=round(PrecursorMz,4),ProductMz=round(ProductMz,4)) %>%
      left_join(NewTransitions_All_Results_OnlyTransitionsInfo_2, by = c("TransitionLocator")) %>%
      mutate(Quantitative=ifelse(is.na(Quantitative),"FALSE",Quantitative)) %>%
      rename(ElementLocator=TransitionLocator)
    
    fwrite(A,paste0(dir.output.results,"/",RefinementTAG,"_Transitions_Annotations_",Name_Tag,".csv"), sep=",", row.names = F)
    
  }
  
  ############ Cut_off Peak Boundaries

  determine_FDR_AvG <- function(NewReScore_All_Results){
    
    NewReScore_All_Results<-NewReScore_All_Results %>% 
      mutate(IsDecoy=gsub(IsDecoy,pattern = "False",replacement = 0)) %>%
      mutate(IsDecoy=gsub(IsDecoy,pattern = "True",replacement = 1))
    
    Q<-NewReScore_All_Results
    
    U<-Q
    
    U<-U %>% arrange(-Skor)
    
    U<-U %>% arrange(IsDecoy) %>% arrange(-Skor)
    
    U2<-data.frame(U$IsDecoy,count=ave(U$IsDecoy==U$IsDecoy,U$IsDecoy, FUN=cumsum))
    
    U3<-cbind(U,U2)
    U3$Num<-1:dim(U3)[1]
    U4<-U3 %>% mutate(FDR=ifelse(IsDecoy==1,count/Num*100,NA))
    
    FDR_1Percent<-as.numeric(U4 %>%  filter(FDR>=1) %>% filter(row_number()==1) %>% select(Skor))
    
    U5<-U4 %>% filter(IsDecoy==1)%>%select(Skor,FDR)
    
    
    b1=ggplot(U5,aes(x= Skor,y=FDR))+
      geom_line(size=1)+
      geom_vline(xintercept = FDR_1Percent,size=1,linetype=2,color="#FF0000")+
      geom_hline(yintercept = 1,linetype=2,size=1,color="black")+
      theme_bw()+scale_x_continuous(limits = c(0,1),breaks = seq(0,1,0.1))+
      theme(legend.position = "none",axis.text.x=element_text(angle = 45,vjust = 0.05))+
      annotate("text", x = 0.175,y=15, label = "FDR=1%",color = "#FF0000",angle = 90)+
      labs(y="FDR (%)",x="AvG score")
    
    
    Z=ggplot(Q,aes(x= Skor,fill=paste0(IsDecoy)))+
      geom_histogram(binwidth = 0.05,alpha=0.5,position="identity",color="black")+
      theme_bw()+
      geom_vline(xintercept = FDR_1Percent,linetype=2,size=1,color="#FF0000")+
      scale_x_continuous(breaks = seq(0,1,0.1))+
      scale_fill_manual(values = c("#5BBCD6", "#FF0000"))+
      theme(legend.position = "none", axis.text.x=element_text(angle = 45,vjust = 0.05))+
      annotate("text", size=3, x = 0.175,y=5000, label = "FDR=1%",color = "#FF0000",angle = 90)+
      theme_classic()+theme(legend.position = "none", axis.text.x=element_text(angle = 45,vjust = 0.05))+
      labs(y="Count",x="AvG score") 
    
    
    FDR_estimation <- grid.arrange(
      Z+labs(y="Number of peptides")+theme_bw()+
        theme(#axis.text.y = element_blank(),
          #axis.ticks.y = element_blank(),
          legend.position = "none",
          axis.text.x=element_text(angle = 45,vjust = 0.05)),
      
      b1+theme_bw()+
        theme(legend.position = "none", axis.text.x=element_text(angle = 45,vjust = 0.05)),
      ncol=1)
    
    ggsave(FDR_estimation,
           filename = paste0(dir.output.results,"/",RefinementTAG,"_FDR_plot_",Name_Tag,".pdf"),
           width = 6,height =9)
    
    return(FDR_1Percent=FDR_1Percent)}

  Cut_off_value = ifelse(length(unique(MetaData_Analytes$IsDecoy))>1,
                       determine_FDR_AvG(NewReScore_All_Results),
                       ifelse(NonZeroBaselineChromatogram==TRUE,0.61,0.1))

  cut_off_tag = gsub(paste0(ifelse(length(unique(MetaData_Analytes$IsDecoy))>1,
                       "FDR_Below1Percent_cutoffAvGScore_",
                       "Fixed_cutoffAvGScore_"), round(Cut_off_value,3)),pattern = "\\.",replacement = "pt")
  
  if (RefinementTAG %in% c("PB","GR")) {
    PeakBoundaries_Final<-NewReScore_All_Results%>% 
      filter(as.numeric(Library.dotp)>=0.7,
             as.numeric(Similarity.Score)>=0.85,
             as.numeric(Score.MassError)>=0.7,
             as.numeric(MPRA.Score)>=0.9,
             as.numeric(Skor)>=Cut_off_value)
    
    Num_Of_ValidReplicate_per_Peptide<-PeakBoundaries_Final %>%
      group_by(ID_Analyte,IsotopeLabelType,Comment,ProteinName,PeptideModifiedSequence, PrecursorCharge,IsDecoy) %>%
      summarise(n=n()) %>% ungroup() %>% rename(Num_Replicates=n)
    
    PeakBoundaries_Final<-PeakBoundaries_Final %>% left_join(Num_Of_ValidReplicate_per_Peptide,by = c("ID_Analyte", "IsotopeLabelType", "Comment", "ProteinName", "PeptideModifiedSequence", "PrecursorCharge", "IsDecoy")) %>%
      select(ID_Rep,ID_Analyte,IsotopeLabelType,FileName,ProteinName,PeptideModifiedSequence,PrecursorCharge,IsDecoy,Num_Replicates) %>%
      distinct() %>% 
      mutate(Keep="Keep")
    
    
    W<-NewPeakBoundaries_All_Results %>% left_join(PeakBoundaries_Final, by = c("ID_Rep", "ID_Analyte", "IsotopeLabelType", "FileName", "ProteinName", "PeptideModifiedSequence","PrecursorCharge", "IsDecoy")) %>% 
      mutate(left2=ifelse(Keep=="Keep",left,NA),right2=ifelse(Keep=="Keep",right,NA)) %>%
      mutate(left=left2,right=right2) %>% select(-left2,-right2) %>%
      mutate(left=ifelse(is.na(left),"#N/A",left),
             right=ifelse(is.na(right),"#N/A",right))
    
    Formatted_Filtered<-W %>%
      select(FileName,PeptideModifiedSequence,left,right,PrecursorCharge,IsDecoy,IsotopeLabelType) %>%
      rename(PrecursorIsDecoy=IsDecoy) %>%
      rename(MinStartTime=left,MaxEndTime=right) %>%
      distinct()
    
    if(UseHeavyPeakBoundariesForLight==TRUE) {Formatted_Filtered<-Formatted_Filtered %>% filter(IsotopeLabelType=="heavy")} ##Only for P100
    
    write.csv(Formatted_Filtered,file=paste0(dir.output.results,"/",RefinementTAG,"_NewPeakBoundaries_",Name_Tag,"_Formatted_Filtered_",cut_off_tag,".csv"),quote=F,row.names=F)}
  
  ## Only Scoring / No optimization
  NoOptimizationList<-list.files(dir.output,pattern = paste0("Report_",RefinementTAG,"_ScoreNonOptimized_"))
  ListFiles_NoOpt<-paste0(dir.output,"/",NoOptimizationList)
  if(length(NoOptimizationList)>=1){
    l_NoOpt <- lapply(ListFiles_NoOpt, fread, header = F,sep=';', stringsAsFactors = FALSE)
    NoOpt_All_Results <- rbindlist(l_NoOpt)
    colnames(NoOpt_All_Results)<-strsplit(ColumnNames[2,],split = ";",fixed = T)[[1]]
    
    NoOpt_All_Results<-NoOpt_All_Results %>% 
      mutate(ID_Analyte=as.character(ID_Analyte),
             ID_Rep=as.character(ID_Rep)) %>%
      left_join(MetaData_Analytes, by = c("ID_Analyte")) %>%
      left_join(MetaData_Replicate, by = c("ID_Rep"))
    
    Exponents<-c(9.5,4.5,2.5,0.5)
    NoOpt_All_Results<- NoOpt_All_Results %>%
      mutate(Skor=Similarity.Score^Exponents[1]*Library.dotp^Exponents[2]*Score.MassError^Exponents[3]*MPRA.Score^Exponents[4])
    
    write.table(NoOpt_All_Results,file=paste0(dir.output.results,"/",RefinementTAG,"_NoOptimization_Score_",Name_Tag,".csv"),quote=F,row.names=F,col.names=T,sep=",")
  }
  
  print("Format Results: Done!")
  
}

#' AvantGardeDIA_DB
#'
#' This function is AvantGardeDIA_DB
#' @param D.file.name Original CSV file containing all Metadata and DIA chromatograms.
#' @param RefinementWorkflow Type of refinement to be run.
#' Possible values are:
#' 'GlobalRefinement'  for 1) transition refinement, 2) peak boundaries refinement and 3)peak rescoring.
#' 'TransitionRefinement'  for 1) transition refinement and 2)peak rescoring.
#' 'PeakBoundariesRefinement'  for 1) peak boundaries refinement and 2)peak rescoring.
#' 'OnlyScoring'  for peak rescoring.
#' @param
#' @param ParamsFile File containing all user-defined parameters.
#' @keywords AvantGardeDIA
#' @export
#' @examples
#' library(AvantGardeDIATest5)
#' ParamsFile="C:/Users/Example/ParamsFile.R"
#' source(ParamsFile)
#' AvantGardeDIA_DB(D.file.name,RefinementWorkflow = "GlobalRefinement",ParamsFile)
AvantGardeDIA_DB<-function(D.file.name,RefinementWorkflow,ParamsFile){
  require(dplyr)
  require(tidyr)
  require(stringr)
  require(data.table)
  require(sqldf)
  require(foreach)
  require(doSNOW)
  require(GA, quietly = T)
  require(parallel)
  AvantGardeDIA_InChunks_DB(D.file.name,RefinementWorkflow)
  LaunchParallelTasks_ReadFromDB(ParamsFile,RefinementWorkflow)
  Read_AndFormatResults()
  print("AvantGardeDIA_DB: Done! Have a nice day!")
}


