args <- commandArgs()
print(args)
setwd(file.path(args[4]))
getwd()
dir.create(file.path(getwd(),"AvantGardeDIA"),showWarnings = F)
setwd("AvantGardeDIA/")

# Create AvantGardeDIA Parameters file
A<-"
Name_Tag='Name_Tag' #Short Name that will be added to all output files. (Must be in quotes)

# I- DataFormating Parameters
MinimalInitialNumberOfTransitions_Filter=5 # Peptides with less than N transitions will be removed and not included in the SQLite Database
RemoveTransitionBelowOrdinal_Filter=3  # Transition ordinal below or equal to this value will be removed (for example if Max_TransitionIonOrdinal=3 transitions b1,b2,b3,y1,y2,y3 are removed) 
RemoveSharedTransitionsBetweenLightAndHeavy_Filter=TRUE # Should be TRUE for DIA,FALSE for PRM.
ReadNumberOfLines=20000 # The CSV file is read in chunks. This parameter sets the number of lines in each chunk.

# II- Optimization Parameters
NonZeroBaselineChromatogram=TRUE # TRUE for data having high noise baseline level (e.g. SWATH data on Q-TOF), FALSE (noise-reduced data, e.g demultiplexed Overlap DIA Data acquired on Q-Exactive series)
MinimalNumberOfTransitionsAfterOptimization=4 # The optimization will find the best solution that will have at least this number of transitions.
KeepPeptidesWithLowerNumberOfTransitions=FALSE # Keep (or not) a peptide if it has a lower number of transitions than the limit above (recomended FALSE for DIA). The transition refinement refinement will not be done on these peptides but the peak boundary refienement and scoring will be performed.
TopN_RankedbyIntensity=6 # Number of most intense transitions in teh spectral library (N) among which at least n transitions need to be present.
MinimalNumberOfTransitionsAmongTopN=2 # Number of transitions (n) that need to be present among the N most intense transitions in the spectral library.
alpha=0.005
Beta=0.05
SpectralLibraryDotProduct_limit=0.9 # 
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

write.table(A,file = "AvG_Params.R",row.names = F,col.names = F,quote = F)

### Installing Latest version of the AvantGardeDIA R package
install_Latest_AvG_RPackage<-file.path(args[5],"Update_AvG_RPackage.r")
source(install_Latest_AvG_RPackage)


print("Create AvantGardeDIA parameters file: Done!")