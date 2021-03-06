% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/20190118_Avant-gardeDIA_v28_CleaningNonUsefulFunctions.R
\name{AvantGardeDIA_DB}
\alias{AvantGardeDIA_DB}
\title{AvantGardeDIA_DB}
\usage{
AvantGardeDIA_DB(D.file.name, RefinementWorkflow, ParamsFile)
}
\arguments{
\item{D.file.name}{Original CSV file containing all Metadata and DIA chromatograms.}

\item{RefinementWorkflow}{Type of refinement to be run.
Possible values are:
'GlobalRefinement'  for 1) transition refinement, 2) peak boundaries refinement and 3)peak rescoring.
'TransitionRefinement'  for 1) transition refinement and 2)peak rescoring.
'PeakBoundariesRefinement'  for 1) peak boundaries refinement and 2)peak rescoring.
'OnlyScoring'  for peak rescoring.}

\item{ParamsFile}{File containing all user-defined parameters.}
}
\description{
This function is AvantGardeDIA_DB
}
\examples{
library(AvantGardeDIATest5)
ParamsFile="C:/Users/Example/ParamsFile.R"
source(ParamsFile)
AvantGardeDIA_DB(D.file.name,RefinementWorkflow = "GlobalRefinement",ParamsFile)
}
\keyword{AvantGardeDIA}
