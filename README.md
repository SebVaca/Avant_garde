# Avant-garde

Avant-garde (AvG) is a new tool to refine Data-Independent Acquisition (and Parallel Reaction Monitoring) by removing interfered transitions, adjusting integration boundaries and scoring peaks to control the FDR. Unlike other tools where MS runs are scored independently from each other, Avant-garde uses a novel data-driven scoring strategy. DIA signals are refined by learning from the data itself, using all measurements in all samples together to achieve the best optimization. Avant-garde evaluates the suitability of a peak to be used for quantification. It is capable of improving the selectivity, accuracy, and reproducibility of the quantification results in very complex biological matrices, reachng the same levels obtained with manual validation.

![AvG](http://drive.google.com/uc?export=view&id=1QOqZKxeFiQYlkPiX-07a4BMpROpuSmyh)

## Getting Started

* Avant-garde is an R-based package that works with [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view)'s reports.
* This github repository contains a [tutorial](HowToRunAvG.pdf) decribing how to install, how to run the package as an External tool in [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view) and how to import the results back into Skyline.
* An example data set in Skyline can be found [here](https://drive.google.com/open?id=1JVoak2CY0lFZ61RWP-PfUk1vCJh5pHxS).

### Installation of the R package
##### From tar.gz file
1. Download tar.gz file located [here](r_package).
2. Follow these [instructions](http://outmodedbonsai.sourceforge.net/InstallingLocalRPackages.html).

##### Using devtools
1. install the devtools package.
```
install.packages("devtools")
```
2. Load the devtools package.
```
library(devtools)
```
3. Install the package using 'install_github'.
```
install_github("SebVaca/Avant_garde_Publication", subdir="r_package/AvantGardeDIA")
```

### Installation of the Skyline external tool

1. Download zip file [here](skyline_external_tool/20190118_AvG_skylinetool/20190118_AvG_skylinetool.zip).
2. Follow instructions found in the [tutorial](HowToRunAvG.pdf).

## Built With

* [GA](https://cran.r-project.org/web/packages/GA/index.html) - R package for genetic algorithms.
* [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html) - Tidyverse collection of R packages for data science.

## Authors

* **Sebastian Vaca** - [Broad Institute](https://www.broadinstitute.org/proteomics), Cambridge, MA
* **Jacob D. Jaffe** - [Broad Institute](https://www.broadinstitute.org/proteomics), Cambridge, MA

## Acknowledgments
* [Skyline development team](https://skyline.ms/project/home/software/Skyline/begin.view)
* **Nick Schulman**, **Brendan MacLean** and **Michael MacCoss** - for their help integrating Avant-garde into Skyline.
