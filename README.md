# Avant-garde

Avant-garde is a new tool to refine Data-Independent Acquisition (and Parallel Reaction Monitoring) by removing interfered transitions, adjusting integration boundaries and scoring peaks to control the FDR. Unlike other tools where MS runs are scored independently from each other, Avant-garde uses a novel data-driven scoring strategy. DIA signals are refined by learning from the data itself, using all measurements in all samples together to achieve the best optimization. Avant-garde evaluates the suitability of a peak to be used for quantification. It is capable of improving the selectivity, accuracy, and reproducibility of the quantification results in very complex biological matrices, reachng the same levels obtained with manual validation.


![AvG](http://drive.google.com/uc?export=view&id=1QOqZKxeFiQYlkPiX-07a4BMpROpuSmyh)

## Getting Started

* Avant-garde is an R-based package that works with [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view)'s reports.
* This github repository contains a [tutorial](HowToRunAvG_20190125.pdf) decribing how to install, how to run the package as an External tool in [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view) and how to import the results back into Skyline.
* An example data set can be found [here](https://drive.google.com/open?id=1JVoak2CY0lFZ61RWP-PfUk1vCJh5pHxS).



## Built With

* [GA](https://cran.r-project.org/web/packages/GA/index.html) - R package for genetic algorithms
* [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html) - Tidyverse collection of R packages for data science



