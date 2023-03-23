# TCRBroom
A package to summarize and plot the pairwise and cross-cohort overlaps for TCR-seq data. 

## Introduction
TCR sequencing (TCR-seq) has been developed to =track the useful information of T cell clonality.
Similar to other high-throughput experiments, contamination can happen in several steps of TCR-seq, 
which creates artifacts in the data, leading to inaccurate or even biased results. 
Most of the existing methods assume “clean” TCR-seq data with minimum contamination as the starting point.
Few of them have specifically addressed the contamination issue in TCR-seq. 

We summarize the contamination into two sources, pairwise overlaps 
(two subjects sharing an excessively high proportion of TCR sequences) 
and cross-cohort overlaps (sequences detected in all or majority of the samples). 
The package can also generate heatmap and OncoPrint to present the pairwise overlaps 
and cross-cohort overlaps separately. 

## Installation
```
devtools::install_github("rx-li/TCRBroom")
```

## Input of the package and your TCR data
To use this package, we expect the TCR data of each sample is saved in ```.tsv``` format in which the first column contains
sequences' names and the other two columns contain the corresponding count and frequency rate.

Two inputs are required for the package:
* A meta table having two columns. The first column contains samples' names named as "sample_name" and the second column contains the corresponding 
subject ids named as . Each subject can have multiple samples. 
* Path to your TCR-seq data.

## Basic usage
```{r}
checkPairwise(meta, path)
checkCrosscohort(meta, path)
```
## Vignettes
A example workflow can be found at https://rx-li.github.io/tcrbroom_vignette.html.
