# Dockerizer container for Seurat
![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)

For consistency across locations, we have built out a quick RStudio container. To use with docker, use the command below to start an instance in your current working directory using whatever password you wish. If you want to use a specific directory in your local system, change the ./ to the appropriate path. 

## Docker command
```bash
docker run --rm -p 8787:8787 -e PASSWORD=$password -v ./:/home/rstudio alemenze/abrfseurat
```
To then load the instance, navigate to [localhost:8787](http://localhost:8787)
Enter rstudio as the username and the password you specified. 

Then enjoy playing in the data :) 

## Updating image
We can update the image as we wish, but unfortunately dockerhub no longer allows for automated builds without a paying account (as far as I know)- so we will need to manually build and push.
Alternatively, for one-off tests or temporary installs, the image includes remotes, devtools, and BiocManager for library installations. 

## Example library load
```r
#Basic manipulations
library(knitr)
library(tidyverse)
library(tibble)
library(dplyr)
library(hdf5r)
library(cowplot)
library(ggplot2)
library(ggalt)
library(ggpubr)
library(Matrix)
library(purrr)
library(reshape2)
library(RColorBrewer)
library(patchwork)
library(magrittr)

#scRNA core packages
library(Seurat)
library(S4Vectors)
library(scRNAseq)
library(SeuratWrappers)
library(SingleCellExperiment)
library(dittoSeq)
library(scDataviz)
library(DropletQC)
library(scCustomize)


#scRNA velocity, pseudotime, annotations, doublets
library(monocle3)
library(velocyto.R)
library(slingshot)
library(SingleR)
library(DoubletFinder)

#Other random packages for pseudobulk or things I like to use
library(apeglm)
library(scater)
library(magrittr)
library(pheatmap)
library(DESeq2)
library(edgeR)
library(EnhancedVolcano)
library(PCAtools)
library(clusterProfiler)
```

