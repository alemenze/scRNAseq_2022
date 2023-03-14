# Container for Seurat
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

For consistency across locations, we have built out a quick RStudio container. To use with docker, use the command below to start an instance in your current working directory using whatever password you wish. If you want to use a specific directory in your local system, change the ./ to the appropriate path. 

## Docker command
```bash
docker run --rm -p 8787:8787 -e PASSWORD=$password -v $localdir:/home/rstudio alemenze/abrfseurat
```
Change the password and localdir to what you wish to use for your settings. 
To then load the instance, navigate to [localhost:8787](http://localhost:8787)
Enter rstudio as the username and the password you specified. 

Note: if you try to load a directory with external connections you may get some funky errors. If you do, try a different local directory that doesnt have outgoing connections. 

Then enjoy playing in the data :) 

### Docker command breakdown:
```bash
docker run \ #starts docker and tells it to run
    --rm \ #removes the image when you are done
    -p 8787:8787 \ #defines the port you wish to use. the image uses an internal port of 8787, but the first value you can change to whatever local port you have free and wish to use. 
    -e PASSWORD=$password \ #recent updates demand adding a password. replace $password with whatever you wish- I usually just use "-e PASSWORD=test" since its easy
    -v $localdir:/home/rstudio \ #defines your local directory that will be mounted in the image. the default image directory is /home/rstudio, so we need to replace the $localdir with the path to your local directory with the data. 
    alemenze/abrfseurat #calls the specific dcoker image we have here!
```

## Singularity command
If you wish to run this on an HPC, you can also use this image there with singularity. Personally I do this in cases of memory issues when I exceed capacity of my local workstation- Ill write it up as a markdown file to be executed through the memory intensive steps (like integration), and save an R image to bring back to local machine. 

### Host your own RStudio container in singularity on an HPC
The HPC_execution.slurm file introduces how we can now run this as an interactive RStudio environment on your HPC, modified from [Rocker](https://rocker-project.org/use/singularity.html) to fit this project! You may need to modify the slurm file to fit your specific HPC- such as the partition you are requesting a node from, size of the node, time etc. The script will generate a random password and grab an available port from your node, but you can manually alter those if you please. From there you submit your job request.
```bash
sbatch HPC_execution.slurm
```
Once your job is accepted by slurm it will output 2 files. I probably should rename these... but the .err fill will output the important information for tunnelling. Read that file real quick and it should give you a command like `ssh -N -L 8787:${HOSTNAME}:${PORT} ${SINGULARITYENV_USER}@LOGIN-HOST` with the relevant information inputted. You will need to change @LOGIN-HOST to your respective cluster ssh access- for ours its just amarel.rutgers.edu for example.
From there on your local machine you should be able to go to [https://localhost:8787](https://localhost:8787), login with the username and password provided in the file, and have fun! This is set up to pull the singularity location from your CWD. 

#### Run just a script
This will run it in the current working directory off a file called "Processing.Rmd". 
```bash
singularity exec docker://alemenze/abrfseurat Rscript -e "rmarkdown::render('Processing.Rmd')"
```

#### CLI interactive singularity command
If your HPC uses slurm you would first have to request your nodes- below will demonstrate asking for 1 single full node without parallelization as our HPC is set up. 
```bash
srun -p 'p_lemenzad' --nodes=1 --time=03:00:00 --ntasks=1 --cpus-per-task=40 --mem=190000 --pty bash -i
```

I have not tested the interactivity, theoretically it will work in an CLI capacity but any sort of Rstudio will be dependent on your HPC structure. 
```bash
singularity shell --cleanenv docker://alemenze/abrfseurat
```
This has been roughly tested with singularity 3.5.2, but depending upon your install for singularity it may throw a user issue. 
Once in the interactive shell you can run things as you normally would for CLI access.

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
library(clustree)
library(Matrix.utils)

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
library(SoupX)

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

#Other additions
library(VAM)
library(gprofiler2)
library(glmGamPoi)
library(openxlsx)
library(kableExtra)
library(GSEABase)
library(ComplexHeatmap)
library(scSorter)
library(UCell)
```