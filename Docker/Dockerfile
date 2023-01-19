FROM rocker/verse:4.2.2
LABEL description="Docker image for Seurat4."

RUN apt-get update && apt-get install -y \ 
    libhdf5-dev build-essential libxml2-dev libssl-dev libv8-dev libsodium-dev libglpk40 libgdal-dev \
    libboost-dev libomp-dev libbamtools-dev libboost-iostreams-dev libboost-log-dev libboost-system-dev libboost-test-dev libcurl4-openssl-dev \
    libz-dev libarmadillo-dev libhdf5-cpp-103

RUN R -e "install.packages(c('Seurat','hdf5r','dplyr','cowplot','knitr','slingshot','msigdbr','remotes','metap','devtools','R.utils','ggalt','ggpubr','BiocManager'),repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('SingleR','slingshot','scRNAseq','celldex','fgsea','multtest','scuttle','BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma', 'S4Vectors', 'SingleCellExperiment','SummarizedExperiment', 'batchelor', 'org.Mm.eg.db', 'AnnotationHub','scater','edgeR','apeglm','DESeq2','pcaMethods','clusterProfiler'))"
RUN R -e "remotes::install_github(c('satijalab/seurat-wrappers','kevinblighe/PCAtools','chris-mcginnis-ucsf/DoubletFinder','velocyto-team/velocyto.R'))"
RUN R -e "devtools::install_github(c('kevinblighe/EnhancedVolcano','cole-trapnell-lab/monocle3'))"

#Later additions
RUN R -e "BiocManager::install(c('edgeR','dittoSeq','scDataviz','glmGamPoi','GSEABase','ComplexHeatmap','UCell'))"
RUN R -e "devtools::install_github(c('hypercompetent/colorway','powellgenomicslab/DropletQC','samuel-marsh/scCustomize'))"
RUN R -e "install.packages(c('clustree','SoupX','gprofiler2','VAM','openxlsx','kableExtra','scSorter'))"

COPY ./Matrix.utils_0.9.8.tar.gz ./Matrix.utils_0.9.8.tar.gz
RUN R -e "install.packages('./Matrix.utils_0.9.8.tar.gz', repos=NULL, type='source')"