#!/bin/bash

Rscript -e 'install.packages("devtools"); source("https://bioconductor.org/biocLite.R"); biocLite("BiocInstaller"); biocLite("BSgenome.Hsapiens.UCSC.hg19"); biocLite("GenomicRanges");  biocLite("RCytoscape"); biocList("VariantAnnotation"); biocList("graph"); devtools::install_github("devtools/gUtils");'
