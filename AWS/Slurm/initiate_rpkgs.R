

options(install.packages.compile.from.source = "always")

install.packages(c("argparser","ggplot2","Cairo","optparse","corrplot","ggfortify","reshape2","proj4","pheatmap","RColorBrewer","scales"),quiet=T,repos = "http://cran.us.r-project.org")

install.packages("/home/jyin/Programs/locfit_1.5-9.4.tar.gz",repos=NULL,type="source")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org",version="3.12")

BiocManager::install(c("DESeq2","limma","EnhancedVolcano"),ask =F)

