

options(install.packages.compile.from.source = "always")

install.packages(c("argparser","ggplot2","Cairo","optparse","corrplot","ggfortify","reshape2","proj4"),quiet=T,repos = "http://cran.us.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos = "http://cran.us.r-project.org")

BiocManager::install(c("DESeq2","limma","EnhancedVolcano"),ask =F)

