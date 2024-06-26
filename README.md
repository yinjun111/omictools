# Omictools

Omictools suite for omics data processing in AWS HPC environment.

omictools

version: $version

Usage: omictools [tool] [parameters]


Parameters:

    ########
    #RNA-Seq
    ########
    rnaseq-process    RNA-seq QC, Align, and Counting for FASTQ files
    rnaseq-merge      Merge rnaseq-process results for downstream analyses
    rnaseq-de         Perform DE/AS analyses
    rnaseq-summary    Summarize RNA-Seq DE results

    prepare-merge     Generate a merge folder from raw count file for rnaseq-de/summary analysis

    rnaseq-var        RNA-seq variant calling pipeline

    ########
    #ChIP-Seq
    ########
    chipseq-process    ChIP-seq/ATAC-seq QC, Align, Peak calling and Counting for FASTQ files
    chipseq-merge      Merge chipseq-process results for downstream analyses
    chipseq-de         Perform DE analyses for ChIP-seq peaks
    chipseq-summary    Summarize ChIP-Seq results

    ########
    #DNA-Seq
    ########
    dnaseq-process    WGS/WES variant calling from FASTQ files using GATK workflow
    dnaseq-summary    Summarize variant calling results
    
    ########
    #AWS
    ########
    parallel-job      Batch job submission to HPC

    ########
    #Data Download
    ########
    bs-fastq          Download and merge FASTQ files from Basespace	
    geo-download      Download raw FASTQ files from GEO

    ########
    #Misc tools
    ########
    mergefiles        Use a model file to merge different files together
    text2excel        Build excel file using text file
    extract-command   Extract a command from omictools script file

