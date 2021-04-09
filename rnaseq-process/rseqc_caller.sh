#!/bin/sh



version="0.1"


#######
#Usage
#sh rseqc_caller.sh -i bamfile -o outputfolder
#######

usage="

rseqqc_caller.sh

version: $version

Usage: sh rseqqc_caller.s -i bamfile -o outputfolder

Description: Use RSeQC to provide QC for bam files

Parameters:

	-i      Input BAM File
	-b      Genome Bed File
	-o      Output folder
	
"



if [ $# -eq 0 ]; then 
	printf "$usage";
	exit 1; 
fi


#####
#Functions needed
#####

realpath() {
    path=`eval echo "$1"`
    folder=$(dirname "$path")
    echo $(cd "$folder"; pwd)/$(basename "$path"); 
}



#######
#Input/Output
#######


#receive options
while getopts "i:o:" opt; do
  case ${opt} in
    i )
		infile=$(realpath $OPTARG)
		infilename=$(basename $infile)
      ;;
    b )
		bedfile=$(realpath $OPTARG)
		bedfilename=$(basename $bedfile)
      ;;	  
    o ) 
		outfolder=$(realpath $OPTARG)
      ;;  
    \? ) 
		printf "ERROR: Unknown options.\n\n$usage"
      ;;
	: ) printf "ERROR: Unknown options.\n\n$usage"
      ;;

  esac
done

#####
#Prerequisite
#####

rseqcfolder=/apps/anaconda3/bin/

#######
#QC begins
#######

mkdir -p $outfolder/


#bam stat
$rseqcfolder/bam_stat.py -i $infile > $outfolder/bam_stat.txt

#gene body
$rseqcfolder/geneBody_coverage.py -i $infile -r $bedfile -o $outfolder/

#infer exp
$rseqcfolder/infer_experiment.py -r $bedfile -i $infile > $outfolder/infer_experiment.txt

#inner distance
$rseqcfolder/inner_distance.py -i $infile -o $outfolder/ -r $bedfile

#read dist
$rseqcfolder/read_distribution.py -i $infile -r $bedfile > $outfolder/read_distribution.txt

#tin
$rseqcfolder/tin.py -i $infile -r $bedfile > $outfolder/read_distribution.txt
