#!/bin/sh

version="1.0"



#######
#Usage
#sh gtf_build.sh -g input.gtf -f genome.fasta -o outputfolder
#######

usage="

gtf_build

version: $version

Usage: omictools gtf-build -g input.gtf -f genome.fasta -o outputfolder

Description: Build indexes for genome annotations for omictools rnaseq and other omics analyses

Parameters:

	-g      GTF file
	-f      Fastq file	
	-o      Output folder
    -n      Name for the build
	
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
while getopts "g:f:o:n:" opt; do
  case ${opt} in
    g )
		gtffile=$(realpath $OPTARG)
		gtffilename=$(basename $fastqfile)
      ;;
    f )
		fastafile=$(realpath $OPTARG)
		fastafilename=$(basename fastafile)
      ;;
    o ) 
		outfolder=$(realpath $OPTARG)
      ;;
    n ) 
		buildname=$OPTARG
      ;;	  
    \? ) 
		printf "ERROR: Unknown options.\n\n$usage"
      ;;
	: ) printf "ERROR: Unknown options.\n\n$usage"
      ;;

  esac
done


if [ -f "$fastafile" ] || [ -f "$gtffile" ];then
	logfile=$outfolder/gtf-build_run.log
else
	echo "ERROR:No input file was found."
fi


#######
#Input/Output
#######

#convert to UCSC format, by adding chr and removing non-primary chrs
perl /apps/omictools/gtf-build/ensembl2ucsc_caller.pl -i Mus_musculus.GRCm38.88.gtf -o Mus_musculus.GRCm38.88_ucsc.gtf -t gtf

perl ~/Pipeline/omictools/gtf-build/ensembl2ucsc_caller.pl -i Mus_musculus.GRCm38.dna.primary_assembly.fa -o Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa -t fasta

#Build STAR Index
 /apps/STAR-2.7.8a/bin/Linux_x86_64/STAR --runMode genomeGenerate --runThreadN 6 --genomeDir /data/jyin/Databases/Genomes/Mouse/Mouse.B38.Ensembl88_STAR/ --genomeFastaFiles /data/jyin/Databases/Genomes/Mouse/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa  --sjdbGTFfile /data/jyin/Databases/Genomes/Mouse/Mus_musculus.GRCm38.88_ucsc.gtf
 
#Build RSEM Index
 /apps/RSEM-1.3.3/rsem-prepare-reference --gtf /data/jyin/Databases/Genomes/Mouse/Mus_musculus.GRCm38.88_ucsc.gtf --star --star-path /apps/STAR-2.7.8a/bin/Linux_x86_64/ -p 6 /data/jyin/Databases/Genomes/Mouse/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa /data/jyin/Databases/Genomes/Mouse/Mouse.B38.Ensembl88_STAR/Mouse_RSEM
 
#Generate gene annotation
perl /home/centos/Pipeline/omictools/gtf-build/gtf2anno.pl -i /data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc.gtf -a /data/jyin/Databases/Genomes/Mouse/mm10/mart_export_v88_mouse.txt -o /data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc_gene_annocombo.txt
 
#Generate tx annotation
perl /home/centos/Pipeline/omictools/gtf-build/tx_anno.pl -i /data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc_gene_annocombo.txt -a /data/jyin/Databases/Genomes/Mouse/mm10/mart_export_v88_mouse_tx.txt -o  /data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc_tx_annocombo.txt
 
 
#gtf2bed 
awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' Mus_musculus.GRCm38.88_ucsc.gtf | gtf2bed - > Mus_musculus.GRCm38.88_ucsc.bed
