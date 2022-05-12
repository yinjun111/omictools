#!/bin/sh

version="1.11"


#v1.1 add samtools index
#v1.11 print arguments

#######
#Usage
#sh gtf_build.sh -g input.gtf -f genome.fasta -o outputfolder -a mart_export_gene.txt -t mart_export_tx.txt -n Human.B38.Ensembl88 -s Human
#######

usage="

gtf_build

version: $version

Usage: omictools gtf-build -g input.gtf -f genome.fasta -o outputfolder

Description: Build indexes for genome annotations for omictools rnaseq and other omics analyses

Parameters:

    -g      GTF file
    -f      Fastq file	
    -a      Gene annotation file
    -t      Tx annotation file	
    -o      Output folder

    -n      Name for the build
    -s      Species Name
	
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
while getopts "g:f:o:a:t:n:s:" opt; do
  case ${opt} in
    g )
		gtffile=$(realpath $OPTARG)
		gtffilename=$(basename $gtffile)
      ;;
    f )
		fastafile=$(realpath $OPTARG)
		fastafilename=$(basename $fastafile)
      ;;
    a )
		geneannofile=$(realpath $OPTARG)
      ;;
    t )
		txannofile=$(realpath $OPTARG)
      ;;	  
	o ) 
		outfolder=$(realpath $OPTARG)
      ;;
    n ) 
		buildname=$OPTARG
      ;;
    s ) 
		species=$OPTARG
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
	
	printf "gtf-build version $version running.\n" > $logfile
	
	echo $@ >> $logfile
	
	date >> $logfile
	
else
	printf "ERROR:No input file was found.\n"
fi


#programs to be used

omictoolsfolder=/apps/omictools/
ensembl2ucsc=$omictoolsfolder/gtf-build/ensembl2ucsc_caller.pl
gtf2anno=$omictoolsfolder/gtf-build/gtf2anno.pl
tx2anno=$omictoolsfolder/gtf-build/tx_anno.pl

star=/apps/STAR-2.7.8a/bin/Linux_x86_64/STAR
starfolder=/apps/STAR-2.7.8a/bin/Linux_x86_64/
rsem=/apps/RSEM-1.3.3/rsem-prepare-reference

samtools=/apps/samtools-1.12/bin/samtools

#######
#Input/Output
#######

#convert to UCSC format, by adding chr and removing non-primary chrs
printf "perl $ensembl2ucsc -i $gtffile -o $outfolder/${gtffilename/.gtf/_ucsc.gtf} -t gtf\n" | tee -a $logfile
perl $ensembl2ucsc -i $gtffile -o $outfolder/${gtffilename/.gtf/_ucsc.gtf} -t gtf


#convert genome fasta file
printf "perl $ensembl2ucsc -i $fastafile -o $outfolder/${fastafilename/.fa/_ucsc.fa} -t fasta\n" | tee -a $logfile
perl $ensembl2ucsc -i $fastafile -o $outfolder/${fastafilename/.fa/_ucsc.fa} -t fasta


#Build STAR Index
printf "$star --runMode genomeGenerate --runThreadN 6 --genomeDir $outfolder/$buildname\_STAR/ --genomeFastaFiles $outfolder/${fastafilename/.fa/_ucsc.fa}  --sjdbGTFfile $outfolder/${gtffilename/.gtf/_ucsc.gtf} >> $logfile 2>&1\n" | tee -a $logfile
$star --runMode genomeGenerate --runThreadN 6 --genomeDir $outfolder/$buildname\_STAR/ --genomeFastaFiles $outfolder/${fastafilename/.fa/_ucsc.fa} --sjdbGTFfile $outfolder/${gtffilename/.gtf/_ucsc.gtf} >> $logfile 2>&1

 
#Build RSEM Index
printf "$rsem --gtf $outfolder/${gtffilename/.gtf/_ucsc.gtf} --star --star-path $starfolder -p 6 $outfolder/${fastafilename/.fa/_ucsc.fa} $outfolder/$buildname\_STAR/$species\_RSEM >> $logfile 2>&1\n" | tee -a $logfile
$rsem --gtf $outfolder/${gtffilename/.gtf/_ucsc.gtf} --star --star-path $starfolder -p 6 $outfolder/${fastafilename/.fa/_ucsc.fa} $outfolder/$buildname\_STAR/$species\_RSEM >> $logfile 2>&1


#Build samtools .fai index
printf "$samtools faidx $outfolder/${fastafilename/.fa/_ucsc.fa}" | tee -a $logfile
$samtools faidx $outfolder/${fastafilename/.fa/_ucsc.fa}

#Build samtools .dict index
printf "$samtools dict $outfolder/${fastafilename/.fa/_ucsc.fa.dict}" | tee -a $logfile
$samtools dict $outfolder/${fastafilename/.fa/_ucsc.fa} > $outfolder/${fastafilename/.fa/_ucsc.dict}


#Generate gene annotation
printf "$gtf2anno -i $outfolder/${gtffilename/.gtf/_ucsc.gtf} -a $geneannofile -o $outfolder/${gtffilename/.gtf/_ucsc_gene_annocombo.txt}\n" | tee -a $logfile
perl $gtf2anno -i $outfolder/${gtffilename/.gtf/_ucsc.gtf} -a $geneannofile -o $outfolder/${gtffilename/.gtf/_ucsc_gene_annocombo.txt}

#Generate tx annotation
printf "perl $tx2anno -i $outfolder/${gtffilename/.gtf/_ucsc_gene_annocombo.txt} -a $txannofile -o  $outfolder/${gtffilename/.gtf/_ucsc_tx_annocombo.txt} | tee -a $logfile\n" | tee -a $logfile
perl $tx2anno -i $outfolder/${gtffilename/.gtf/_ucsc_gene_annocombo.txt} -a $txannofile -o  $outfolder/${gtffilename/.gtf/_ucsc_tx_annocombo.txt} 
 
#gtf2bed for RSeQC
#awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' Mus_musculus.GRCm38.88_ucsc.gtf | gtf2bed - > Mus_musculus.GRCm38.88_ucsc.bed

#Generate exon annotation
printf "perl $exon2anno -i $outfolder/${gtffilename/.gtf/_ucsc.gtf} -o  $outfolder/${gtffilename/.gtf/_ucsc_exon_annocombo.txt} \n" | tee -a $logfile
perl $exon2anno -i $outfolder/${gtffilename/.gtf/_ucsc.gtf} -o  $outfolder/${gtffilename/.gtf/_ucsc_exon_annocombo.txt} 
 

printf "gtf-build done.\n" >> $logfile
date >> $logfile
