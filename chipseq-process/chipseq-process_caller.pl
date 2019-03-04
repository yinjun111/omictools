#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

#CutAdapt+FASTQC+RSEM+STAR


########
#Prerequisites
########

my $cutadapt="/home/jyin/.local/bin/cutadapt";
my $fastqc="/apps/FastQC/fastqc";

#star
my $star="/apps/STAR-master/bin/Linux_x86_64/STAR";

#homer
my $homer="/home/jyin/Programs/Homer/bin";
my $maketagdirectory="$homer/makeTagDirectory";
my $annotatepeaks="$homer/annotatePeaks.pl";
my $findpeaks="$homer/findPeaks";
my $makeucscfile="$homer/makeUCSCfile";
my $pos2bed="$homer/pos2bed.pl";

#UCSC
my $bedtobigbed="/apps/ucsc/bedToBigBed";

########
#Interface
########


my $version="0.1";

my $usage="

chipseq-processs
version: $version
Usage: sbptools chipseq-processs [parameters]

Description: 


Parameters:

    --config|-c       Configuration file
                           first column as sample name.
                           Controled headers include:
                               fastq1,fastq2 columns for fastqs
                               index for sample index used for folder names
    --output|-o       Output folder

    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84

    Parallel computating parameters
    --core            No. of cores or threads used by each task [4]
    --task            No. of tasks for parallel computing [4]

    --runmode|-r      Where to run the scripts, local, server or none [none]
    --verbose|-v      Verbose
	
	
";


unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts



########
#Parameters
########


my $configfile;
my $outputfolder;
my $verbose;
my $tx;
my $runmode="none";
my $core=4;

GetOptions(
	"config|c=s" => \$configfile,
	"output|o=s" => \$outputfolder,
	"tx|t=s" => \$tx,	
	"core=s" => \$core,	
	"runmode|r=s" => \$runmode,		
	"verbose|v" => \$verbose,
);

$outputfolder = abs_path($outputfolder);

if(!-e $outputfolder) {
	mkdir($outputfolder);
}


my $scriptfolder="$outputfolder/scripts";

if(!-e $scriptfolder) {
	mkdir($scriptfolder);
}


my $logfile="$outputfolder/chipseq-processs_run.log";
my $newconfigfile="$outputfolder/chipseq-processs_config.txt";

my $scriptfile1="$scriptfolder/chipseq-processs_run1.sh";
my $scriptfile2="$scriptfolder/chipseq-processs_run2.sh";

#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "$cutadapt version:", getsysoutput("$cutadapt --version"),"\n";
print LOG "$fastqc version:", getsysoutput("$fastqc --version"),"\n";
print LOG "$star version:", getsysoutput("$star --version"),"\n";
print LOG "homer veresion:",get_homer_version(),"\n";
print LOG "\n";


#test tx option

my %tx2ref=(
	"Human.B38.Ensembl84"=> { 
		"star"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR",
		"chrsize"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_homeranno.txt",},
	"Mouse.B38.Ensembl84"=>{ 
		"star"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR",
		"chrsize"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_homeranno.txt",}
);


if(defined $tx2ref{$tx}) {
	print STDERR "Starting analysis using $tx.\n\n" if $verbose;
	print LOG "Starting analysis using $tx.\n\n";
}
else {
	print STDERR "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n" if $verbose;
	print LOG "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n";
}


########
#Process
########

print STDERR "\nsbptools chipseq-processs running ...\n\n" if $verbose;
print LOG "\nsbptools chipseq-processs running ...\n\n";



###-------------------------
#From here, exactly the same with RNA-Seq ...		
#open config files to find fastqs

my %sample2fastq;
my %sample2indexname;
my %configattrs;
my %sample2chippedsample;
my %chippedsample2sample;
my %chippedsample2input;
my %inputsamples;

my $fileline=0;
my $newconfigfiletitle;
open(IN,$configfile) || die "Error reading $configfile. $!";
open(OUT,">$newconfigfile") || die "Error reading $newconfigfile. $!";

while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	if($fileline==0) {
		#record title into hash
		for(my $num=0;$num<@array;$num++) {
			#case sensitive match to key words, Sample, FASTQ1, FASTQ2, Index (for foldernames)
			$configattrs{uc $array[$num]}=$num;
		}
		
		#should always use first column as sample name and index to maintain continuous workflow
		
		#first column as sample
		
		#check title
		unless(defined $configattrs{"FASTQ1"}) {
			#SAMPLE and FASTQ1 have to be defined. The others are optional
			print STDERR "ERROR: FASTQ1 need to be defined in $configfile. Input config includes:",join(",",map {uc $_} @array),"\n";
			print LOG "ERROR: FASTQ1 need to be defined in $configfile. Input config includes:",join(",",map {uc $_} @array),"\n";
			exit;
		}
		else {
			print STDERR "Input config $configfile includes:",join(",",map {uc $_} @array),"\n\n" if $verbose;
			print LOG "Input config $configfile includes:",join(",",map {uc $_} @array),"\n\n";
		}
		
		#version 0.2, removed self defined INDEX
		
		if(defined $configattrs{"INDEX"}) {
			$newconfigfiletitle=join("\t","INDEX",map {uc $_} (splice @array,$configattrs{"INDEX"},1));
		}
		else {
			$newconfigfiletitle=join("\t","INDEX",map {uc $_} @array);
		}
		
		print OUT $newconfigfiletitle,"\n";
		
		print STDERR "New config $newconfigfile includes:",join(",",split("\t",$newconfigfiletitle)),"\n\n" if $verbose;
		print LOG "New config $newconfigfile includes:",join(",",split("\t",$newconfigfiletitle)),"\n\n";		
	
	}
	else {
		
		my $indexname;
		my $newconfigline;
		
		#need to change non word chars to words
		#if defined INDEX, use index, otherwise, use first column
		if(defined $configattrs{"INDEX"}) {
			$indexname=$array[$configattrs{"INDEX"}];
		}
		else {
			$indexname=$array[0]; #change first column as sample column in v0.2
			$indexname=~s/[^\w-]/_/g;
		}
		
		
		$sample2indexname{$array[0]}=$indexname;

		
		#fastq1
		push @{$sample2fastq{$indexname}},$array[$configattrs{"FASTQ1"}];
		
		#fastq2 (Read2)
		if(defined $configattrs{"FASTQ2"}) {
			push @{$sample2fastq{$indexname}},$array[$configattrs{"FASTQ2"}];
		}
		
		if(defined $configattrs{"INDEX"}) {
			$newconfigline=join("\t",$indexname,(splice @array,$configattrs{"INDEX"},1));
		}
		else {
			$newconfigline=join("\t",$indexname,@array);
		}

		#unique for chip-seq
		#INPUT and CHIPPEDSAMPLE as retained terms to match input
		if(defined $configattrs{"INPUT"} && defined $configattrs{"CHIPPEDSAMPLE"}) {
			if(defined $array[$configattrs{"CHIPPEDSAMPLE"}] && length($array[$configattrs{"CHIPPEDSAMPLE"}] )>0) {
				#chippedsamplegroup to match input
				$sample2chippedsample{$indexname}=$array[$configattrs{"CHIPPEDSAMPLE"}];
				$chippedsample2sample{$array[$configattrs{"CHIPPEDSAMPLE"}]}{$indexname}++;
				
				#find input sample
				if(defined $array[$configattrs{"INPUT"}] && length($array[$configattrs{"INPUT"}])>0 && $array[$configattrs{"INPUT"}]=~/Y/i) {
					$chippedsample2input{$array[$configattrs{"CHIPPEDSAMPLE"}]}=$indexname;
					$inputsamples{$indexname}++;
				}
			}
		}
		elsif(defined $configattrs{"INPUT"} || defined $configattrs{"CHIPPEDSAMPLE"}) {
			print STDERR "ERROR:INPUT and CHIPPEDSAMPLE are retained terms for ChIP-Seq. Both need to be defined.\nCurrent title includes:",join(",",sort keys %configattrs),"\n";
			exit;
		}
		
		print OUT $newconfigline,"\n";
		
	}
	
	$fileline++;
}

close IN;
close OUT;



#find input for chip-seq experiments

if(defined $configattrs{"INPUT"} && defined $configattrs{"CHIPPEDSAMPLE"}) {
	foreach my $chippedsample (sort keys %chippedsample2sample) {
		print STDERR "CHIPPEDSAMPLE $chippedsample includes:",join(",",sort keys %{$chippedsample2sample{$chippedsample}}),"\n" if $verbose;
		print LOG "CHIPPEDSAMPLE $chippedsample includes:",join(",",sort keys %{$chippedsample2sample{$chippedsample}}),"\n";
		
		if(defined $chippedsample2input{$chippedsample}) {
			print STDERR "CHIPPEDSAMPLE $chippedsample INPUT is:",$chippedsample2input{$chippedsample},"\n" if $verbose;
			print LOG "CHIPPEDSAMPLE $chippedsample INPUT is:",$chippedsample2input{$chippedsample},"\n";
		}
		else {
			print STDERR "ERROR:CHIPPEDSAMPLE $chippedsample doesn't have INPUT.\n\n";
			print LOG "ERROR:CHIPPEDSAMPLE $chippedsample doesn't have INPUT.\n\n";
			exit;
		}
	}
}


#----------------
#create folder
print STDERR scalar(keys %sample2indexname)," samples identified from $configfile, including:\nINDEX\tSAMPLE\n",join("\n",map {$sample2indexname{$_}."\t".$_} sort keys %sample2indexname),"\n\n" if $verbose;
print LOG scalar(keys %sample2indexname)," samples identified from $configfile, including:\nINDEX\tSAMPLE\n",join("\n",map {$sample2indexname{$_}."\t".$_} sort keys %sample2indexname),"\n\n" if $verbose;

#print STDERR scalar(keys %sample2fastq)," samples identified from $configfile, including:",join("\n",sort keys %sample2fastq),"\n\n" if $verbose;
#print LOG scalar(keys %sample2fastq)," samples identified from $configfile, including:",join("\n",sort keys %sample2fastq),"\n\n";


print STDERR "Make folders for different samples.\n\n" if $verbose;
print LOG "Make folders for different samples.\n\n";

foreach my $sample (sort keys %sample2fastq) {
	my $samplefolder="$outputfolder/$sample"; #edited to enable space in folder name
	if(!-e	$samplefolder) {
		mkdir $samplefolder;
	}
}

#----------------
#prepare fastq file for cutadapt
#assume truseq for now




my %sample2workflow;

#FASTQ files should all be .fastq.gz

if(defined $configattrs{"FASTQ2"}) {
	#PE
	print STDERR "Printing cutadapt PE script.\n\n" if $verbose;
	print LOG "Printing cutadapt PE script.\n\n";
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		
		my $fastq1=$sample2fastq{$sample}[0];
		my $fastq2=$sample2fastq{$sample}[1]; #
		
		my $fastq1trim="$samplefolder/".basename($fastq1);
		$fastq1trim=~s/\.fastq\.gz/_trimmed.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq1trim;

		my $fastq2trim="$samplefolder/".basename($fastq2);
		$fastq2trim=~s/\.fastq\.gz/_trimmed.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq2trim;
		
		my $cutadaptlog="$samplefolder/$sample\_cutadapt.log";

		$sample2workflow{$sample}.="$cutadapt -j 4 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $fastq1trim -p $fastq2trim $fastq1 $fastq2 > $cutadaptlog;";
	}
}
else {
	#SE
	print STDERR "Printing cutadapt SE script.\n\n" if $verbose;
	print LOG "Printing cutadapt SE script.\n\n";
	
	
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		
		my $cutadaptlog="$samplefolder/$sample\_cutadapt.log";
		
		my $fastq1=$sample2fastq{$sample}[0];
		
		my $fastq1trim="$samplefolder/".basename($fastq1);
		$fastq1trim=~s/\.fastq\.gz/_trimmed.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq1trim;

		$sample2workflow{$sample}.="$cutadapt -j 4 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o $fastq1trim $fastq1 > $cutadaptlog;";
	}
}
		

#----------------
#FASTQC for trimmed reads


if(defined $configattrs{"FASTQ2"}) {
	#PE
	print STDERR "Printing FASTQC PE script.\n\n" if $verbose;
	print LOG "Printing FASTQC PE script.\n\n";
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		
		$sample2workflow{$sample}.="$fastqc --nogroup -o $samplefolder -f fastq ".$sample2fastq{$sample}[2]." ".$sample2fastq{$sample}[3].";";
	}
}
else {
	#SE
	print STDERR "Printing FASTQC SE script.\n\n" if $verbose;
	print LOG "Printing FASTQC SE script.\n\n";
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		
		$sample2workflow{$sample}.="$fastqc --nogroup -o $samplefolder -f fastq ".$sample2fastq{$sample}[1].";";
	}
}	
		
		
#Till here, exactly the same with RNA-Seq ...		
###-------------------------



#----------------
#STAR

if(defined $configattrs{"FASTQ2"}) {
	#PE
	print STDERR "Printing STAR PE script.\n\n" if $verbose;
	print LOG "Printing STAR PE script.\n\n";

	#     rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name 
	#STAR --runThreadN 10 --genomeDir /home/jyin/Projects/Databases/Ensembl/v84/Human_UCSC_STAR/ --readFilesIn /home/jyin/Data/XuLab/Processed/ATAC-Seq/WT_1_R1.fastq.gz  /home/jyin/Data/XuLab/Processed/ATAC-Seq/WT_1_R2.fastq.gz --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix WT_1/WT_1_
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample/";

		#need to move to samplefolder, because star will genearte temp files 
		$sample2workflow{$sample}.="cd $samplefolder;$star --runThreadN $core --genomeDir ".$tx2ref{$tx}{"star"}." --readFilesIn ".$sample2fastq{$sample}[2]." ".$sample2fastq{$sample}[3]." --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix $samplefolder/$sample\_;";
	}
}
else {
	#SE
	print STDERR "Printing RSEM SE script.\n\n" if $verbose;
	print LOG "Printing RSEM SE script.\n\n";

	#     rsem-calculate-expression [options] upstream_read_file(s) reference_name sample_name 
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample/";
				
		$sample2workflow{$sample}.="$star --runThreadN $core --genomeDir ".$tx2ref{$tx}{"star"}." --readFilesIn ".$sample2fastq{$sample}[1]." --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix $samplefolder/$sample\_;";
	}
}


#----------------
#Peak Calling

foreach my $sample (sort keys %sample2fastq) {
	#tag
	
	#make sure each homer running in its own folder
	$sample2workflow{$sample}.="cd $outputfolder/$sample;";
	
	#makeTagDirectory
	#makeTagDirectory WT_1_TagDir -sspe WT_1/WT_1_Aligned.sortedByCoord.out.bam
	
	if(defined $configattrs{"FASTQ2"}) {
		$sample2workflow{$sample}.="$maketagdirectory $outputfolder/$sample/$sample\_TagDir -sspe $outputfolder/$sample/$sample\_Aligned.sortedByCoord.out.bam >& $outputfolder/$sample/$sample\_maketagdirectory.log;";
	}
	else {
		$sample2workflow{$sample}.="$maketagdirectory $outputfolder/$sample/$sample\_TagDir $outputfolder/$sample/$sample\_Aligned.sortedByCoord.out.bam >& $outputfolder/$sample/$sample\_maketagdirectory.log;";
	}
	
	unless(defined $inputsamples{$sample}) {

		#findPeaks
		#findPeaks WT_1_TagDir/ -style histone -o WT_1_Peaks.txt >& WT_1_Peaks.log
		
		if(defined $configattrs{"INPUT"} && defined $configattrs{"CHIPPEDSAMPLE"}) {
			my $inputsample=$chippedsample2input{$sample2chippedsample{$sample}};
			$sample2workflow{$sample}.="$findpeaks $outputfolder/$sample/$sample\_TagDir -i $outputfolder/$inputsample/$inputsample\_TagDir -style histone -o $outputfolder/$sample/$sample\_Peaks.txt >& $outputfolder/$sample/$sample\_findpeaks.log;";
		}
		else {
			$sample2workflow{$sample}.="$findpeaks $outputfolder/$sample/$sample\_TagDir -style histone -o $outputfolder/$sample/$sample\_Peaks.txt >& $outputfolder/$sample/$sample\_findpeaks.log;";
		}

		#Annotation
		#annotatePeaks.pl WT_1_Peaks.txt hg38 -gtf /home/jyin/Projects/Databases/Ensembl/v84/Homo_sapiens.GRCh38.84_ucsc.gtf -ann /home/jyin/Projects/Databases/Ensembl/v84/Human_Homer/Homo_sapiens.GRCh38.84_ucsc_anno.txt > WT_1_Peaks_anno2.txt

		my $genomeversion;
		
		if($tx=~/Human.B38/) {
			$genomeversion="hg38";
		}
		elsif($tx=~/Mouse.B38/) {
			$genomeversion="mm10";
		}
	
		$sample2workflow{$sample}.="$annotatepeaks $outputfolder/$sample/$sample\_Peaks.txt $genomeversion -gtf ".$tx2ref{$tx}{"gtf"}." -ann ".$tx2ref{$tx}{"homeranno"}." > $outputfolder/$sample/$sample\_Peaks_Anno.txt 2> $outputfolder/$sample/$sample\_annotatepeaks.log;";
	}	
	
	#make bed file
	$sample2workflow{$sample}.="$pos2bed $outputfolder/$sample/$sample\_Peaks.txt > $outputfolder/$sample/$sample\_Peaks.bed;";
	
	#bedtobigbed
	$sample2workflow{$sample}.="cat $outputfolder/$sample/$sample\_Peaks.bed | grep -v \"#\" | sort -k1,1 -k2,2n > $outputfolder/$sample/$sample\_Peaks_sorted.bed;$bedtobigbed $outputfolder/$sample/$sample\_Peaks_sorted.bed ".$tx2ref{$tx}{"chrsize"}." $outputfolder/$sample/$sample\_Peaks_sorted.bb;";
	
	#makeUCSCfile
	$sample2workflow{$sample}.="$makeucscfile $outputfolder/$sample/$sample\_TagDir -o $outputfolder/$sample/$sample.bw -bigWig ".$tx2ref{$tx}{"chrsize"}." > $outputfolder/$sample/$sample\_trackinfo.txt;";
	
}


########
#Print out commands, for local and server run
########


open(S1,">$scriptfile1") || die "Error writing $scriptfile1. $!";
open(S2,">$scriptfile2") || die "Error writing $scriptfile2. $!";

foreach my $sample (sort keys %sample2workflow) {
	if(defined $configattrs{"INPUT"} && defined $configattrs{"CHIPPEDSAMPLE"}) {
		if(defined $inputsamples{$sample}) {
			print S1 $sample2workflow{$sample},"\n";
		}
		else {
			print S2 $sample2workflow{$sample},"\n";
		}
	}
	else {
		print S1 $sample2workflow{$sample},"\n";
	}
}

close S1;
close S2;


#local mode
if(defined $configattrs{"INPUT"} && defined $configattrs{"CHIPPEDSAMPLE"}) {
	print STDERR "To run locally, in shell type: sh $scriptfile1;sh $scriptfile2\n\n";
	print LOG "To run locally, in shell type: sh $scriptfile1;sh $scriptfile2\n\n";
}
else {
	print STDERR "To run locally, in shell type: sh $scriptfile1\n\n";
	print LOG "To run locally, in shell type: sh $scriptfile1\n\n";
}


#whether to run it instantly
if($runmode eq "local") {
	system("sh $scriptfile1");
}


close LOG;
########
#Functions
########

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}

sub getsysoutput {
	my $command=shift @_;
	my $output=`$command`;
	$output=~tr/\r\n//d;
	return $output;
}


sub get_homer_version {

#process this
#/home/jyin/Programs/Homer/config.txt

	return("homer   v4.10.3");

}









