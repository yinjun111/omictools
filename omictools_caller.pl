#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

my $version="0.22";

#v0.11, human & mouse annotation added
#v0.12, rnaseq-var by GATK3 is added
#v0.13, some code cleaning, and add prepare-merge
#v0.14, add rat annotation, changed default i/o folders
#v0.2, slurm ready
#v0.21, add AS calculations for different functions
#v0.22, exon junction functions added

my $usage="

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


";


#    rnaseq-motif      RNA-seq TFBS motif finding pipeline
#    rnaseq-motif-summary  RNA-seq TFBS motif finding results summary	
#    motif-finder      Transcription factor binding motif prediction

unless (@ARGV) {
	print STDERR $usage;
	exit;
}


#functions

my ($command,@params)=@ARGV;

####
#check whether to use --dev version
####
my $params_list=join(" ",@params);

my $dev=0; #developmental version

if($params_list=~/--dev/) {
	$dev=1;
}


my $omictoolsfolder="/apps/omictools/";

#adding --dev switch for better development process
if($dev) {
	$omictoolsfolder="/home/jyin/Pipeline/omictools/";
}
else {
	#the tools called will be within the same folder of the script
	$omictoolsfolder=dirname(abs_path($0));
}



#####
#then call different scripts
#####


my $bs_fastq="perl $omictoolsfolder/bs-fastq/bs-fastq_caller.pl";
my $geo_download="perl $omictoolsfolder/geo-download/geo-download_caller.pl";

my $rnaseq_process="perl $omictoolsfolder/rnaseq-process/rnaseq-process_caller.pl";
my $rnaseq_merge="perl $omictoolsfolder/rnaseq-merge/rnaseq-merge_caller.pl";
my $rnaseq_de="perl $omictoolsfolder/rnaseq-de/rnaseq-de_caller.pl";
my $rnaseq_summary="perl $omictoolsfolder/rnaseq-summary/rnaseq-summary_caller.pl";

my $prepare_merge="perl $omictoolsfolder/prepare-merge/prepare-merge_caller.pl";

my $rnaseq_var="sh $omictoolsfolder/rnaseq-var/rnaseq-var_caller.sh";
my $rnaseq_var_summary="perl $omictoolsfolder/rnaseq-var/rnaseq-var_summary.pl";
my $rnaseq_motif="perl $omictoolsfolder/rnaseq-motif/rnaseq-motif_caller.pl";
my $rnaseq_motif_summary="perl $omictoolsfolder/rnaseq-motif-summary/rnaseq-motif-summary.pl";


my $chipseq_process="perl $omictoolsfolder/chipseq-process/chipseq-process_caller.pl";
my $chipseq_merge="perl $omictoolsfolder/chipseq-merge/chipseq-merge_caller.pl";
my $chipseq_de="perl $omictoolsfolder/chipseq-de/chipseq-de_caller.pl";
my $chipseq_summary="perl $omictoolsfolder/chipseq-summary/chipseq-summary_caller.pl";

my $dnaseq_process="perl $omictoolsfolder/dnaseq-process/dnaseq-process_caller.pl";

my $motif_finder="perl $omictoolsfolder/motif-finder/motif-finder_caller.pl";

my $mergefiles="perl $omictoolsfolder/mergefiles/mergefiles_caller.pl";
my $text2excel ="perl $omictoolsfolder/text2excel/text2excel.pl";
my $extract_command="perl $omictoolsfolder/misc/extract_command.pl";

my $parallel_job ="perl $omictoolsfolder/parallel-job/parallel-job_caller.pl";

my %commands2program=(
    "bs-fastq"=>$bs_fastq,
	"geo-download"=>$geo_download,
	
    "rnaseq-process"=>$rnaseq_process,
    "rnaseq-merge"=>$rnaseq_merge,
    "rnaseq-de"=>$rnaseq_de,
	"rnaseq-summary"=>$rnaseq_summary,

    "prepare-merge"=>$prepare_merge,

	"rnaseq-var"=>$rnaseq_var,
	"rnaseq-var-summary"=>$rnaseq_var_summary,
	"rnaseq-motif"=>$rnaseq_motif,
	"rnaseq-motif-summary"=>$rnaseq_motif_summary,
	
    "chipseq-process"=>$chipseq_process,
    "chipseq-merge"=>$chipseq_merge,
    "chipseq-de"=>$chipseq_de,
	"chipseq-summary"=>$chipseq_summary,	

    "dnaseq-process"=>$dnaseq_process,

	"motif-finder"=>$motif_finder,
    
	"mergefiles"=>$mergefiles,
    "text2excel"=>$text2excel,
	"extract-command"=>$extract_command,
	
	"parallel-job"=>$parallel_job,
);



if(defined $commands2program{$command}) {
	system($commands2program{$command}." ".join(" ",@params));
}
else {
	print STDERR "ERORR $command not found in omictools.\n\n";
	system("omictools");
}
