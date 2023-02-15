#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);

#CutAdapt+FASTQC+RSEM+STAR


########
#Interface
########


my $version="1.3";

#v0.2, add filter function to get files for PCA
#v0.3, removed -v, add -r implementation for local
#v0.31, solves screen envinroment problem
#v0.32, add CPM
#v0.4, major updates, to add parallel job, and --dev option
#v0.41, versioning
#v0.42, process multiqc
#v0.5, add expr-qc
#v0.51, --group for expr-qc
#v0.6, updated to AWS and v88
#v0.61, log for expr-qc
#v0.7, add rat annotation, change default I/O names
#v1.0, qstat to squeue, slurm tested
#v1.1, add AS calculation
#v1.2, add exon support. add submit_job function
#v1.3, add exon junc support
#v1.31, turn as T as default

my $usage="

rnaseq-merge
version: $version
Usage: omictools rnaseq-merge [parameters]

Description: Merge rnaseq-process folder to get summarized QC, count, TPM etc.

Parameters:

    --in|-i           Input folder(s), separated by \",\" [01.Process]
	
	#two ways of retrieving samples, either by names or by config files
    --samples|-s      Samples
    --config|-c       Configuration file

    --output|-o       Output folder [02.Merge]


    --tx|-t           Transcriptome
                        Currently support Human.B38.Ensembl88,Mouse.B38.Ensembl88,Rat.Rn6.Ensembl88

    --group|-g        Group name in config file for expr-qc plots [Group]

    --anno|-a         Add annotation

    --as              Prepare files for alternative splicing [T]

    --filter          Signal filter [auto]
                         automatically defined signal cutoff as
                           Count >= No. of samples * 5
                         or can define a count number
	
    --runmode|-r      Where to run the scripts, local, cluster or none [none]
	
	
    #Parameters for HPC

    --task            Number of tasks to be paralleled. By default by the number of commands in submitted jobs.

    --mem|-m          Deprecated for Slurm. Memory usage for each process, e.g. 100mb, 100gb [40gb]
	
    For each task, there are two ways of specifying the computing resource,
      but you can't mix --nodes and --ncpus together.
	A) by specifying number of nodes and process (default)
    --nodes           The value can be A) No. of nodes for each task
                                       B) Name of the node, e.g. n001.cluster.com                        
    --ppn             No. of processes for each task. By default [4]
                         Default for rnaseq-process
	
	B) by specifying the total number of cpus
    --ncpus           No. of cpus for each task for tasks can't use multiple nodes


";

#    --verbose|-v      Verbose

unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts



########
#Parameters
########

my $samples;
my $inputfolders="01.Process";

my $configfile;
my $outputfolder="02.Merge";
my $verbose=1;
my $tx;
my $group="Group";
my $as="T";
my $task;
my $ncpus=4;
my $ppn=6;
my $nodes;
my $mem;
my $filter="auto";
my $runmode="none";

my $dev=0; #developmental version

GetOptions(
	"in|i=s" => \$inputfolders,
	"samples|s=s" => \$samples,	
	"config|c=s" => \$configfile,
	"output|o=s" => \$outputfolder,
	"group|g=s" => \$group,
	"tx|t=s" => \$tx,	
	"as=s" => \$as,
	"filter=s" => \$filter,
	
	"task=s" => \$task,
	"ncpus=s" => \$ncpus,
	"ppn=s" => \$ppn,
	"nodes=s" => \$nodes,	
	"mem=s" => \$mem,	

	"runmode|r=s" => \$runmode,		

	"verbose|v" => \$verbose,
	
	"dev" => \$dev,		
);





########
#Prerequisites
########

my $omictoolsfolder="/apps/omictools/";

#adding --dev switch for better development process
if($dev) {
#	$omictoolsfolder="/home/jyin/Projects/Pipeline/omictools/";
#}
#else {
	#the tools called will be within the same folder of the script
	$omictoolsfolder=get_parent_folder(abs_path(dirname($0)));
}

#used programs
my $multiqc=find_program("/apps/anaconda3/bin/multiqc");
my $Rscript=find_program("/apps/R-4.0.2/bin/Rscript");


my $mergefiles="perl $omictoolsfolder/mergefiles/mergefiles_caller.pl";
my $parallel_job="perl $omictoolsfolder/parallel-job/parallel-job_caller.pl";
my $rnaseqmergefilter="$Rscript $omictoolsfolder/rnaseq-merge/rnaseq-merge_filter.R";
my $count2cpm="$Rscript $omictoolsfolder/rnaseq-merge/count2cpm.R";
my $process_multiqc="perl $omictoolsfolder/rnaseq-merge/process_multiqc_summary.pl";
my $expr_qc="$Rscript $omictoolsfolder/expr-qc/expr-qc.R";
my $tx2si="$Rscript $omictoolsfolder/rnaseq-merge/tx2si.R";
my $merge_exonjunc="perl $omictoolsfolder/rnaseq-merge/merge_exonjunc.pl";
my $anno_exonjunc="perl $omictoolsfolder/rnaseq-merge/anno_exonjunc.pl";



########
#default ouputs
########


my $genecountmerged="gene.results.merged.count.txt";
my $genetpmmerged="gene.results.merged.tpm.txt";
my $genefpkmmerged="gene.results.merged.fpkm.txt";
my $genecpmmerged="gene.results.merged.cpm.txt";

my $txcountmerged="tx.results.merged.count.txt";
my $txtpmmerged="tx.results.merged.tpm.txt";
my $txfpkmmerged="tx.results.merged.fpkm.txt";
my $txcpmmerged="tx.results.merged.cpm.txt";

my $txtpmsi="tx.results.merged.tpm.si.txt";

my $exoncountmerged="exon.results.merged.count.txt";
my $exoncountsi="exon.results.merged.count.si.txt";

my $exonjunccountmerged="exonjunc.results.merged.count.txt";
my $exonjuncanno="exonjuncanno.txt";
my $exonjunccountsi="exonjunc.results.merged.count.si.txt";

#Create folders

unless(defined $outputfolder && length($outputfolder)>0 ) {
	print STDERR "\nERROR: -o outputfolder needs to be defined without default value.\n\n";
	exit;
}

if(!-d $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);


my $scriptfolder="$outputfolder/scripts";
my $tempfolder="$outputfolder/temp";

if(!-d $scriptfolder) {
	mkdir($scriptfolder);
}

if(!-d $tempfolder) {
	mkdir($tempfolder);
}

my $logfile="$outputfolder/rnaseq-merge_run.log";

my $scriptfile1="$scriptfolder/rnaseq-merge_run1.sh";
my $scriptfile2="$scriptfolder/rnaseq-merge_run2.sh";

my $scriptlocalrun="$outputfolder/rnaseq-merge_local_submission.sh";
my $scriptclusterrun="$outputfolder/rnaseq-merge_cluster_submission.sh";


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");


print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "$multiqc version:", getsysoutput("$multiqc --version"),"\n";
print LOG "\n";


#test tx option

#my %tx2ref=(
#	"Human.B38.Ensembl88"=>"/home/jyin/Projects/Databases/Ensembl/v88/Human_STAR/Human_RSEM",
#	"Mouse.B38.Ensembl88"=>"/home/jyin/Projects/Databases/Ensembl/v88/Mouse_STAR/Mouse_RSEM",
#);


my %tx2ref=(
	"Human.B38.Ensembl88"=> { 
		"star"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl88_STAR",
		"rsem"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl88_STAR/Human_RSEM",
		"chrsize"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl88_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.88_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.88_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.88_ucsc_gene_annocombo.txt",
		"txanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.88_ucsc_tx_annocombo.txt",
		"exonanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.88_ucsc_exon_annocombo.txt",},
	"Mouse.B38.Ensembl88"=>{ 
		"star"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl88_STAR",
		"rsem"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl88_STAR/Mouse_RSEM",
		"chrsize"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl88_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc_gene_annocombo.txt",
		"txanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc_tx_annocombo.txt",
		"exonanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.88_ucsc_exon_annocombo.txt"},
	"Rat.Rn6.Ensembl88"=>{ 
		"star"=>"/data/jyin/Databases/Genomes/Rat/rn6/Rat.Rn6.Ensembl88_STAR",
		"rsem"=>"/data/jyin/Databases/Genomes/Rat/rn6/Rat.Rn6.Ensembl88_STAR/Rat_RSEM",
		"chrsize"=>"/data/jyin/Databases/Genomes/Rat/rn6/Rat.Rn6.Ensembl88_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Rat/rn6/Rattus_norvegicus.Rnor_6.0.dna.toplevel_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Rat/rn6/Rattus_norvegicus.Rnor_6.0.88_ucsc.gtf",
		"geneanno"=>"/data/jyin/Databases/Genomes/Rat/rn6/Rattus_norvegicus.Rnor_6.0.88_ucsc_gene_annocombo.txt",
		"txanno"=>"/data/jyin/Databases/Genomes/Rat/rn6/Rattus_norvegicus.Rnor_6.0.88_ucsc_tx_annocombo.txt",
		"exonanno"=>"/data/jyin/Databases/Genomes/Rat/rn6/Rattus_norvegicus.Rnor_6.0.88_ucsc_exon_annocombo.txt"}
);


########
#Process
########

print STDERR "\nomictools rnaseq-merge $version running ...\n\n" if $verbose;
print LOG "\nomictools rnaseq-merge $version running ...\n\n";


if(defined $tx2ref{$tx}) {
	print STDERR "Starting analysis using $tx.\n\n" if $verbose;
	print LOG "Starting analysis using $tx.\n\n";
}
else {
	print STDERR "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n" if $verbose;
	print LOG "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n";
	exit;
}



if($as eq "T") {
	print STDERR "\n--as T defined. Perform exon/exonjunc merging and SI calculation for Alternative Splicing analyses.\n\n" if $verbose;
	print LOG "\n--as T defined. Perform exon/exonjunc merging and SI calculation for Alternative Splicing analyses.\n\n";
}


#open config files to find samples

my @samples_array;
my %samples_hash; #check no dup samples


if(defined $samples && length($samples)>0) {
	foreach my $sample (split(",",$samples)) {
		#keep the order of samples
		push @samples_array, $sample;
		unless (defined $samples_hash{$sample}) {
			$samples_hash{$sample}++
		}
		else {
			print STDERR "ERROR:$sample is defined multiple times in $samples.\n";
			print LOG "ERROR:$sample is defined multiple times in $samples.\n";
			exit;
		}
	}
}
elsif(defined $configfile && length($configfile)>0) {
	#another way of getting samples
	
	$configfile=abs_path($configfile);
	
	open(IN,$configfile) || die "Error reading $configfile. $!";
	#first column should be sample, with a header

	my $fileline=0;
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		
		if($fileline>0) {
			my $sample=$array[0];
			push @samples_array, $sample;
			unless (defined $samples_hash{$sample}) {
				$samples_hash{$sample}++
			}
			else {
				print STDERR "ERROR:$sample is defined multiple times in $configfile.\n";
				print LOG "ERROR:$sample is defined multiple times in $configfile.\n";
				exit;
			}
		}
		
		$fileline++;
	}
	close IN;
}
else {
	#adding error message for no samples defined
	print STDERR "ERROR:Either --config or --samples needs to be defined.\n";
	print LOG "ERROR:Either --config or --samples needs to be defined.\n";

	exit;

}

print STDERR scalar(@samples_array)," samples identified: ",join(",",@samples_array),"\n\n" if $verbose;
print LOG scalar(@samples_array)," samples identified: ",join(",",@samples_array),"\n\n";

#----------------
#Find files for merging
#Examples
#JY_315_Isl.genes.results
#JY_315_Isl.isoforms.results

print STDERR "\nReading sample folders.\n" if $verbose;
print LOG "\nReading sample folders.\n";

my %sample2gene;
my %sample2tx;
my %sample2exon;
my %sample2exonjunc;
my %sample2folder;

foreach my $infolder (split(",",$inputfolders)) {
	#samples in different folders
	#if different folders have the same sample name, use the first found sample
	
	my @infiles=glob("$infolder/*");
	
	foreach my $file (@infiles) {
		if(-d $file) {
			my $samplename=basename($file);
			
			#see whether it is a defined sample
			if(defined $samples_hash{$samplename}) {
				print STDERR $samplename," is found in ",$file,".\n" if $verbose;
				print LOG $samplename," is found in ",$file,".\n" if $verbose;
				
				$sample2folder{$samplename}=abs_path($file);
				
				my @samplefiles=glob("$file/*");
				
				foreach  my $samplefile (@samplefiles) {
					#gene rsem
					if($samplefile=~/.genes.results/) {
						$sample2gene{$samplename}=abs_path($samplefile);
					}
					
					#tx rsem
					if($samplefile=~/.isoforms.results/) {
						$sample2tx{$samplename}=abs_path($samplefile);
					}
					
					if($as eq "T") {
						#exon count
						if($samplefile=~/featurecounts_exon.txt$/) {
							$sample2exon{$samplename}=abs_path($samplefile);
						}
						if($samplefile=~/exon.txt.jcounts$/) {
							$sample2exonjunc{$samplename}=abs_path($samplefile);
						}						
					}
					
				}
			}
		}
	}
}



print STDERR scalar(keys %sample2gene)," samples identified with gene results.\n" if $verbose;
print LOG scalar(keys %sample2gene)," samples identified with gene results.\n";

print STDERR scalar(keys %sample2tx)," samples identified with isoform results.\n" if $verbose;
print LOG scalar(keys %sample2tx)," samples identified with isoform results.\n";


if( scalar(@samples_array) != scalar(keys %sample2gene) || scalar(@samples_array) != scalar(keys %sample2tx)) {
	print STDERR "ERROR:Not all samples have gene or isoform results.\n\n";
	print LOG "ERROR:Not all samples have gene or isoform results.\n\n";
	exit;
}

if($as eq "T") {
	print STDERR scalar(keys %sample2exon)," samples identified with exon results.\n" if $verbose;
	print LOG scalar(keys %sample2exon)," samples identified with exon results.\n";
	
	if( scalar(@samples_array) != scalar(keys %sample2exon) ) {
		print STDERR "\nERROR:Not all samples have exon results. You need to run rnaseq-process with --as T.\n\n";
		print LOG "\nERROR:Not all samples have exon results. You need to run rnaseq-process with --as T.\n\n";
		exit;
	}
	
	print STDERR scalar(keys %sample2exonjunc)," samples identified with exon juncs results.\n" if $verbose;
	print LOG scalar(keys %sample2exonjunc)," samples identified with exon juncs results.\n";
	
	if( scalar(@samples_array) != scalar(keys %sample2exonjunc) ) {
		print STDERR "\nERROR:Not all samples have exon juncs results. You need to run rnaseq-process with --as T.\n\n";
		print LOG "\nERROR:Not all samples have exon juncs results. You need to run rnaseq-process with --as T.\n\n";
		exit;
	}	
}



#print OUT model to combine gene and tx results

system("cut -f 1 ".$sample2gene{$samples_array[0]}." > $tempfolder/genes.txt");
system("cut -f 1 ".$sample2tx{$samples_array[0]}." > $tempfolder/txs.txt");

if($as eq "T") {
	system("cut -f 1 ".$tx2ref{$tx}{"exonanno"}."> $tempfolder/exons.txt");
}

#print title for gene and txs
open(OUT,">$tempfolder/gene_title.txt") || die $!;
print OUT "Gene\t",join("\t",@samples_array),"\n";
close OUT;

open(OUT,">$tempfolder/tx_title.txt") || die $!;
print OUT "Tx\t",join("\t",@samples_array),"\n";
close OUT;

if($as eq "T") {
	open(OUT,">$tempfolder/exon_title.txt") || die $!;
	print OUT "Exon\t",join("\t",@samples_array),"\n";
	close OUT;
}



#get all files
my @genefiles;
my @txfiles;
my @exonfiles;
my @exonjuncfiles;

foreach my $sample (@samples_array) {
	if(defined $sample2gene{$sample}) {
		push @genefiles,$sample2gene{$sample};
	}
	else {
		print STDERR "ERROR: $sample gene.results not defined in $inputfolders.\n\n";
		print LOG "ERROR: $sample gene.results not defined in $inputfolders.\n\n";
		exit;
	}

	if(defined $sample2tx{$sample}) {
		push @txfiles,$sample2tx{$sample};
	}
	else {
		print STDERR "ERROR: $sample isoform.results not defined in $inputfolders.\n\n";
		print LOG "ERROR: $sample isoform.results not defined in $inputfolders.\n\n";
		exit;
	}
	
	if($as eq "T") {
		if(defined $sample2exon{$sample}) {
			push @exonfiles,$sample2exon{$sample};
		}
		else {
			print STDERR "ERROR: $sample exon results not defined in $inputfolders.\n\n";
			print LOG "ERROR: $sample exon results not defined in $inputfolders.\n\n";
			exit;
		}

		if(defined $sample2exonjunc{$sample}) {
			push @exonjuncfiles,$sample2exonjunc{$sample};
		}
		else {
			print STDERR "ERROR: $sample exon junc results not defined in $inputfolders.\n\n";
			print LOG "ERROR: $sample exon junc results not defined in $inputfolders.\n\n";
			exit;
		}
		
	}
	
}


########
#Print out commands, for local and cluster run
########

open(S1,">$scriptfile1") || die "Error writing $scriptfile1. $!";


#----------------
#Multiqc with specific folders


#file for multiqc folder search
open(OUT,">$tempfolder/samplefolders.txt") || die $!;
foreach my $sample (@samples_array) {
	print OUT $sample2folder{$sample},"\n";
}
close OUT;

print S1 "$multiqc -l $tempfolder/samplefolders.txt -o $outputfolder/multiqc;$process_multiqc -i $outputfolder/multiqc/multiqc_data/multiqc_general_stats.txt -o $outputfolder/multiqc/multiqc_data/multiqc_general_stats_rev.txt -c $configfile\n";

#------------------
#may need to implement annotation ...
#need to implement filter step for PCA ready file

#copy gene annotation
print S1 "cp ".$tx2ref{$tx}{"geneanno"}." $outputfolder/geneanno.txt;";
print S1 "cp ".$tx2ref{$tx}{"txanno"}." $outputfolder/txanno.txt;";


#Gene count
print S1 "$mergefiles -m $tempfolder/genes.txt -i ",join(",",@genefiles)," -l 5 -o $tempfolder/$genecountmerged\_wrongtitle;";
print S1 "tail -n +2 $tempfolder/$genecountmerged\_wrongtitle > $tempfolder/$genecountmerged\_notitle;";
print S1 "cat $tempfolder/gene_title.txt $tempfolder/$genecountmerged\_notitle > $outputfolder/$genecountmerged;";
print S1 "$count2cpm --count $outputfolder/$genecountmerged --cpm $outputfolder/$genecpmmerged;";
print S1 "rm $tempfolder/$genecountmerged\_wrongtitle;rm $tempfolder/$genecountmerged\_notitle;";
print S1 "\n";

#Gene tpm
print S1 "$mergefiles -m $tempfolder/genes.txt -i ",join(",",@genefiles)," -l 6 -o $tempfolder/$genetpmmerged\_wrongtitle -v;";
print S1 "tail -n +2 $tempfolder/$genetpmmerged\_wrongtitle > $tempfolder/$genetpmmerged\_notitle;";
print S1 "cat $tempfolder/gene_title.txt $tempfolder/$genetpmmerged\_notitle > $outputfolder/$genetpmmerged;";
print S1 "rm $tempfolder/$genetpmmerged\_wrongtitle;rm $tempfolder/$genetpmmerged\_notitle;";
print S1 "\n";

#Gene fpkm
print S1 "$mergefiles -m $tempfolder/genes.txt -i ",join(",",@genefiles)," -l 7 -o $tempfolder/$genefpkmmerged\_wrongtitle -v;";
print S1 "tail -n +2 $tempfolder/$genefpkmmerged\_wrongtitle > $tempfolder/$genefpkmmerged\_notitle;";
print S1 "cat $tempfolder/gene_title.txt $tempfolder/$genefpkmmerged\_notitle > $outputfolder/$genefpkmmerged;";
print S1 "rm $tempfolder/$genefpkmmerged\_wrongtitle;rm $tempfolder/$genefpkmmerged\_notitle;";
print S1 "\n";

#tx count
print S1 "$mergefiles -m $tempfolder/txs.txt -i ",join(",",@txfiles)," -l 5 -o $tempfolder/$txcountmerged\_wrongtitle -v;";
print S1 "tail -n +2 $tempfolder/$txcountmerged\_wrongtitle > $tempfolder/$txcountmerged\_notitle;";
print S1 "cat $tempfolder/tx_title.txt $tempfolder/$txcountmerged\_notitle > $outputfolder/$txcountmerged;";
print S1 "$count2cpm --count $outputfolder/$txcountmerged --cpm $outputfolder/$txcpmmerged;";
print S1 "rm $tempfolder/$txcountmerged\_wrongtitle;rm $tempfolder/$txcountmerged\_notitle;";
print S1 "\n";

#tx tpm
print S1 "$mergefiles -m $tempfolder/txs.txt -i ",join(",",@txfiles)," -l 6 -o $tempfolder/$txtpmmerged\_wrongtitle -v;";
print S1 "tail -n +2 $tempfolder/$txtpmmerged\_wrongtitle > $tempfolder/$txtpmmerged\_notitle;";
print S1 "cat $tempfolder/tx_title.txt $tempfolder/$txtpmmerged\_notitle > $outputfolder/$txtpmmerged;";
print S1 "rm $tempfolder/$txtpmmerged\_wrongtitle;rm $tempfolder/$txtpmmerged\_notitle;";
print S1 "\n";

#tx fpkm
print S1 "$mergefiles -m $tempfolder/txs.txt -i ",join(",",@txfiles)," -l 7 -o $tempfolder/$txfpkmmerged\_wrongtitle -v;";
print S1 "tail -n +2 $tempfolder/$txfpkmmerged\_wrongtitle > $tempfolder/$txfpkmmerged\_notitle;";
print S1 "cat $tempfolder/tx_title.txt $tempfolder/$txfpkmmerged\_notitle > $outputfolder/$txfpkmmerged;";
print S1 "rm $tempfolder/$txfpkmmerged\_wrongtitle;rm $tempfolder/$txfpkmmerged\_notitle;";
print S1 "\n";

#exon count

if($as eq "T") {

	print S1 "cp ".$tx2ref{$tx}{"exonanno"}." $outputfolder/exonanno.txt;";
	
	#exon count
	print S1 "$mergefiles -m $tempfolder/exons.txt -i ",join(",",@exonfiles)," -l 7 -o $tempfolder/$exoncountmerged\_wrongtitle;";
	print S1 "tail -n +2 $tempfolder/$exoncountmerged\_wrongtitle > $tempfolder/$exoncountmerged\_notitle;";
	print S1 "cat $tempfolder/exon_title.txt $tempfolder/$exoncountmerged\_notitle > $outputfolder/$exoncountmerged;";
	print S1 "rm $tempfolder/$exoncountmerged\_wrongtitle;rm $tempfolder/$exoncountmerged\_notitle;";
	print S1 "\n";
	
	#exon junc count
	print S1 "$merge_exonjunc -i ",join(",",@exonjuncfiles)," -o $outputfolder/$exonjunccountmerged;";
	print S1 "$anno_exonjunc -i $outputfolder/$exonjunccountmerged --tx $tx -o $outputfolder/$exonjuncanno;\n";
	
}


close S1;



#tasks, local 4, cluster by number of samples
unless(defined $task && length($task)>0) {
	if($runmode eq "cluster") {
		open(IN,$scriptfile1);
		my @script1cont=<IN>;
		my $scriptlines=scalar(@script1cont);
		close IN;
		print STDERR "\n$scriptlines commands detected in $scriptfile1. Use $scriptlines tasks in HPC.\n\n";
		$task=$scriptlines;
	}
	else {
		$task=4;
	}
}




##########
#filter for PCA
##########

open(S2,">$scriptfile2") || die "Error writing $scriptfile2. $!";

#gene
print S2 "$rnaseqmergefilter --count $outputfolder/$genecountmerged --fpkm $outputfolder/$genefpkmmerged --tpm $outputfolder/$genetpmmerged --filter $filter;";

my $genecountmerged_filtered="gene.results.merged.count.filtered.$filter.txt";
my $genecpmmerged_filtered="gene.results.merged.cpm.filtered.$filter.txt";
my $genetpmmerged_filtered="gene.results.merged.tpm.filtered.$filter.txt";

print S2 "$count2cpm --count $outputfolder/$genecountmerged_filtered --cpm $outputfolder/$genecpmmerged_filtered;";


#Expr-qc

print S2 "$expr_qc --input $outputfolder/$genetpmmerged_filtered --config $configfile --geneanno $outputfolder/geneanno.txt --out $outputfolder/expr-qc --group $group > $outputfolder/expr_qc_log.txt 2>&1;";

print S2 "\n";


#tx filter
print S2 "$rnaseqmergefilter --count $outputfolder/$txcountmerged --fpkm $outputfolder/$txfpkmmerged --tpm $outputfolder/$txtpmmerged --filter $filter;";

my $txcountmerged_filtered="tx.results.merged.count.filtered.$filter.txt";
my $txcpmmerged_filtered="tx.results.merged.cpm.filtered.$filter.txt";

print S2 "$count2cpm --count $outputfolder/$txcountmerged_filtered --cpm $outputfolder/$txcpmmerged_filtered;";

print S2 "\n";



########
#AS calculation
########

#to be implemented

if($as eq "T") {
	#tx SI
	#filter to have at least sum of n*0.1 tpm from all samples
	print S2 "$tx2si --gene $outputfolder/$genetpmmerged --tx $outputfolder/$txtpmmerged --anno ",$tx2ref{$tx}{"txanno"}," --out $outputfolder/$txtpmsi --filter auto*0.1\n";

	#exon SI
	#filter to have at least n*3 reads from all samples
	print S2 "$tx2si --gene $outputfolder/$genecountmerged --tx $outputfolder/$exoncountmerged --anno ",$tx2ref{$tx}{"exonanno"}," --out $outputfolder/$exoncountsi  --filter auto*3\n";


	#exon junc SI
	print S2 "$tx2si --gene $outputfolder/$genecountmerged --tx $outputfolder/$exonjunccountmerged --anno $outputfolder/$exonjuncanno --out $outputfolder/$exonjunccountsi  --filter auto*3\n";

	
}

close S2;


#######
#Run mode
#######

submit_job($scriptfile1,$scriptfile2);

close LOG;


########
#Functions
########

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}

sub build_timestamp {
	my ($now,$opt)=@_;
	
	if($opt eq "long") {
		$now=~tr/ /_/;
		$now=~tr/://d;
	}
	else {
		$now=substr($now,0,10);
	}
	
	return $now;
}

sub getsysoutput {
	my $command=shift @_;
	my $output=`$command`;
	$output=~tr/\r\n//d;
	return $output;
}

sub basename_short {
	my $filename=shift @_;
	
	my $basename;
	
	if($filename=~/([^\/]+)\/?$/) {
		$basename=$1;
	}
	
	return $basename;

}


sub find_program {
	my $fullprogram=shift @_;
	
	my $program;
	if($fullprogram=~/([^\/]+)$/) {
		$program=$1;
	}
	
	if(-e $fullprogram) {
		return $fullprogram;
	}
	else {
		my $sysout=`$program`;
		if($sysout) {
			my $location=`which $program`;
			return $location;
		}
		else {
			print STDERR "ERROR:$fullprogram or $program not found in your system.\n\n";
			exit;
		}
	}
}


sub get_parent_folder {
	my $dir=shift @_;
	
	if($dir=~/^(.+\/)[^\/]+\/?/) {
		return $1;
	}
}


sub submit_job {
	
	#lots of global variables. May need to be passed on in future implementation
	#$scriptlocalrun, $scriptclusterrun
	#$ppn, $task, $nodes, $runmode, LOG
	
	
	#only line needs to be changed for other scripts
	my @scripts_all=@_;


	open(LOUT,">$scriptlocalrun") || die "ERROR:can't write to $scriptlocalrun. $!";
	open(SOUT,">$scriptclusterrun") || die "ERROR:can't write to $scriptclusterrun. $!";
	
	#####
	#print out command for local parallel runs
	#####
	
	my $jobnumber=0;
	my $jobname="omictools-$timestamp";

	if($task eq "auto") {
		$jobnumber=0; #0 means as many as possible
	}
	else {
		$jobnumber=$task;
	}

	my @local_runs;
	my @script_names;

	foreach my $script (@scripts_all) {
		push @local_runs,"cat $script | parallel -j $jobnumber";

		if($script=~/([^\/]+)\.\w+$/) {
			push @script_names,$1;
		}
	}

	my $localcommand="screen -S $jobname -dm bash -c \"source ~/.bashrc;".join(";",@local_runs).";\"";

	print LOUT $localcommand,"\n";
	close LOUT;
	
	#####
	#print out command for cluster parallel runs
	#####
	
	my $clustercommand="$parallel_job -i ".join(",", @scripts_all)." -o $scriptfolder -n ".join(",",@script_names)." --tandem -t $task --env -r "; #changed here for none version


	if(defined $ppn && length($ppn)>0) {
		if(defined $nodes && length($nodes)>0) {
			$clustercommand.=" --nodes $nodes --ppn $ppn";		
		}
		else {
			$clustercommand.=" --nodes 1 --ppn $ppn";	
		}
	}
	else {
		$clustercommand.=" --ncpus $ncpus";	
	}

	if(defined $mem && length($mem)>0) {
		$clustercommand.=" -m $mem";	
	}

	print SOUT $clustercommand,"\n";
	close SOUT;



	if($runmode eq "none") {
		print STDERR "\nTo run locally, in shell type: sh $scriptlocalrun\n";
		print STDERR "To run in cluster, in shell type: sh $scriptclusterrun\n";
		
		print LOG "\nTo run locally, in shell type: sh $scriptlocalrun\n";
		print LOG "To run in cluster, in shell type: sh $scriptclusterrun\n";
	}
	elsif($runmode eq "local") {
		#local mode
		#implemented for local execution
		
		system("sh $scriptlocalrun");
		print LOG "sh $scriptlocalrun;\n\n";

		print STDERR "Starting local paralleled processing using $jobnumber tasks. To monitor process, use \"screen -r $jobname\".\n\n";
		print LOG "Starting local paralleled processing using $jobnumber tasks. To monitor process, use \"screen -r $jobname\".\n\n";
		
	}
	elsif($runmode eq "cluster") {
		#cluster mode
		#implement for HPC execution
		
		system("sh $scriptclusterrun");
		print LOG "sh $scriptclusterrun;\n\n";

		print STDERR "Starting cluster paralleled processing using $jobnumber tasks. To monitor process, use \"squeue\".\n\n";

	}
}
