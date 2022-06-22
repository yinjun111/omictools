#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);



########
#Interface
########


my $version="1.2";

#version 0.2a, add r version log
#v0.3 add runmode
#v0.31, solves screen envinroment problem
#v0.4, add server option, updating r script
#v0.5, use R 4.0.2, add comparison argument for multiple comparisons
#v0.51, versioning
#v0.6, major updates planned for R4.0, comparisons of multiple groups. turn off as
#v0.61, adding -s for --comparisons
#v0.7, v88 and AWS
#v0.8, support complicated GLM
#v0.81, add idf and cc options
#v0.9, add rat annotation, change default I/O names
#v1.0, fix for slurm envir, slurm tested
#v1.1, add AS calculation
#v1.2, add exon AS. add submit_job

my $usage="

rnaseq-de
version: $version
Usage: omictools rnaseq-de [parameters]

Description: Differential Expression (DE) and Alternative Splicing (AS) analyses. The DE analysis script works for most of the counting based data, e.g. RNA-Seq, ChIP-Seq, ATAC-Seq


Mandatory Parameters:
    --in|-i           Input folder from rnaseq-merge [02.Merge]

    --output|-o       Output folder [03.DE]
                         Changes since v0.5. If your output folder in -o DE/,
                         the comparison results will be save in DE/treatment_vs_reference folder

    --config|-c       Configuration file match the samples in the rnaseq-merge folder
                           first column as sample name.

    --tx|-t           Transcriptome
                        Currently support Human.B38.Ensembl88,Mouse.B38.Ensembl88,Rat.Rn6.Ensembl88

    --formula|-f      Formula for GLM [~Group]
                          the last factor of the formula is used for comparison

    #if you have multiple comparisons to perform in a project
    --comparisons|-s  Tab delimited file with first column as treatment groups
                        and second column as reference groups for multiple pairwise comparisons

    #if you only have one comparison
    --treatment       Treatment group name
    --reference       Reference group name

						
Optional Parameters:
    --pmethod         DESeq2 method, default as Wald test [DESeq2-Wald]
    --qmethod         Multiple testing correction method [BH]

    --useallsamples   Use all samples in config file for 
                           gene dispersion calcultation [F]

    --filter          Signal filter [auto]
                         automatically defined signal cutoff as
                           Count >= No. of samples * 5
                         or can define a count number

    --fccutoff        Log2 FC cutoff [1]
    --qcutoff         Corrected P cutoff [0.05]

    #DESeq2 specific
    --independentfiltering   Apply independentfiltering in DESeq2 [T]
    --cookscutoff            Apply Cook\'s cutoff  in DESeq2 [T]

Alternative Splicing Analyses Parameters:
    --as              Run Tx DE/AS tests for alternative splicing analyses for tx and exon [F]

    --sifccutoff      Log2 FC cutoff, Log2(1.2)=0.263 [0.263]
    --siqcutoff       Corrected P cutoff [0.05]


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



    --runmode|-r      Where to run the scripts, local or none [none]
                           Note: cluster doesn't work for now for rnaseq-de, due to R pkg dependency


";

#    --verbose|-v      Verbose

#add --comparison argument to read comparisons.txt for pairwise comparisons...
#may need to turn off -o, and use comparison name as output folder name to be compatible with gsea-gen


#R parameters
my $rparams="
  -i, --in IN			Expr input file
  -a, --anno ANNO			Annotation file
  -o, --out OUT			Output file
  -f, --formula FORMULA			DESeq formula
  -t, --treat TREAT			treatment name
  -r, --ref REF			reference name
  --fccutoff FCCUTOFF			Log2 FC cutoff [default: 1]
  -q, --qcutoff QCUTOFF			qcutoff [default: 0.05]
  -p, --pmethod PMETHOD			Method used in DESeq2 [default: Wald]
  --qmethod QMETHOD			FDR Method [default: BH]
  --filter FILTER			Count filter for all the exps [default: 10]
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
	
my $inputfolder="02.Merge";	
my $configfile;
my $outputfolder="03.DE";

my $formula="~Group";
my $comparisons;
my $treatment;
my $reference;

my $pmethod="Wald";
my $qmethod="BH";
my $useallsamples="F";
my $filter="auto";

my $as="F";
my $fccutoff=1;
my $qcutoff=0.05;
my $sifccutoff=0.263;
my $siqcutoff=0.05;
my $independentfiltering="T";
my $cookscutoff="T";


my $verbose=1;
my $task;
my $ncpus=2;
my $ppn=6;
my $nodes;
my $mem;

my $dev=0; #developmental version 


my $tx;
my $runmode="none";

GetOptions(
	"in|i=s" => \$inputfolder,
	"config|c=s" => \$configfile,
	"out|o=s" => \$outputfolder,
	"formula|f=s" => \$formula,
	
	"comparisons|s=s" => \$comparisons,
	"treatment=s" => \$treatment,
	"reference=s" => \$reference,
	
	"pmethod=s" => \$pmethod,
	"qmethod=s" => \$qmethod,
	"useallsamples=s" => \$useallsamples,
	"filter=s" => \$filter,
	
	"fccutoff=s" => \$fccutoff,
	"qcutoff=s" => \$qcutoff,
	"sifccutoff=s" => \$sifccutoff,
	"siqcutoff=s" => \$siqcutoff,
	
	"as=s" => \$as,
	"independentfiltering=s" => \$independentfiltering,
	"cookscutoff=s" => \$cookscutoff,	
	
	"tx|t=s" => \$tx,	
	"runmode|r=s" => \$runmode,		
	"task=s" => \$task,	
	"ppn=s" => \$ppn,
	"nodes=s" => \$nodes,
	"ncpus=s" => \$ncpus,	
	"verbose|v" => \$verbose,

	"dev" => \$dev,		
);

########
#Prerequisites
########
my $R=find_program("/apps/R-4.0.2/bin/R");
my $Rscript=find_program("/apps/R-4.0.2/bin/Rscript");


my $omictoolsfolder="/apps/omictools/";

#adding --dev switch for better development process
if($dev) {
#	$omictoolsfolder="/home/centos/Projects/Pipeline/omictools/";
#}
#else {
	#the tools called will be within the same folder of the script
	$omictoolsfolder=get_parent_folder(dirname(abs_path($0)));
}


#omictools
my $descript="$Rscript $omictoolsfolder/rnaseq-de/de_test_caller.R";
my $summarize_txsi="perl $omictoolsfolder/rnaseq-de/summarize_txsi.pl";
my $summarize_exonsi="perl $omictoolsfolder/rnaseq-de/summarize_exonsi.pl";
my $mergefiles="perl $omictoolsfolder/mergefiles/mergefiles_caller.pl";
my $parallel_job="perl $omictoolsfolder/parallel-job/parallel-job_caller.pl";



#my $r=find_program("/apps/R-3.4.1/bin/R");
#my $Rscript=find_program("/apps/R-3.4.1/bin/Rscript");


#######
#Input/Output
#######

$inputfolder = abs_path($inputfolder);

unless(defined $outputfolder && length($outputfolder)>0 ) {
	print STDERR "\nERROR: -o outputfolder needs to be defined without default value.\n\n";
	exit;
}

if(!-d $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);

my $scriptfolder="$outputfolder/scripts";

if(!-d $scriptfolder) {
	mkdir($scriptfolder);
}


my $logfile="$outputfolder/rnaseq-de_run.log";

my $rlogfile="$outputfolder/rnaseq-de_r_env.log";

my $scriptfile1="$scriptfolder/rnaseq-de_run1.sh";
my $scriptfile2="$scriptfolder/rnaseq-de_run2.sh";

my $scriptlocalrun="$outputfolder/rnaseq-de_local_submission.sh";
my $scriptclusterrun="$outputfolder/rnaseq-de_cluster_submission.sh";



#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");


print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";

print STDERR "\nomictools rnaseq-de $version running ...\n\n" if $verbose;
print LOG "\nomictools rnaseq-de $version running ...\n\n";



#Report R package version here !!!

open(RLOG,"|$R --no-restore --no-save --slave") || die $!;
select RLOG;
print << "CODE";

rinfo<-c()

Rversion<-getRversion()

rinfo<-rbind(rinfo,c("R",R.Version()\$version.string))
rinfo<-rbind(rinfo,c("Rscript","$Rscript"))
rinfo<-rbind(rinfo,c("R library",paste(.libPaths(), collapse=",")))

for (package in c("BiocManager","DESeq2","argparser","ggplot2","EnhancedVolcano")) {
	rinfo<-rbind(rinfo,c(package,packageDescription(package,fields="Version")))
}

colnames(rinfo)<-c("Package","Version")

write.table(file="$rlogfile",rinfo,sep="\t",quote=F,row.names=F)

q()
CODE
close RLOG;


##test tx option

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
	print STDERR "\n--as T defined. Perform Alternative Splicing analyses.\n\n" if $verbose;
	print LOG "\n--as T defined. Perform Alternative Splicing analyses.\n\n";
}


#RNA-Seq

#my $genecountmerged="gene.results.merged.count.txt";
#my $txcountmerged="tx.results.merged.count.txt";
#my $genederesult="gene.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt";
#my $asresult="tx.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt";


my %rnaseq2files=(
	"gene"=> { 
		"count"=> "gene.results.merged.count.txt",
		"selected"=> "gene.results.merged.count.selected.txt",
		"result"=> "gene.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"resultanno"=> "gene.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
	},
	"tx"=> { 
		"count"=> "tx.results.merged.count.txt",
		"si"=> "tx.results.merged.tpm.si.txt",		
		"selected"=> "tx.results.merged.count.selected.txt",
		"siselected"=> "tx.results.merged.tpm.si.selected.txt",				
		"result"=> "tx.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"siresult"=> "tx.results.merged.tpm.si.Lm.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"resultanno"=> "tx.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
		"siresultanno"=> "tx.results.merged.tpm.si.Lm.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
		"siresultsummary"=> "tx.results.merged.tpm.si.summary.txt",
		"siresultsummaryanno"=> "tx.results.merged.tpm.si.summary.anno.txt",
	},
	"exon"=> { 
		"count"=> "exon.results.merged.count.txt",
		"si"=> "exon.results.merged.count.si.txt",		
		"selected"=> "exon.results.merged.count.selected.txt",
		"siselected"=> "exon.results.merged.count.si.selected.txt",				
		"result"=> "exon.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"siresult"=> "exon.results.merged.count.si.Lm.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"resultanno"=> "exon.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
		"siresultanno"=> "exon.results.merged.count.si.Lm.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
		"siresultsummary"=> "exon.results.merged.count.si.summary.txt",
		"siresultsummaryanno"=> "exon.results.merged.count.si.summary.anno.txt",
	},
	"exonjunc"=> { 
		"count"=> "exonjunc.results.merged.count.txt",
		"si"=> "exonjunc.results.merged.count.si.txt",		
		"anno"=> "exonjuncanno.txt",				
		"selected"=> "exonjunc.results.merged.count.selected.txt",
		"siselected"=> "exonjunc.results.merged.count.si.selected.txt",				
		"result"=> "exonjunc.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"siresult"=> "exonjunc.results.merged.count.si.Lm.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"resultanno"=> "exonjunc.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
		"siresultanno"=> "exonjunc.results.merged.count.si.Lm.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
		"siresultsummary"=> "exonjunc.results.merged.count.si.summary.txt",
		"siresultsummaryanno"=> "exonjunc.results.merged.count.si.summary.anno.txt",
	}		
);


#ChIP/ATAC-Seq

my %chip2files=(
	"gene"=> { 
		"count"=> "gene.results.merged.count.txt",
		"selected"=> "gene.results.merged.count.selected.txt",
		"result"=> "gene.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt"
	},
	"tx"=> { 
		"count"=> "tx.results.merged.count.txt",
		"selected"=> "tx.results.merged.count.selected.txt",
		"result"=> "tx.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"siresult"=> "tx.results.merged.tpm.si.Lm.FC$fccutoff.$qmethod.P$qcutoff.txt"
	}
);


########
#Process
########


#open config files to find fastqs
#----------------
#read design

my %factors;
my @factors_array;

if($formula=~/\~(.+)/) {
	my $cont=$1;
	#a bit loose on factor name
	while($cont=~/([^\+\:]+)/g) {
		
		#remove trailing white spaces
		my $factor=$1;
		$factor=~s/^\s+//;
		$factor=~s/\s+$//;
		
		$factors{$factor}++;
		push @factors_array,$factor;
	}
}

print STDERR join(",",sort keys %factors)," factors are identified from -f $formula\n\n" if $verbose;
print LOG join(",",sort keys %factors)," factors are identified from -f $formula\n\n";


#---------------
#read comparisons

my %comparisons_all;

my @trts;
my @refs;

if(defined $comparisons && length($comparisons)>0) {
	open(IN,$comparisons) || die $!;
	while(<IN>) {
		tr/\r\n//d;
		#no header
		
		my @array=split/\t/;
		
		push @trts,$array[0];
		push @refs,$array[1];
	}
}
else {
	push @trts,$treatment;
	push @refs,$reference;
}

#-------------
#sanity check for --as

if($as eq "T") {
	unless(-e "$inputfolder/".$rnaseq2files{"exon"}{"count"}) {
		print STDERR "ERROR:no exon count file found: ","$inputfolder/".$rnaseq2files{"exon"}{"count"},". You will need to run rnaseq-merge using --as T.\n";
		print LOG "ERROR:no exon count file found: ","$inputfolder/".$rnaseq2files{"exon"}{"count"},". You will need to run rnaseq-merge using --as T.\n";
		exit;
	}

	unless(-e "$inputfolder/".$rnaseq2files{"exonjunc"}{"count"}) {
		print STDERR "ERROR:no exonjunc count file found: ","$inputfolder/".$rnaseq2files{"exonjunc"}{"count"},". You will need to run rnaseq-merge using --as T.\n";
		print LOG "ERROR:no exonjunc count file found: ","$inputfolder/".$rnaseq2files{"exonjunc"}{"count"},". You will need to run rnaseq-merge using --as T.\n";
		exit;
	}

}

		
#----------------
open(S1,">$scriptfile1") || die "Error writing $scriptfile1. $!";
open(S2,">$scriptfile2") || die "Error writing $scriptfile2. $!";

for(my $compnum=0;$compnum<@trts;$compnum++) {

	my $trt=$trts[$compnum];
	my $ref=$refs[$compnum];


	my $outputfolder_de="$outputfolder/$trt\_vs\_$ref";

	if(!-d $outputfolder_de) {
		mkdir($outputfolder_de);
	}

	
	my $newconfigfile="$outputfolder_de/rnaseq-de_config.txt";
	
	#read config file
	my %sample2fastq;
	my %sample2indexname;
	my %configattrs;
	my %attr2name;

	my @configsamples;

	my $fileline=0;
	my @attrselcols;
	my @sampleselrows;

	open(IN,$configfile) || die "Error reading $configfile. $!";
	open(OUT,">$newconfigfile") || die "Error reading $newconfigfile. $!";

	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		if($fileline==0) {
			for(my $num=0;$num<@array;$num++) {
				#case insensitive match to key words
				$configattrs{uc $array[$num]}=$num;
			}
			
			foreach my $factor (sort keys %factors) {
				if(defined $configattrs{uc $factor}) {
					print STDERR "$factor is identified at the ",$configattrs{uc $factor}+1,"(th) column of $configfile.\n" if $verbose;
					print LOG "$factor is identified at the ",$configattrs{uc $factor}+1,"(th) column of $configfile.\n";
					push @attrselcols,$configattrs{uc $factor};
				}
				else {
					print STDERR "ERROR:$factor is not defined in $configfile.\n";
					print LOG "ERROR:$factor is not defined in $configfile.\n";
					exit;
				}
			}
			
			#print out new config for DE
			print OUT "Sample\t",join("\t",sort keys %factors),"\n";
		
		}
		else {
			
			push @configsamples,$array[0];
			
			#print out config file for DE
			if($useallsamples eq "T") {
				print OUT join("\t",@array[0,@attrselcols]),"\n";
			}
			else {
				#only print out used samples by last column
				if($array[$configattrs{uc $factors_array[$#factors_array]}] eq $trt || $array[$configattrs{uc $factors_array[$#factors_array]}] eq $ref) {
					print OUT join("\t",@array[0,@attrselcols]),"\n";
					push @sampleselrows,$fileline+1; #sample row number in config is the same with sample col number in count file
				}
			}
			
			foreach my $factor (sort keys %factors) {
				$attr2name{$factor}{$array[$configattrs{uc $factor}]}++;
			}
		}
		$fileline++;
	}

	close IN;
	close OUT;

	#----------------
	#test treat and ref
	if(defined $attr2name{$factors_array[$#factors_array]}{$trt}) {
		print STDERR $attr2name{$factors_array[$#factors_array]}{$trt}," samples identified for $trt in $configfile.\n" if $verbose;
		print LOG $attr2name{$factors_array[$#factors_array]}{$trt}," samples identified for $trt in $configfile.\n";
	}
	else {
		print STDERR "ERROR:$trt not defined in $configfile.\n";
		print LOG "ERROR:$trt not defined in $configfile.\n";	
		exit;
	}

	if(defined $attr2name{$factors_array[$#factors_array]}{$ref}) {
		print STDERR $attr2name{$factors_array[$#factors_array]}{$ref}," samples identified for $ref in $configfile.\n" if $verbose;
		print LOG $attr2name{$factors_array[$#factors_array]}{$ref}," samples identified for $ref in $configfile.\n";
	}
	else {
		print STDERR "ERROR:$ref not defined in $configfile.\n";
		print LOG "ERROR:$ref not defined in $configfile.\n";	
		exit;
	}


	#----------------
	#check input folder

	if(-e "$inputfolder/".$rnaseq2files{"gene"}{"count"}) {
		open(IN,"$inputfolder/".$rnaseq2files{"gene"}{"count"}) || die $!;
		while(<IN>) {
			tr/\r\n//d;
			my @array=split/\t/;
			my @mergesamples=@array[1..$#array];
			
			if(join(",",@configsamples) ne join(",",@mergesamples)) {
				print STDERR "ERROR:Sample order different.\n";
				print STDERR "ERROR:Configure file $configfile sample order:",join(",",@configsamples),"\n";
				print STDERR "ERROR:Merged file $inputfolder/",$rnaseq2files{"gene"}{"count"}," sample order:",join(",",@mergesamples),"\n";
				
				print LOG "ERROR:Sample order different.\n";
				print LOG "ERROR:Configure file $configfile sample order:",join(",",@configsamples),"\n";
				print LOG "ERROR:Merged file $inputfolder/",$rnaseq2files{"gene"}{"count"}," sample order:",join(",",@mergesamples),"\n";

				exit;
			}
			else {
				print STDERR "Sample order matched.\n" if $verbose;
				print STDERR "Configure file $configfile sample order:",join(",",@configsamples),"\n" if $verbose;
				print STDERR "Merged file $inputfolder/",$rnaseq2files{"gene"}{"count"}," sample order:",join(",",@mergesamples),"\n\n" if $verbose;
				
				print LOG "Sample order matched.\n";
				print LOG "Configure file $configfile sample order:",join(",",@configsamples),"\n";
				print LOG "Merged file $inputfolder/",$rnaseq2files{"gene"}{"count"}," sample order:",join(",",@mergesamples),"\n\n";			
			}
			last;
		}
		close IN;
	}
	else {
		print STDERR "ERROR:$inputfolder/",$rnaseq2files{"gene"}{"count"}," doesn't exist. You need to provide a rnaseq-merge folder.\n";
		print LOG "ERROR:$inputfolder/",$rnaseq2files{"gene"}{"count"}," doesn't exist. You need to provide a rnaseq-merge folder.\n";
		exit;
	}

	########
	#Filter samples
	########

	if($useallsamples eq "T") {
		#copy counting files to de folder
		print STDERR "--useallsamples T defined. All samples are used for DE test.\n\n" if $verbose;
		print LOG "--useallsamples T defined. All samples are used for DE test.\n\n";
		
		system("cp $inputfolder/".$rnaseq2files{"gene"}{"count"}." $outputfolder_de/".$rnaseq2files{"gene"}{"selected"});
		
		if($as eq "T") {
			system("cp $inputfolder/".$rnaseq2files{"tx"}{"count"}." $outputfolder_de/".$rnaseq2files{"tx"}{"selected"});
			system("cp $inputfolder/".$rnaseq2files{"tx"}{"si"}." $outputfolder_de/".$rnaseq2files{"tx"}{"siselected"});
			system("cp $inputfolder/".$rnaseq2files{"exon"}{"count"}." $outputfolder_de/".$rnaseq2files{"exon"}{"selected"});
			system("cp $inputfolder/".$rnaseq2files{"exon"}{"si"}." $outputfolder_de/".$rnaseq2files{"exon"}{"siselected"});
			system("cp $inputfolder/".$rnaseq2files{"exonjunc"}{"count"}." $outputfolder_de/".$rnaseq2files{"exonjunc"}{"selected"});
			system("cp $inputfolder/".$rnaseq2files{"exonjunc"}{"si"}." $outputfolder_de/".$rnaseq2files{"exonjunc"}{"siselected"});			
		}
	}
	else {
		#copy counting files to de folder, using selected samples
		print STDERR "--useallsamples F defined. Only selected samples are used for DE test.\n" if $verbose;
		print STDERR "Sample columns ".join(",",@sampleselrows)." are used.\n\n" if $verbose;
		
		print LOG "--useallsamples F defined. Only selected samples are used for DE test.\n";
		print LOG "Sample columns ".join(",",@sampleselrows)." are used.\n\n";

		system("cut -f 1,".join(",",@sampleselrows)." $inputfolder/".$rnaseq2files{"gene"}{"count"}." > $outputfolder_de/".$rnaseq2files{"gene"}{"selected"});
		
		if($as eq "T") {
			#for DE
			system("cut -f 1,".join(",",@sampleselrows)." $inputfolder/".$rnaseq2files{"tx"}{"count"}." > $outputfolder_de/".$rnaseq2files{"tx"}{"selected"});
			
			#for SI
			system("cut -f 1,".join(",",@sampleselrows)." $inputfolder/".$rnaseq2files{"tx"}{"si"}." > $outputfolder_de/".$rnaseq2files{"tx"}{"siselected"});	

			#for DE
			system("cut -f 1,".join(",",@sampleselrows)." $inputfolder/".$rnaseq2files{"exon"}{"count"}." > $outputfolder_de/".$rnaseq2files{"exon"}{"selected"});
			
			#for SI
			system("cut -f 1,".join(",",@sampleselrows)." $inputfolder/".$rnaseq2files{"exon"}{"si"}." > $outputfolder_de/".$rnaseq2files{"exon"}{"siselected"});				

			#for DE
			system("cut -f 1,".join(",",@sampleselrows)." $inputfolder/".$rnaseq2files{"exonjunc"}{"count"}." > $outputfolder_de/".$rnaseq2files{"exonjunc"}{"selected"});
			
			#for SI
			system("cut -f 1,".join(",",@sampleselrows)." $inputfolder/".$rnaseq2files{"exonjunc"}{"si"}." > $outputfolder_de/".$rnaseq2files{"exonjunc"}{"siselected"});				

			
		}
	}



	#Print out script
	
	#foreach my $sample (sort keys %sample2workflow) {
	#	print S1 $sample2workflow{$sample},"\n";
	#}

	#Gene #changed after v0.6
	print S1 "$descript -i $outputfolder_de/",$rnaseq2files{"gene"}{"selected"}," -c $newconfigfile -o $outputfolder_de/",$rnaseq2files{"gene"}{"result"}," -f \"$formula\" -t $trt -r $ref --fccutoff $fccutoff --qcutoff $qcutoff --qmethod $qmethod --pmethod $pmethod --filter $filter -a ",$tx2ref{$tx}{"geneanno"}," --independentfiltering $independentfiltering --cookscutoff $cookscutoff > $outputfolder_de/gene_de_test_run.log 2>&1;"; #need to check here #add r script running log

	#Gene anno
	print S1 "$mergefiles -m $outputfolder_de/",$rnaseq2files{"gene"}{"result"}," -i ".$tx2ref{$tx}{"geneanno"}." -o $outputfolder_de/",$rnaseq2files{"gene"}{"resultanno"},";\n";


	if($as eq "T") {
		#Tx DE #whether to lower TX DE cutoff?
		print S1 "$descript -i $outputfolder_de/",$rnaseq2files{"tx"}{"selected"}," -c $newconfigfile -o $outputfolder_de/",$rnaseq2files{"tx"}{"result"}," -f \"$formula\" -t $trt -r $ref --fccutoff $fccutoff --qcutoff $qcutoff --qmethod $qmethod --pmethod $pmethod --filter $filter --independentfiltering $independentfiltering --cookscutoff $cookscutoff > $outputfolder_de/tx_de_test_run.log 2>&1;";

		#tx anno
		print S1 "$mergefiles -m $outputfolder_de/",$rnaseq2files{"tx"}{"result"}," -i ".$tx2ref{$tx}{"txanno"}." -o $outputfolder_de/",$rnaseq2files{"tx"}{"resultanno"},";\n";
		
		#Tx AS
		print S1 "$descript -i $outputfolder_de/",$rnaseq2files{"tx"}{"siselected"}," -c $newconfigfile -o $outputfolder_de/",$rnaseq2files{"tx"}{"siresult"}," -f \"$formula\" -t $trt -r $ref --fccutoff $sifccutoff --qcutoff $siqcutoff --qmethod $qmethod --pmethod Lm --filter 0.001 > $outputfolder_de/tx_as_test_run.log 2>&1;";

		#tx AS anno
		print S1 "$mergefiles -m $outputfolder_de/",$rnaseq2files{"tx"}{"siresult"}," -i ".$tx2ref{$tx}{"txanno"}." -o $outputfolder_de/",$rnaseq2files{"tx"}{"siresultanno"},";\n";		
	
		#tx AS summary
		print S2 "$summarize_txsi --tx $tx -g $outputfolder_de/",$rnaseq2files{"gene"}{"result"}," --td $outputfolder_de/",$rnaseq2files{"tx"}{"result"}," --ts $outputfolder_de/",$rnaseq2files{"tx"}{"siresultanno"}," --ed $outputfolder_de/",$rnaseq2files{"exon"}{"result"}," --es $outputfolder_de/",$rnaseq2files{"exon"}{"siresultanno"}, " -o $outputfolder_de/",$rnaseq2files{"tx"}{"siresultsummary"}," > $outputfolder_de/tx_as_summary_run.log 2>&1;\n";


		#exon DE #whether to lower exon DE cutoff?
		print S1 "$descript -i $outputfolder_de/",$rnaseq2files{"exon"}{"selected"}," -c $newconfigfile -o $outputfolder_de/",$rnaseq2files{"exon"}{"result"}," -f \"$formula\" -t $trt -r $ref --fccutoff $fccutoff --qcutoff $qcutoff --qmethod $qmethod --pmethod $pmethod --filter $filter --independentfiltering $independentfiltering --cookscutoff $cookscutoff > $outputfolder_de/exon_de_test_run.log 2>&1;";

		#exon anno
		print S1 "$mergefiles -m $outputfolder_de/",$rnaseq2files{"exon"}{"result"}," -i ".$tx2ref{$tx}{"exonanno"}." -o $outputfolder_de/",$rnaseq2files{"exon"}{"resultanno"},";\n";
		
		#exon AS
		print S1 "$descript -i $outputfolder_de/",$rnaseq2files{"exon"}{"siselected"}," -c $newconfigfile -o $outputfolder_de/",$rnaseq2files{"exon"}{"siresult"}," -f \"$formula\" -t $trt -r $ref --fccutoff $sifccutoff --qcutoff $siqcutoff --qmethod $qmethod --pmethod Lm --filter 0.001 > $outputfolder_de/exon_as_test_run.log 2>&1;";

		#exon AS anno
		print S1 "$mergefiles -m $outputfolder_de/",$rnaseq2files{"exon"}{"siresult"}," -i ".$tx2ref{$tx}{"exonanno"}." -o $outputfolder_de/",$rnaseq2files{"exon"}{"siresultanno"},";\n";		
	
		#exon AS summary
		print S2 "$summarize_exonsi --tx $tx -g $outputfolder_de/",$rnaseq2files{"gene"}{"result"}," --td $outputfolder_de/",$rnaseq2files{"tx"}{"result"}," --ts $outputfolder_de/",$rnaseq2files{"tx"}{"siresultanno"}," --ed $outputfolder_de/",$rnaseq2files{"exon"}{"result"}," --es $outputfolder_de/",$rnaseq2files{"exon"}{"siresultanno"}, " -o $outputfolder_de/",$rnaseq2files{"exon"}{"siresultsummary"}," > $outputfolder_de/exon_as_summary_run.log 2>&1;\n";


		#exon junc DE
		print S1 "$descript -i $outputfolder_de/",$rnaseq2files{"exonjunc"}{"selected"}," -c $newconfigfile -o $outputfolder_de/",$rnaseq2files{"exonjunc"}{"result"}," -f \"$formula\" -t $trt -r $ref --fccutoff $fccutoff --qcutoff $qcutoff --qmethod $qmethod --pmethod $pmethod --filter $filter --independentfiltering $independentfiltering --cookscutoff $cookscutoff > $outputfolder_de/exonjunc_de_test_run.log 2>&1;";

		#exonjunc anno
		print S1 "$mergefiles -m $outputfolder_de/",$rnaseq2files{"exonjunc"}{"result"}," -i $inputfolder/".$rnaseq2files{"exonjunc"}{"anno"}." -o $outputfolder_de/",$rnaseq2files{"exonjunc"}{"resultanno"},";\n";
		
		#exonjunc AS
		print S1 "$descript -i $outputfolder_de/",$rnaseq2files{"exonjunc"}{"siselected"}," -c $newconfigfile -o $outputfolder_de/",$rnaseq2files{"exonjunc"}{"siresult"}," -f \"$formula\" -t $trt -r $ref --fccutoff $sifccutoff --qcutoff $siqcutoff --qmethod $qmethod --pmethod Lm --filter 0.001 > $outputfolder_de/exonjunc_as_test_run.log 2>&1;";

		#exonjunc AS anno
		print S1 "$mergefiles -m $outputfolder_de/",$rnaseq2files{"exonjunc"}{"siresult"}," -i $inputfolder/".$rnaseq2files{"exonjunc"}{"anno"}." -o $outputfolder_de/",$rnaseq2files{"exonjunc"}{"siresultanno"},";\n";		
	
		#exonjunc AS summary
		#print S2 "$summarize_exonjuncsi --tx $tx -g $outputfolder_de/",$rnaseq2files{"gene"}{"result"}," --td $outputfolder_de/",$rnaseq2files{"tx"}{"result"}," --ts $outputfolder_de/",$rnaseq2files{"tx"}{"siresultanno"}," --ed $outputfolder_de/",$rnaseq2files{"exonjunc"}{"result"}," --es $outputfolder_de/",$rnaseq2files{"exonjunc"}{"siresultanno"}, " -o $outputfolder_de/",$rnaseq2files{"exonjunc"}{"siresultsummary"}," > $outputfolder_de/exonjunc_as_summary_run.log 2>&1;\n";
			
	}
}

close S1;
close S2;




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



#######
#Run mode
#######

if($as eq "F") {
	submit_job($scriptfile1);
}
else {
	submit_job($scriptfile1,$scriptfile2);
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


