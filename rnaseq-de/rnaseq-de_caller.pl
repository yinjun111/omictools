#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);


########
#Prerequisites
########


my $rscript="/apps/R-3.4.1/bin/Rscript";

my $descript="/home/jyin/Projects/Pipeline/sbptools/rnaseq-de/de_test_caller.R";


########
#Interface
########


my $version="0.1";

my $usage="

rnaseq-de
version: $version
Usage: sbptools rnaseq-de [parameters]

Description: Differential Expression (DE) tests using DESeq2. This script works for most of the counting based data, e.g. RNA-Seq, ChIP-Seq, ATAC-Seq


Mandatory Parameters:
    --in|-i           Input folder from rnaseq-merge
    --output|-o       Output folder
    --config|-c       Configuration file match the samples in the rnaseq-merge folder
                           first column as sample name.

    --formula|-f      Formula for GLM, e.g. ~Group.
                          the last factor of the formula is used for comparison
    --treatment       Treatment group name
    --reference       Reference group name

    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84
						
Optional Parameters:
    --pmethod         DESeq2 method, default as Wald test [Wald]
    --qmethod         Multiple testing correction method [BH]

    --useallsamples   Use all samples in config file for 
                           gene dispersion calcultation [F]

    --filter          Signal filter [auto]
                         automatically defined signal cutoff as
                           Count >= No. of samples * 2
                         or can define a count number

    --fccutoff        Log2 FC cutoff [1]
    --qcutoff         Corrected P cutoff [0.05]

    --runmode|-r      Where to run the scripts, local, server or none [none]
    --verbose|-v      Verbose
	
";


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
	
my $inputfolder;	
my $configfile;
my $outputfolder;

my $formula;
my $treatment;
my $reference;

my $pmethod="Wald";
my $qmethod="BH";
my $useallsamples="F";
my $filter="auto";

my $fccutoff=1;
my $qcutoff=0.05;

my $verbose;
my $tx;
my $runmode="none";

GetOptions(
	"in|i=s" => \$inputfolder,
	"config|c=s" => \$configfile,
	"out|o=s" => \$outputfolder,
	"formula|f=s" => \$formula,
	
	"treatment=s" => \$treatment,
	"reference=s" => \$reference,
	
	"pmethod=s" => \$pmethod,
	"qmethod=s" => \$qmethod,
	"useallsamples=s" => \$useallsamples,
	"filter=s" => \$filter,
	
	"fccutoff=s" => \$fccutoff,
	"qcutoff=s" => \$qcutoff,
	
	"tx|t=s" => \$tx,	
	"runmode|r=s" => \$runmode,		
	"verbose|v" => \$verbose,
);

$inputfolder = abs_path($inputfolder);
$outputfolder = abs_path($outputfolder);

if(!-e $outputfolder) {
	mkdir($outputfolder);
}


my $scriptfolder="$outputfolder/scripts";

if(!-e $scriptfolder) {
	mkdir($scriptfolder);
}


my $logfile="$outputfolder/rnaseq-de_run.log";
my $newconfigfile="$outputfolder/rnaseq-de_config.txt";

my $scriptfile1="$scriptfolder/rnaseq-de_run.sh";

#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

#Report R package version here !!!



print LOG "\n";

print STDERR "\nsbptools rnaseq-de $version running ...\n\n" if $verbose;
print LOG "\nsbptools rnaseq-de $version running ...\n\n";


#test tx option

my %tx2ref=(
	"Human.B38.Ensembl84"=> { 
		"star"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR",
		"chrsize"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Human/hg38/",
		"txanno"=>"/data/jyin/Databases/Genomes/Human/hg38/"},
	"Mouse.B38.Ensembl84"=>{ 
		"star"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR",
		"chrsize"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/",
		"txanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/"}
);


if(defined $tx2ref{$tx}) {
	print STDERR "Starting analysis using $tx.\n\n" if $verbose;
	print LOG "Starting analysis using $tx.\n\n";
}
else {
	print STDERR "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n" if $verbose;
	print LOG "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n";
}


#RNA-Seq
my $genecountmerged="gene.results.merged.count.txt";
my $txcountmerged="tx.results.merged.count.txt";

my $genederesult="gene.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt";
my $txderesult="tx.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt";


#ChIP/ATAC-Seq



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
	while($cont=~/([^\+]+)/g) {
		
		#remove trailing white spaces
		my $factor=$1;
		$factor=~s/^\s+//;
		$factor=~s/\s+$//;
		
		$factors{$factor}++;
		push @factors_array,$factor;
	}
}

print STDERR join(",",@factors_array)," factors are identified from -f $formula\n\n" if $verbose;
print LOG join(",",@factors_array)," factors are identified from -f $formula\n\n";
		
#----------------
#read config file

my %sample2fastq;
my %sample2indexname;
my %configattrs;
my %attr2name;

my @configsamples;

my $fileline=0;
my @attrselcols;

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
		
		foreach my $factor (@factors_array) {
			if(defined $configattrs{uc $factor}) {
				print STDERR "$factor is identified at the ",$configattrs{uc $factor},"(th) column of $configfile.\n" if $verbose;
				print LOG "$factor is identified at the ",$configattrs{uc $factor},"(th) column of $configfile.\n";
				push @attrselcols,$configattrs{uc $factor};
			}
			else {
				print STDERR "ERROR:$factor is not defined in $configfile.\n";
				print LOG "ERROR:$factor is not defined in $configfile.\n";
				exit;
			}
		}
		
		#print out new config for DE
		print OUT "Sample\t",join("\t",@factors_array),"\n";
	
	}
	else {
		
		push @configsamples,$array[0];
		
		print OUT join("\t",@array[0,@attrselcols]),"\n";
		
		foreach my $factor (@factors_array) {
			$attr2name{$factor}{$array[$configattrs{uc $factor}]}++;
		}
	}
	$fileline++;
}

close IN;
close OUT;

#----------------
#test treat and ref
if(defined $attr2name{$factors_array[$#factors_array]}{$treatment}) {
	print STDERR $attr2name{$factors_array[$#factors_array]}{$treatment}," samples identified for $treatment in $configfile.\n" if $verbose;
	print LOG $attr2name{$factors_array[$#factors_array]}{$treatment}," samples identified for $treatment in $configfile.\n";
}
else {
	print STDERR "ERROR:$treatment not defined in $configfile.\n";
	print LOG "ERROR:$treatment not defined in $configfile.\n";	
	exit;
}

if(defined $attr2name{$factors_array[$#factors_array]}{$reference}) {
	print STDERR $attr2name{$factors_array[$#factors_array]}{$reference}," samples identified for $reference in $configfile.\n" if $verbose;
	print LOG $attr2name{$factors_array[$#factors_array]}{$reference}," samples identified for $reference in $configfile.\n";
}
else {
	print STDERR "ERROR:$reference not defined in $configfile.\n";
	print LOG "ERROR:$reference not defined in $configfile.\n";	
	exit;
}


#----------------
#check input folder

if(-e "$inputfolder/$genecountmerged") {
	open(IN,"$inputfolder/$genecountmerged") || die $!;
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		my @mergesamples=@array[1..$#array];
		
		if(join(",",@configsamples) ne join(",",@mergesamples)) {
			print STDERR "ERROR:Sample order different.\n";
			print STDERR "ERROR:Configure file $configfile sample order:",join(",",@configsamples),"\n";
			print STDERR "ERROR:Merged file $inputfolder/$genecountmerged sample order:",join(",",@mergesamples),"\n";
			
			print LOG "ERROR:Sample order different.\n";
			print LOG "ERROR:Configure file $configfile sample order:",join(",",@configsamples),"\n";
			print LOG "ERROR:Merged file $inputfolder/$genecountmerged sample order:",join(",",@mergesamples),"\n";

			exit;
		}
		else {
			print STDERR "Sample order matched.\n" if $verbose;
			print STDERR "Configure file $configfile sample order:",join(",",@configsamples),"\n" if $verbose;
			print STDERR "Merged file $inputfolder/$genecountmerged sample order:",join(",",@mergesamples),"\n" if $verbose;
			
			print LOG "Sample order matched.\n";
			print LOG "Configure file $configfile sample order:",join(",",@configsamples),"\n";
			print LOG "Merged file $inputfolder/$genecountmerged sample order:",join(",",@mergesamples),"\n";			
		}
		last;
	}
	close IN;
}
else {
	print STDERR "ERROR:$inputfolder/$genecountmerged doesn't exist. You need to provide a rnaseq-merge folder.\n";
	print LOG "ERROR:$inputfolder/$genecountmerged doesn't exist. You need to provide a rnaseq-merge folder.\n";
	exit;
}

########
#Filter samples
########

if($useallsamples eq "T") {
	#do something

}

else {
	#do something


}



########
#Print out commands, for local and server run
########


open(S1,">$scriptfile1") || die "Error writing $scriptfile1. $!";

#foreach my $sample (sort keys %sample2workflow) {
#	print S1 $sample2workflow{$sample},"\n";
#}

#Gene
print S1 "$rscript $descript -i $inputfolder/$genecountmerged -a $newconfigfile -o $outputfolder/$genederesult -f \"$formula\" -t $treatment -r $reference --fccutoff $fccutoff --qcutoff $qcutoff --qmethod $qmethod --pmethod $pmethod --filter $filter;\n";

#add annotation here!

#Gene
print S1 "$rscript $descript -i $inputfolder/$txcountmerged -a $newconfigfile -o $outputfolder/$txderesult -f \"$formula\" -t $treatment -r $reference --fccutoff $fccutoff --qcutoff $qcutoff --qmethod $qmethod --pmethod $pmethod --filter $filter;\n";

#add annotation here!

close S1;


#local mode
print STDERR "To run locally, in shell type: sh $scriptfile1\n\n";
print LOG "To run locally, in shell type: sh $scriptfile1\n\n";


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












