#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

########
#Updates
########



########
#Prerequisites
########

#none


########
#Interface
########


my $version="0.31";

#0.2, add FILTER tag
#0.21, change dp to <, not <=. Same as VariantFiltration
#0.3, add AF filter
#0.31, add log file

my $usage="

rnaseq-var_filter
version: $version
Usage: 

Description: 

Parameters:
    --in|-i           Input file, a vcf file
    --out|-o          Output file, a filtered vcf file
    --com|-c          Common SNP vcf
    --dp|-d           DP filter, recommended 5
    --af|-a           AF filter, recommended 0.05
    --filter|-f       Only retain SNPs by PASS tag [PASS]
	
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

my $infile;
my $outfile;
my $commonsnp;
my $dp;
my $af;
my $filter="PASS";

my $runmode=0;
my $verbose=1;


GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfile,
	"com|c=s" => \$commonsnp,
	"dp|d=s" => \$dp,
	"af|a=s" => \$af,
	"filter|f=s" => \$filter,

	"verbose" => \$verbose,
	"help|h" => sub {print STDERR $usage;exit;}
	
);

my $outfilefolder=abs_path_dir(abs_path($outfile));

my $logfile="$outfilefolder/rnaseq-var_filter_run.log";

my %filters=map {$_,1} split(",",$filter);


######
#Process input file
######

print STDERR "\nrnaseq-var_filter version $version running...\n";

#Read common snp inforamtion
my %commonsnps;
if(defined $commonsnp && length($commonsnp)>0) {
	open(IN,$commonsnp) || die $!;
	while(<IN>) {
		tr/\r\n//d;
		next if ($_=~/^\#/);

		my @array=split/\t/;
		
		my $snpname=join(":",@array[0,1,3,4]); #index by chr:coord:ref:alt

		$commonsnps{$snpname}++;
	}
	close IN;
}


open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;

my $anum=0;#all snps
#my $cnum=0;#common snps
#my $dnum=0;#dp snps
#my $fnum=0;#filtered snps

my %filternums;

my $finalnum=0; #all the filtered snps

while(<IN>) {
	tr/\r\n//d;
	if ($_=~/^\#/) {
		print OUT $_,"\n";
	}
	else {	
		$anum++;
		
		my @array=split/\t/;
		
		my $snpname=join(":",@array[0,1,3,4]); #index by chr:coord:ref:alt
		
		#filter common snp
		if(defined $commonsnps{$snpname}) {
			#$cnum++;			
			$filternums{"CommonSNP"}++;		
			next;
		}
		
		
		#filter existing DP
		if( (!defined $dp || length($dp)==0) && $array[6] =~/DP/) {
			#$dnum++;
			$filternums{"DP"}++;		
			next;
		}
		
		#filter by FILTER column
		unless(defined $filters{$array[6]}) {
			#$fnum++;
			$filternums{"GATK-Filter:".$array[6]}++;		
			next;
		}
		
		
		#use new dp filter
		if(defined $dp && length($dp)>0) {
			if($array[7]=~/DP=(\d+)/) {
				if($1<$dp) {
					#0.21, change dp to <, not <=. Same as VariantFiltration
					#$dnum++;
					$filternums{"DP<".$dp}++;		
					next;
				}
			}
		}

		#use af filter
		my @cates=split(":",$array[8]);
		my @vals=split(":",$array[9]);
		my %cate2val;
		
		for(my $num=0;$num<@cates;$num++) {
			$cate2val{$cates[$num]}=$vals[$num];
		}
		
		if(defined $af && length($af)>0) {
			if($cate2val{"AF"}<$af) {
				$filternums{"AF<".$af}++;		
				next;
			}
		}
		
		print OUT $_,"\n";
		
		$finalnum++;
		
	}
}
close IN;
close OUT;

#print STDERR "$anum SNPs are processed using rnaseq-var_filter.\n$cnum SNPs are filtered by matching common SNPs.\n$dnum SNPs are filtered by matching DP filter.\n$fnum SNPs are filtered by matching FILTER $filter tag.\n$finalnum SNPs passed all the filters.\n\n";

open(LOG,">$logfile") || die $!;

print LOG "$anum SNPs are processed using rnaseq-var_filter.\n";
print LOG "Variants filtered:\n";
print LOG join("\n",map {$_."\t".$filternums{$_}} sort keys %filternums),"\n\n";
print LOG "$finalnum SNPs passed all the filters.\n\n";

close LOG;



######
#Functions
######


sub abs_path_dir {
	my $path=shift @_;
	my $path_dir;
	if($path=~/[^\/]+$/) {
		$path_dir=$`;
	}
	return $path_dir;
}