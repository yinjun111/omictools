#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);
use Text::CSV;


########
#Version
########

my $version="0.2";

#v0.2, add arguments

########
#Interface
########

my $usage="

check_download_results
version: $version
Usage: perl check_download_results.pl --in infolder/ --script sra_download.sh --seq -o sra_download_1.sh

Description: Check SRA download results and make new script for lost files

Parameters:

    --in             Infolder containing SRA files
    --script         sra_download.sh
    --out|-o         output file, sra_download_1.sh

	
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


my $infolder;
my $inscript;
my $outscript;
my $verbose=1;

GetOptions(

	"in|i=s" => \$infolder,
	"script|s=s" => \$inscript,
	"out|o=s" => \$outscript,

	"verbose|v" => \$verbose,
);


#check fastq files by name, and match to download list, generate new sh for files not downloaded
#my ($infolder,$inscript, $outscript)=@ARGV;

my @fastqfiles=glob("$infolder/*.fastq.gz");

my %srrs;
foreach my $fastqfile (@fastqfiles) {
	if($fastqfile=~/(SRR\d+)/) {
		$srrs{$1}++;
	}
}

open(IN,$inscript) || die $!;
open(OUT,">$outscript") || die $!;
while(<IN>) {
	tr/\r\n//d;
	if($_=~/(SRR\d+)/) {
		unless(defined $srrs{$1}) {
			print STDERR $1," is not found in --in ",abs_path($infolder),".\n";
			print OUT $_,"\n";
		}
	}
}
close IN;
close OUT;
