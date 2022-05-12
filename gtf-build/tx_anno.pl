#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);


#version
my $version=0.2;

#0.2, Use gene symbol as index. Add .number for duplicated gene names. If duplicated, use the longest as the default one


my $usage="

gtf2anno
version: $version
Usage: perl gtf2anno.pl -i genome.fa -o genome_ucsc.fa

Description: Produce gene annotation used by omictools. Rename duplicated genes by giving .number for duplicated gene names. If duplicated, use the longest as the default one.

Parameters:

    --in|-i           input file
    --out|-o          out file
    --anno|-a         annotation file

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

my $infile;
my $outfile;
my $annofile;
my $verbose;

GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfile,
	"anno|a=s" => \$annofile,
	"verbose|v" => \$verbose,
);


my %tx2anno;
my %tx2name;
my $txtitle;
my $txannocol;

open(IN,$annofile) || die $!;

while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;

	if ($_=~/^Transcript/) {
		$txtitle=join("\t",@array[1..$#array]);
		$txannocol=@array-2;
	}
	else {
		$tx2name{$array[0]}=$array[1]; #name
		$tx2anno{$array[0]}=join("\t",@array[2..$#array]); #others
		
	}
}
close IN;
close OUT;




my $txcol=5;

#from existing gene anno, build transcript anno

open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;

my $linenum=0;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;

	if($linenum==0) {
		print OUT "Transcript\ttranscript_name\t",join("\t",@array[0..($txcol-2),$txcol..$#array],$txtitle),"\n";
	}
	else {
		foreach my $tx (split(",",$array[$txcol-1])) {			
			if(defined $tx2anno{$tx}) {
				print OUT $tx,"\t",$tx2name{$tx},"\t",join("\t",@array[0..($txcol-2),$txcol..$#array]),"\t";
				print OUT $tx2anno{$tx},"\n";
			}
			else {
				print OUT $tx,"\t \t",join("\t",@array[0..($txcol-2),$txcol..$#array]),"\t";
				print OUT join("\t",(" ") x $txannocol),"\n";
			}
		}
	}
	$linenum++;
}
close IN;
close OUT;
