#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);


#version
my $version=0.21;

#0.2, Use gene symbol as index. Add .number for duplicated gene names. If duplicated, use the longest as the default one
#0.21, Dont' take genes not in complete chromosomes

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
my $type;
my $annofile;
my $verbose;

GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfile,
	"type|t=s" => \$type,
	"anno|a=s" => \$annofile,
	"verbose|v" => \$verbose,
);


my $logfile=$outfile;
$logfile=~s/\.\w+$/_ensembl2ucsc.log/;


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "\n";



########
#Process
########

print STDERR "\nomictools gtf2anno $version\n\n";
print LOG "\nomictools gtf2anno $version\n\n";


#attr used from gtf
my @attrs=qw(gene_name description gene_biotype transcript_id strand chromosome start end gc_content);

my %gene2info;

#read file
open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^#/;
	my @array=split/\t/;
	if($array[2] eq "gene") {
		my %terms=process_info($array[8]);
		
		my $geneid=$terms{"gene_id"};
		
		foreach  my $attr (sort keys %terms) {
			$gene2info{$geneid}{$attr}{$terms{$attr}}++;
		}
		
		#str
		$gene2info{$geneid}{"strand"}=$array[6];
		#chr
		$gene2info{$geneid}{"chromosome"}=$array[0];
		#start, end
		$gene2info{$geneid}{"start"}=$array[3];
		$gene2info{$geneid}{"end"}=$array[4];
	}
	
	
	if($array[2] eq "transcript") {
		my %terms=process_info($array[8]);
		
		my $geneid=$terms{"gene_id"};
		
		foreach  my $attr (sort keys %terms) {
			$gene2info{$geneid}{$attr}{$terms{$attr}}++;
		}
	}
}

#annofile

my $linenum=0;
my @annotitle; #include Gene Description and %GC
open(IN,$annofile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^#/;
	my @array=split/\t/;

	if($linenum==0) {
		@annotitle=@array;
	}
	else {
		if(defined $gene2info{$array[0]}) {
			#don't need genes not defined above, which are genes not in complete chromosomes
			$gene2info{$array[0]}{"description"}=$array[1];
			$gene2info{$array[0]}{"gc_content"}=$array[2];
		}
	}
	
	$linenum++;
}



open(OUT,">$outfile") || die $!;
print OUT "Gene\t",join("\t",@attrs),"\n";
foreach my $geneid (sort keys %gene2info) {
	print OUT $geneid,"\t";
	my @contents;
	
	foreach my $attr (@attrs) {
		if(defined $gene2info{$geneid}{$attr}) {
			if(ref($gene2info{$geneid}{$attr}) eq "HASH") {
				push @contents,join(",",sort keys %{$gene2info{$geneid}{$attr}});
			}
			else {
				push @contents,$gene2info{$geneid}{$attr};
			}
		}
		else {
			push @contents," ";
		}
	}
	
	print OUT join("\t",@contents),"\n";
}

close OUT;
		


########
#Functions
########

sub process_info {
	my $info=shift @_;
	
	my %infos;
	
	while($info=~/(\w+) ([^;]+);/g) {
		my $attr=$1;
		my $value=$2;
		$value=~tr/"//d;
		$infos{$attr}=$value;
	}
	
	return %infos;
}


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
	
	#use defined program as default, otherwise search for this program in PATH
	
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

