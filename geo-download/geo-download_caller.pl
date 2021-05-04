#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);
use Text::CSV;

#CutAdapt+FASTQC+RSEM+STAR



########
#Interface
########


my $version="0.3";

#v0.2, Fastq files for PE
#v0.3, change to fasterq-dump

my $usage="

geo-download
version: $version
Usage: omictools geo-download --sra SraRunTable.txt --geo GSE78220_series_matrix.txt --seq PE -o ./

Description: Download GEO SRA fastq files.

Parameters:

    --sra             SRA file, SraRunTable.txt (csv). SRA should be manually checked and corrected for annotations.
    --geo             GEO series matrix file(s), separate by \",\"
    --seq             SE or PE [SE]
    --out|-o          output folder

	
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


my $srafile;
my $geofile;
my $outputfolder;
my $seq="SE";
my $verbose=1;

GetOptions(

	"sra=s" => \$srafile,
	"geo=s" => \$geofile,
	"seq=s" => \$seq,
	"out|o=s" => \$outputfolder,

	"verbose|v" => \$verbose,
);


########
#Prerequisites
########

my $fasterq_dump="/apps/sratoolkit.2.11.0-centos_linux64/bin/fasterq-dump";
my $fastq_dump="/apps/sratoolkit.2.11.0-centos_linux64/bin/fastq-dump";


########
#I/O
########

$outputfolder=abs_path($outputfolder);

unless(-e $outputfolder) {
	mkdir $outputfolder;
}

my $logfile="$outputfolder/geo-download_run.log";

#write log file
open(LOG, ">$logfile") || die "Error write $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";


########
#Process
########

print STDERR "\nomictools geo-download $version\n\n";
print LOG "\nomictools geo-download $version running ...\n\n";


#SRA file, read GSM ID
my %sra2gsm;
my %sra2info;
my $sratitle;
my %title2col;

my %gsm2sra; #gsm may have multiple sras

open(IN,$srafile) || die "ERROR:Can't read $srafile.\n\n";

my $csv=Text::CSV->new({ sep_char => ',' });

while(<IN>) {
	tr/\r\n//d;
	
	if($csv->parse($_)) {
		my @array=$csv->fields();
	
	#my @array=split/\t/;
	
		if($_=~/^Run/) {
				
			for(my $num=0;$num<@array;$num++) {
				$title2col{$array[$num]}=$num;
			}
			
			$sratitle=$_;
			
		}
		else {
			$sra2gsm{$array[0]}=$array[$title2col{"GEO_Accession (exp)"}];
			$sra2info{$array[0]}=[@array];
			$gsm2sra{$array[$title2col{"GEO_Accession (exp)"}]}{$array[0]}++; #One GSM ID may contain multiple SRAs
		}
	}
	else {
		warn "Line could not be parsed: $_\n";
	}
}
close IN;

print STDERR scalar(keys %sra2gsm)," SRA IDs were read.\n";
print STDERR join(",",map {$_."-".$sra2gsm{$_}} sort keys %sra2gsm),"\n\n";

print LOG scalar(keys %sra2gsm)," SRA IDs were read.\n";
print LOG join(",",map {$_."-".$sra2gsm{$_}} sort keys %sra2gsm),"\n\n";

#GEO Series, read sample title

my %gsm2title;

my @geotitles;
my @geogsms;

foreach my $gfile (split(",",$geofile)) {
	open(IN,$gfile) || die "ERROR:Can't read $gfile.\n\n";

	while(<IN>) {
		tr/\r\n\"//d;
		my @array=split/\t/;
		
		if($_=~/^\!Sample_title/) {
			push @geotitles, @array[1..$#array];	
		}

		if($_=~/^\!Sample_geo_accession/) {
			push @geogsms, @array[1..$#array];	
		}
	}
	close IN;
}


my %knowngsms;
for(my $num=0;$num<@geotitles;$num++) {
	$gsm2title{$geogsms[$num]}=$geotitles[$num];
	if(defined $gsm2sra{$geogsms[$num]}) {
		$knowngsms{$geogsms[$num]}++;
	}
}

print STDERR scalar(@geotitles)," GEO IDs were read.",scalar(keys %knowngsms)," matched to SRA file.\n";
print LOG scalar(@geotitles)," GEO IDs were read.",scalar(keys %knowngsms)," matched to SRA file.\n";

print STDERR join(",",@geogsms),"\n\n";
print LOG join(",",@geogsms),"\n\n";

#Write out

my $foldername=basename($outputfolder);

open(OUT,">$outputfolder/sra_download.sh") || die $!;

foreach my $sra (sort keys %sra2gsm) {
	if(defined $gsm2title{$sra2gsm{$sra}}) {
		#print OUT "fastq-dump --gzip -F --split-3 $sra -O $outputfolder\n";
		print OUT "$fasterq_dump $sra -O $outputfolder;gzip $outputfolder/$sra*;\n";	
	}
	else {
		print STDERR "ERROR:$sra2gsm{$sra} not defined.\n";
	}
}

close OUT;

open(OUT,">$outputfolder/sra_merge.sh") || die $!;

foreach my $gsm (sort keys %gsm2sra) {
	my @srafiles;
	my @srafiles_rev;
	
	foreach my $sra (sort keys %{$gsm2sra{$gsm}}) {
		if($seq eq "SE") {
			push @srafiles,"$outputfolder/$sra.fastq.gz";
		}
		else {
			#for PE, the files are SRR3184279_1.fastq.gz, SRR3184279_2.fastq.gz
			push @srafiles,"$outputfolder/$sra\_1.fastq.gz";
			push @srafiles_rev,"$outputfolder/$sra\_2.fastq.gz";
			
		}
	}
	
	if($seq eq "SE") {
		if(@srafiles==1) {
			print OUT "mv ",join(" ",@srafiles)," $outputfolder/$gsm\_merged.fastq.gz\n";			
		} else {
			print OUT "cat ",join(" ",@srafiles)," > $outputfolder/$gsm\_merged.fastq.gz\n";
		}
	}
	else {
		if(@srafiles==1) {
			print OUT "mv ",join(" ",@srafiles)," $outputfolder/$gsm\_merged_1.fastq.gz\n";
			print OUT "mv ",join(" ",@srafiles_rev)," $outputfolder/$gsm\_merged_2.fastq.gz\n";
		}
		else {	
			print OUT "cat ",join(" ",@srafiles)," > $outputfolder/$gsm\_merged_1.fastq.gz\n";
			print OUT "cat ",join(" ",@srafiles_rev)," > $outputfolder/$gsm\_merged_2.fastq.gz\n";
		}
	}
}
close OUT;


open(ANNO,">$outputfolder/$foldername\_samples_anno.txt") || die $!;

if($seq eq "SE") {
	print ANNO "Sample\tFastq1\tGSM\tSRAs\t$sratitle\n"; #now SE RNASeq ...

	foreach my $gsm (sort keys %gsm2sra) {
		print ANNO $gsm2title{$gsm},"\t$outputfolder/$gsm\_merged.fastq.gz\t$gsm\t";
		print ANNO join(",",sort keys %{$gsm2sra{$gsm}}),"\t";
		
		#print the first SRA info
		print ANNO join("\t",@{$sra2info{(sort keys %{$gsm2sra{$gsm}})[0]}}),"\n";
	}

	close ANNO;
}
else {
	print ANNO "Sample\tFastq1\tFastq2\tGSM\tSRAs\t$sratitle\n"; #now SE RNASeq ...

	foreach my $gsm (sort keys %gsm2sra) {
		print ANNO $gsm2title{$gsm},"\t$outputfolder/$gsm\_merged_1.fastq.gz\t$outputfolder/$gsm\_merged_2.fastq.gz\t$gsm\t";
		print ANNO join(",",sort keys %{$gsm2sra{$gsm}}),"\t";
		
		#print the first SRA info
		print ANNO join("\t",@{$sra2info{(sort keys %{$gsm2sra{$gsm}})[0]}}),"\n";
	}

	close ANNO;

}


#####
#Functions
######

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}
