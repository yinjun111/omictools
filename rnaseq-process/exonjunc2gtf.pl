#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);

#CutAdapt+FASTQC+RSEM+STAR


########
#Interface
########


my $version="0.1";


my $usage="

exonjunc2gtf
version: $version
Usage: exonjunc2gtf.pl [parameters]

Description: Annotate the exon juncs

Parameters:

    --in|-i           Inputfile
    --anno|-a         Annotation file
    --output|-o       Outputfile


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


my $inputfile;
my $outputfile;
my $annofile;

my $verbose=1;

my $dev=0; #developmental version

GetOptions(
	"in|i=s" => \$inputfile,
	"output|o=s" => \$outputfile,
	"anno|a=s" => \$annofile,
	
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


########
#default ouputs
########




########
#Process
########


print STDERR "\nexonjunc2gtf $version running ...\n\n" if $verbose;


#Read exon junc anno	
print STDERR "Reading exon junc annotation file:$annofile\n\n";

my %exonjunc2name;

open(IN,$annofile) || die $!;
while(<IN>) {
	tr/\r\n//d;

	my @array=split/\t/;
	
	$exonjunc2name{$array[0]}=$array[1];
		
}
close IN;


print STDERR "Processing exon junction file:$inputfile\n\n";

#generate annotation for exons
open(IN,$inputfile) || die $!;
open(OUT,">$outputfile") || die $!;

my $samplename;
if($inputfile=~/([^\/]+)_featurecounts_exon.txt.jcounts/) {
	$samplename=$1;
}
else {
	$samplename=basename($inputfile);
}

#print GTF header
print OUT "track name=$samplename\_exon_junctions\n";


while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Exon|PrimaryGene/;
	
	my @array=split/\t/;

	#GTF example
	#<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
	#AB000381 Twinscan  CDS          380   401   .   +   0  gene_id "001"; transcript_id "001.1";

	#decide exon junc location
	my ($ejchr,$ejstart,$ejend,$ejstr)=@array[2,3,6,4];

	if($ejstr eq "NA") {
		$ejstr=".";
	}
	
	my $exonjuncid=join("_",$ejchr,$ejstart,$ejend,$ejstr);
	
	my ($newstart,$newend);
	
	#first part
	print OUT "$ejchr\tfeaturecounts\texon\t";
	if($ejstart>=10) {
		print OUT $ejstart-9,"\t$ejstart\t$array[8]\t$ejstr\t0\t";
		$newstart=$ejstart-9;
	}
	else {
		print OUT "1\t$ejstart\t$array[8]\t$ejstr\t0\t";
		$newstart=1;
	}
	
	#name of exon junc
	if(defined $exonjunc2name{$exonjuncid}) {
		print OUT "transcript_id \"",$exonjunc2name{$exonjuncid},"\";\n";
	}
	else {
		print OUT "transcript_id \"$exonjuncid\";\n";
	}
		
	#second part
	print OUT "$ejchr\tfeaturecounts\texon\t";
	print OUT "$ejend\t",$ejend+9,"\t$array[8]\t$ejstr\t0\t";
	$newend=$ejend+9;
	
	#name of exon junc
	if(defined $exonjunc2name{$exonjuncid}) {
		print OUT "transcript_id \"",$exonjunc2name{$exonjuncid},"\";\n";
	}
	else {
		print OUT "transcript_id \"$exonjuncid\";\n";
	}

}
close OUT;

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

