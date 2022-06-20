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

rnaseq-merge
version: $version
Usage: merge_exonjunc.pl [parameters]

Description: Merge exonjunc count

Parameters:

    --in|-i           Inputfile
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


my $inputfiles;
my $outputfile;
my $verbose=1;

my $dev=0; #developmental version

GetOptions(
	"in|i=s" => \$inputfiles,
	"output|o=s" => \$outputfile,

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

print STDERR "\nmerge_juncs $version running ...\n\n" if $verbose;

my %file2exonjuncs;
my $samplename;
my %exonjuncs;
my @samplenames;

foreach my $file (split(",",$inputfiles)) {
	
	if($file=~/([^\/]+)_featurecounts_exon.txt.jcounts/) {
		$samplename=$1;
		push @samplenames, $samplename;
	}
	
	open(IN,$file) || die "ERROR:Can't read $file.$!";
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		unless($_=~/PrimaryGene/) {
			#find count
			my $strand;
			if($array[4] eq "NA") {
				$strand=".";
			}
			else {
				$strand=$array[4];
			}
						
			my $exonjuncname=join("_",$array[2],$array[3],$array[6],$strand);
			
			$exonjuncs{$exonjuncname}++;
			
			$file2exonjuncs{$samplename}{$exonjuncname}=$array[8];
		}
	}
	close IN;
}
			
open(OUT,">$outputfile") || die "ERROR:Can't write to $outputfile.$!";
print OUT "Exonjunc\t",join("\t",@samplenames),"\n";
foreach my $exonjunc (sort keys %exonjuncs) {
	print OUT $exonjunc,"\t";
	my @counts;
	
	foreach my $samplename (@samplenames) {
		if(defined $file2exonjuncs{$samplename}{$exonjunc}) {
			push @counts,$file2exonjuncs{$samplename}{$exonjunc};
		}
		else {
			push @counts,0;
		}
	}
	
	print OUT join("\t",@counts),"\n";
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

