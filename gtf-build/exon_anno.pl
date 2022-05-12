#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);


#version
my $version=0.1;

#0.1, exon annotaiton based on GTF


my $usage="

exon_anno
version: $version
Usage: perl exon_anno.pl -i ucsc.gtf -o Homo_sapiens.GRCh38.88_ucsc_exon_annocombo.txt

Description: 

Parameters:

    --in|-i           input file
    --out|-o          out file

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
my $verbose;

GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfile,
	"verbose|v" => \$verbose,
);





########
#Program starts
########

#GTF file
open(IN,$infile) || die $!;

my %exon2gene;
my %exon2tx;

my %gene2name;
my %tx2name;
my %exon2name;

my %exon2info;

#record, exonid, tx, gene, chr, start, end, direction
while(<IN>) {
	next if $_=~/^#/;
	tr/\r\n//d;
	my @array=split/\t/;

	if($array[2] eq "exon") {
		#exon info
		my %attrs;
		
		while($array[8]=~/(\w+) \"([^\"]+)\"/g) {
			$attrs{$1}=$2;
		}
	
		$gene2name{$attrs{"gene_id"}}=$attrs{"gene_name"};
		$tx2name{$attrs{"transcript_id"}}=$attrs{"transcript_name"};		
		#$exon2name{$attrs{"exon_id"}}{$attrs{"transcript_name"}.".E".$attrs{"exon_number"}}++;
		
		$exon2gene{$attrs{"exon_id"}}{$attrs{"gene_id"}}++;
		$exon2tx{$attrs{"exon_id"}}{$attrs{"transcript_id"}}=$attrs{"transcript_name"}.".E".$attrs{"exon_number"};
		
		$exon2info{$attrs{"exon_id"}}=join("\t",@array[0,3,4,6],$array[4]-$array[3]+1);
	}
}
close IN;



#output
open(OUT,">$outfile") || die $!;
print OUT "Exon\texon_name\tGene\tgene_name\tTranscript\ttranscript_name\tChr\tStart\tEnd\tStrand\tExonLength\n";

foreach my $exon (sort keys %exon2gene) {
	print OUT $exon,"\t";
	
	#exon name #use the same order by transcript
	my @exonnames;
	foreach my $tx (sort keys %{$exon2tx{$exon}}) {
		push @exonnames,$exon2tx{$exon}{$tx};
	}
	print OUT join(",",@exonnames),"\t";
	
	#gene
	print OUT join(",",sort keys %{$exon2gene{$exon}}),"\t \t";
	my @genenames;
	foreach my $gene (sort keys %{$exon2gene{$exon}}) {
		push @genenames,$gene2name{$gene};
	}
	print OUT join(",",@genenames),"\t";
	
	#tx
	print OUT join(",",sort keys %{$exon2tx{$exon}}),"\t \t";
	my @txnames;
	foreach my $tx (sort keys %{$exon2tx{$exon}}) {
		push @txnames,$tx2name{$tx};
	}
	print OUT join(",",@txnames),"\t";	
	
	#anno
	print OUT $exon2info{$exon},"\n";
}

close OUT;

