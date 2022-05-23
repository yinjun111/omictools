#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);



########
#Interface
########


my $version="0.1";

#0.1 summarize exonde, exonsi, txde,txsi,genede



my $usage="

summarize_exonsi
version: $version
Usage: perl summarize_exonsi [parameters]

Description: summarize exon alternative splicing results


Mandatory Parameters:
    --tx|-t           Transcriptome
                        Currently support Human.B38.Ensembl88,Mouse.B38.Ensembl88,Rat.Rn6.Ensembl88

    -g|-g             genede
    --td              txde
    --ts              txsi
    --ed              exonde
    --es              exonsi
	
    --output|-o       Output file 

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

my $tx;	
my $genede;	
my $txde;	
my $txsi;	
my $exonde;
my $exonsi;
my $outputfile;	
my $dev=0;

GetOptions(
	"tx|t=s" => \$tx,
	"g=s" => \$genede,
	"td=s" => \$txde,
	"ts=s" => \$txsi,
	"ed=s" => \$exonde,
	"es=s" => \$exonsi,	
	"o=s" => \$outputfile,
	
	"dev" => \$dev,			
);

########
#Prerequisites
########

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


#######
#Input/Output
#######


##test tx option
#may need to change for different annotation versions

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
	print STDERR "Starting analysis using $tx.\n\n";
}
else {
	print STDERR "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n";
	exit;
}



########
#Process
########


my %exon2gene;
my %exon2tx;
my %exon2anno;
my $annotitle;
open(IN,$tx2ref{$tx}{"exonanno"}) || die $!;
while(<IN>) {

	tr/\r\n//d;
	my @array=split/\t/;

	if ($_=~/^Exon/) {
		$annotitle=join("\t",@array[1..$#array]);
	}
	else {
		$exon2gene{$array[0]}=$array[2];
		
		foreach my $transcript (split(",",$array[4])) {
			push @{$exon2tx{$array[0]}},$transcript;
		}
		
		$exon2anno{$array[0]}=join("\t",@array[1..$#array]);
	}
}
close IN;

#gene de
my %gene2de;
open(IN,$genede) || die $!;
while(<IN>) {
	next if $_=~/^Feature/;
	tr/\r\n//d;
	my @array=split/\t/;
	
	$gene2de{$array[0]}=$array[5];
}
close IN;


#tx de
my %tx2de;
open(IN,$txde) || die $!;
while(<IN>) {
	next if $_=~/^Feature/;
	tr/\r\n//d;
	my @array=split/\t/;
	
	$tx2de{$array[0]}=$array[5];
}
close IN;


#tx si
my %tx2si;
open(IN,$txsi) || die $!;
while(<IN>) {
	next if $_=~/^Feature/;
	tr/\r\n//d;
	my @array=split/\t/;
	
	$tx2si{$array[0]}=$array[5];
}
close IN;

#exon de
my %exon2de;
open(IN,$exonde) || die $!;
while(<IN>) {
	next if $_=~/^Feature/;
	tr/\r\n//d;
	my @array=split/\t/;
	
	$exon2de{$array[0]}=$array[5];
}
close IN;


#exon si
my %exon2si;
open(IN,$exonsi) || die $!;
while(<IN>) {
	next if $_=~/^Feature/;
	tr/\r\n//d;
	my @array=split/\t/;
	
	$exon2si{$array[0]}=$array[5];
}
close IN;

#summarize #all exons
open(OUT,">$outputfile") || die "ERROR:Can't write to $outputfile.\n";
print OUT "Feature\tExonDE\tExonSi\tTxDE\tTxSi\tGeneDE\t$annotitle\n";

foreach my $exon (sort keys %exon2gene) {
	if(defined $exon2de{$exon} || $exon2si{$exon}) {
		print OUT $exon,"\t";
		#exonde
		if(defined $exon2de{$exon}) {
			print OUT $exon2de{$exon},"\t";
		}
		else {
			print OUT " \t";
		}
		#exonsi
		if(defined $exon2si{$exon}) {
			print OUT $exon2si{$exon},"\t";
		}
		else {
			print OUT " \t";
		}


		#exonde
		
		my @txdemarks;
		my @txsimarks;
		
		foreach my $transcript (@{$exon2tx{$exon}}) {
			if(defined $tx2de{$transcript}) {
				push @txdemarks, $tx2de{$transcript};
			}
			else {
				push @txdemarks," ";
			}
			#txsi
			if(defined $tx2si{$transcript}) {
				push @txsimarks, $tx2si{$transcript};
			}
			else {
				push @txsimarks, " ";
			}
		}
		
		
		if(@txdemarks) {
			print OUT join(",",@txdemarks),"\t";
		}
		else {
			print OUT " \t";		
		}

		if(@txsimarks) {
			print OUT join(",",@txsimarks),"\t";
		}
		else {
			print OUT " \t";		
		}

		#genede
		if(defined $gene2de{$exon2gene{$exon}}) {
			print OUT $gene2de{$exon2gene{$exon}},"\t";
		}
		else {
			print OUT " \t";
		}
		#anno
		print OUT $exon2anno{$exon},"\n";
	}
}
close OUT;
