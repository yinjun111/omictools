#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);



########
#Interface
########


my $version="0.2";

#0.1 only summarize txde,txsi,genede
#0.2 add exon


my $usage="

summarize_txsi
version: $version
Usage: perl summarize_txsi [parameters]

Description: summarize tx alternative splicing results


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

#tx anno
my %tx2gene;
my %tx2anno;
my %tx2name;
open(IN,$tx2ref{$tx}{"txanno"}) || die $!;
while(<IN>) {
	next if $_=~/^Transc/;
	tr/\r\n//d;
	my @array=split/\t/;
	
	$tx2gene{$array[0]}=$array[2];
	$tx2anno{$array[0]}=join("\t",@array[1,2,3]);
	$tx2name{$array[0]}=$array[1];
}
close IN;

#exon anno
my %tx2exon;
my %exon2name;
open(IN,$tx2ref{$tx}{"exonanno"}) || die $!;
while(<IN>) {
	next if $_=~/^Exon/;
	tr/\r\n//d;
	my @array=split/\t/;
	
	foreach my $transcript (split(",",$array[4])) {
		
		my $txname=$tx2name{$transcript};
		
		if($array[1]=~/$txname\.E0*(\d+)/) {
			#print $1,"\t",$&,"\n";
			$tx2exon{$transcript}{$1}=$array[0];
			$exon2name{$array[0]}=$&;
		}
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


#summarize #all txs
open(OUT,">$outputfile") || die "ERROR:Can't write to $outputfile.\n";
print OUT "Feature\tTxDE\tTxSI\tExonDE\tExonSI\tGeneDE\tTranscript_name\tGene\tgene_name\tExon\texon_name\n";

foreach my $txn (sort keys %tx2gene) {
	if(defined $tx2de{$txn} || $tx2si{$txn}) {
		print OUT $txn,"\t";
		#txde
		if(defined $tx2de{$txn}) {
			print OUT $tx2de{$txn},"\t";
		}
		else {
			print OUT " \t";
		}
		#txsi
		if(defined $tx2si{$txn}) {
			print OUT $tx2si{$txn},"\t";
		}
		else {
			print OUT " \t";
		}
		
		my @exondemarks;
		my @exonsimarks;
		my @exons;
		my @exonnames;
		
		if(defined $tx2exon{$txn}) {
			foreach my $exonnum (sort {$a<=>$b} keys %{$tx2exon{$txn}}) {
				
				my $exon=$tx2exon{$txn}{$exonnum};
				
				if(defined $exon2de{$exon}) {
					push @exondemarks,$exon2de{$exon};
				}
				else {
					push @exondemarks," ";
				}
				
				if(defined $exon2si{$exon}) {
					push @exonsimarks,$exon2si{$exon};
				}
				else {
					push @exonsimarks," ";
				}
				
				push @exonnames,$exon2name{$exon};
			}
		}
		
		if(@exondemarks) {
			print OUT join(",",@exondemarks),"\t";
		}
		else {
			print OUT " \t";
		}
		
		if(@exonsimarks) {
			print OUT join(",",@exonsimarks),"\t";
		}
		else {
			print OUT " \t";
		}
		
		
		#genede
		if(defined $gene2de{$tx2gene{$txn}}) {
			print OUT $gene2de{$tx2gene{$txn}},"\t";
		}
		else {
			print OUT " \t";
		}
		
		#anno
		print OUT $tx2anno{$txn},"\t";
		
		if(@exons) {
			print OUT join(",",@exons),"\t";
		}
		else {
			print OUT " \t";
		}
		
		if(@exonnames) {
			print OUT join(",",@exonnames),"\n";
		}
		else {
			print OUT " \n";
		}
		
	}
}
close OUT;
