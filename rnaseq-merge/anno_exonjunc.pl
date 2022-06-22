#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);

#CutAdapt+FASTQC+RSEM+STAR


########
#Interface
########


my $version="0.2";

#v0.2, fix bugs for names

my $usage="

anno_exonjunc
version: $version
Usage: anno_exonjunc.pl [parameters]

Description: Annotate the exon juncs

Parameters:

    --in|-i           Inputfile
    --output|-o       Outputfile

    --type            Input type, featurecounts or rnaseq-merge [rnaseq-merge]

    --tx|-t           Transcriptome
                        Currently support Human.B38.Ensembl88,Mouse.B38.Ensembl88,Rat.Rn6.Ensembl88


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
my $type="rnaseq-merge";
my $tx;
my $verbose=1;

my $dev=0; #developmental version

GetOptions(
	"in|i=s" => \$inputfile,
	"output|o=s" => \$outputfile,
	"tx|t=s" => \$tx,
	
	"type=s" => \$type,
	
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
	print STDERR "Starting analysis using $tx.\n\n" if $verbose;
	#print LOG "Starting analysis using $tx.\n\n";
}
else {
	print STDERR "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n" if $verbose;
	#print LOG "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n";
	exit;
}


########
#default ouputs
########




########
#Process
########


print STDERR "\nanno_exonjunc $version running ...\n\n" if $verbose;


#Read exon anno	
print STDERR "Reading exon annotation file:",$tx2ref{$tx}{"exonanno"},"\n\n";

my %exon2anno;
my %coord2exon;
my %tx2str;
my %tx2id;

open(IN,$tx2ref{$tx}{"exonanno"}) || die $!;
while(<IN>) {
	tr/\r\n//d;

	my @array=split/\t/;
	
	#remmember exon junc name by chr+coord
	if($array[9] eq "+") {
		$coord2exon{join("_",$array[6],$array[7])}{$array[0]."_S"}++;
		$coord2exon{join("_",$array[6],$array[8])}{$array[0]."_E"}++;
	}
	else {
		$coord2exon{join("_",$array[6],$array[7])}{$array[0]."_E"}++;
		$coord2exon{join("_",$array[6],$array[8])}{$array[0]."_S"}++;
	}
	
	$exon2anno{$array[0]}=[@array[1,2,3,4,5]];
	
	foreach my $txname (split(",",$array[5])) {
		$tx2str{$txname}=$array[9];
	}
	
	my @txids=split(",",$array[4]);
	my @txnames=split(",",$array[5]);
	
	for(my $num=0;$num<@txids;$num++){
		$tx2id{$txnames[$num]}=$txids[$num];
	}
		
}
close IN;

print STDERR "Processing exon junction file:$inputfile\n\n";

#generate annotation for exons
open(IN,$inputfile) || die $!;
open(OUT,">$outputfile") || die $!;
print OUT "Exonjunc\tExonjuncNames\tGeneID\tGeneName\tEvidences\tType\tStartExon\tEndExon\tTxID\tTxName\n";

while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Exon|PrimaryGene/;
	my @array=split/\t/;

	#decide exon junc location
	my ($ejchr,$ejstart,$ejend,$ejstr,$exonjuncid);
	
	if($type eq "rnaseq-merge") {
		($ejchr,$ejstart,$ejend,$ejstr)=split("_",$array[0]);
		$exonjuncid=$array[0];
	}
	else {
		($ejchr,$ejstart,$ejend,$ejstr)=@array[2,3,6,4];

		if($ejstr eq "NA") {
			$ejstr=".";
		}
		
		$exonjuncid=join("_",$ejchr,$ejstart,$ejend,$ejstr);
	}
	
	my %sexon;
	my %eexon;
	
	if(defined $coord2exon{join("_",$ejchr,$ejstart)}) {
		foreach my $exon (sort keys %{$coord2exon{join("_",$ejchr,$ejstart)}}) {
			$sexon{$exon}++;
		}
	}

	if(defined $coord2exon{join("_",$ejchr,$ejend)}) {
		foreach my $exon (sort keys %{$coord2exon{join("_",$ejchr,$ejend)}}) {
			$eexon{$exon}++;
		}
	}
	
	my %stx;
	my %etx;
	my %txs;
	my %genes;
	
	#decide common tx
	if(keys %sexon && keys %eexon) {
		foreach my $exon (sort keys %sexon) {
			my ($exonname,$exondir)=split("_",$exon);
			foreach my $txname (split(",",$exon2anno{$exonname}[4])) {
				if($exon2anno{$exonname}[0]=~/$txname\.(E\d+)/) {
					$stx{$txname}=[$1,$exondir];
				}
				$txs{$tx2id{$txname}}=$txname;
			}
			$genes{$exon2anno{$exonname}[1]}=$exon2anno{$exonname}[2];
		}

		foreach my $exon (sort keys %eexon) {
			my ($exonname,$exondir)=split("_",$exon);
			foreach my $txname (split(",",$exon2anno{$exonname}[4])) {
				if($exon2anno{$exonname}[0]=~/$txname\.(E\d+)/) {
					$etx{$txname}=[$1,$exondir];
				}
				$txs{$tx2id{$txname}}=$txname;
			}
			$genes{$exon2anno{$exonname}[1]}=$exon2anno{$exonname}[2];
		}			
	}
	
	my %exonjuncnames; #known exonjunc anno
	my %exonjuncevidences; #known exonjunc anno
	my %types;
	
	if(keys %stx && keys %etx) {
		#find common tx
		foreach my $txname (sort keys %stx) {
			if(defined $etx{$txname}) {
				#found common tx
				if($tx2str{$txname} eq "+") {
					$exonjuncevidences{$txname.".".$stx{$txname}[0]."_".$stx{$txname}[1]."-".$etx{$txname}[0]."_".$etx{$txname}[1]}++;
					
					if($stx{$txname}[1] eq "E" && $etx{$txname}[1] eq "S") {
						#normal junc with first end and second start
						$exonjuncnames{$txname.".".$stx{$txname}[0]."-".$etx{$txname}[0]}++;
					}
					else {
						$exonjuncnames{$txname.".".$stx{$txname}[0]."_".$stx{$txname}[1]."-".$etx{$txname}[0]."_".$etx{$txname}[1]}++;
					}
				}
				else {
					$exonjuncevidences{$txname.".".$etx{$txname}[0]."_".$etx{$txname}[1]."-".$stx{$txname}[0]."_".$stx{$txname}[1]}++;
					
					if($etx{$txname}[1] eq "E" && $stx{$txname}[1] eq "S") {
						#normal junc with first end and second start
						$exonjuncnames{$txname.".".$etx{$txname}[0]."-".$stx{$txname}[0]}++;
					}
					else {
						$exonjuncnames{$txname.".".$etx{$txname}[0]."_".$etx{$txname}[1]."-".$stx{$txname}[0]."_".$stx{$txname}[1]}++;
					}
				}
			}
		}
	}
	
	#output file
	#print OUT "Exonjunc\tExonjuncNames\tType\tStartExon\tEndExon\tTxID\tTxName\tGeneID\tGeneName\n";
	
	print OUT $exonjuncid,"\t";
	
	#assigned junc names
	if(keys %exonjuncnames) {
		print OUT join(",",sort keys %exonjuncnames),"\t";
	}
	else {
		print OUT "$exonjuncid\t";
	}
	
	#gene names
	if(keys %genes) {
		print OUT join(",",sort keys %genes),"\t";
		print OUT join(",",map {$genes{$_}} sort keys %genes),"\t";
	}
	else {
		print OUT " \t \t";
	}
	
	#evidences
	if(keys %exonjuncevidences) {
		print OUT join(",",sort keys %exonjuncevidences),"\t";
	}
	else {
		print OUT " \t";
	}
	
	#type, Known (known exons), Novel (no exonjuncnames assigned)
	if(keys %exonjuncnames) {
		print OUT "KnownJunctions\t";
	}
	else {
		print OUT "Novel\t";
	}
	
	
	
	if(keys %sexon) {
		print OUT join(",",sort keys %sexon),"\t";
	}
	else {
		print OUT " \t";
	}
	
	if(keys %eexon) {
		print OUT join(",",sort keys %eexon),"\t";
	}
	else {
		print OUT " \t";
	}
	
	if(keys %txs) {
		#sort by txnames
		print OUT join(",",sort {$txs{$a} cmp $txs{$b}} keys %txs),"\t";
		print OUT join(",",map {$txs{$_}} sort {$txs{$a} cmp $txs{$b}} keys %txs),"\n";
	}
	else {
		print OUT " \t \n";
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

