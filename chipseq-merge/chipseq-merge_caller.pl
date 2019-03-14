#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

#CutAdapt+FASTQC+RSEM+STAR


########
#Prerequisites
########

my $multiqc="/home/jyin/.local/bin/multiqc";
my $mergefiles="perl /home/jyin/Projects/Pipeline/sbptools/mergefiles/mergefiles_caller.pl";

my $mergebed="/apps/bedtools2-2.26.0/bin/bedtools merge";
my $renamebed="perl /home/jyin/Projects/Pipeline/sbptools/chipseq-merge/rename_bed.pl";
my $selectbed="perl /home/jyin/Projects/Pipeline/sbptools/chipseq-merge/select_mergebed.pl";

#homer location
my $homer="/home/jyin/Programs/Homer/bin";
my $findpeaks="$homer/findPeaks";
my $makeucscfile="$homer/makeUCSCfile";
my $pos2bed="$homer/pos2bed.pl";
my $bed2pos="$homer/bed2pos.pl";
my $annotatepeaks="$homer/annotatePeaks.pl";

#UCSC
my $bedtobigbed="/apps/ucsc/bedToBigBed";


########
#Interface
########


my $version="0.1";

my $usage="

chipseq-merge
version: $version
Usage: sbptools chipseq-merge [parameters]

Description: Merge chipseq-process folder to get summarized QC, counting for promoter and merged peak, and so on

Parameters:

    --in|-i           Input folder(s)
	
    --config|-c       Configuration file #at this stage, configuration only, to make merged peaks
    --group|-g        Column name for sample groups to merge peaks

    --output|-o       Output folder


    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84
    --anno|-a         Add annotation

    --runmode|-r      Where to run the scripts, local, server or none [none]
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

my $samples;
my $inputfolders;

my $group;

my $configfile;
my $outputfolder;
my $verbose;
my $tx;
my $runmode="none";

GetOptions(
	"in|i=s" => \$inputfolders,
	#"samples|s=s" => \$samples,	
	"config|c=s" => \$configfile,
	"group|g=s" => \$group,
	"output|o=s" => \$outputfolder,
	"tx|t=s" => \$tx,
	"runmode|r=s" => \$runmode,		
	"verbose|v" => \$verbose,
);


#default ouputs

my $promotermergedcount_norm="promoter.merged.norm.count.txt"; #count ,annotated count, 
my $promotermergedcount_raw="promoter.merged.raw.count.txt"; #count ,annotated count, 
#my $samplemergedcount="promoter.merged.count.txt"; #peaks for samples from the same group are merged 
my $allmergedcount_norm="all.reprod.peak.merged.norm.count.txt"; #peaks for all samples are merged, based on group merged
my $allmergedcount_raw="all.reprod.peak.merged.raw.count.txt"; #peaks for all samples are merged, based on group merged


my $allmergedbed="all.reprod.peak.merged.bed";
my $allmergedbed_sorted="all.reprod.peak.merged_sorted.bed";
my $allmergedpos="all.reprod.peak.merged.pos";
my $allmergedbb="all.reprod.peak.merged_sorted.bb";

#Create folders

$outputfolder = abs_path($outputfolder);

if(!-e $outputfolder) {
	mkdir($outputfolder);
}


my $scriptfolder="$outputfolder/scripts";
my $tempfolder="$outputfolder/temp";

if(!-e $scriptfolder) {
	mkdir($scriptfolder);
}

if(!-e $tempfolder) {
	mkdir($tempfolder);
}

my $logfile="$outputfolder/chipseq-merge_run.log";

my $scriptfile1="$scriptfolder/chipseq-merge_run1.sh";
my $scriptfile2="$scriptfolder/chipseq-merge_run2.sh";


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "$multiqc version:", getsysoutput("$multiqc --version"),"\n";
print LOG "\n";


print STDERR "\nsbptools chipseq-merge version $version running ...\n\n" if $verbose;
print LOG "\nsbptools chipseq-merge version $version running ...\n\n";

#test tx option

my %tx2ref=(
	"Human.B38.Ensembl84"=> { 
		"star"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR",
		"chrsize"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_homeranno.txt",},
	"Mouse.B38.Ensembl84"=>{ 
		"star"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR",
		"chrsize"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_homeranno.txt",}
);



#to genome version
my %tx2gv=(
	"Human.B38.Ensembl84"=>"hg38",
	"Mouse.B38.Ensembl84"=>"mm10",
);

my %tx2promoter=(
	"Human.B38.Ensembl84"=> {
		"1000u0d_longest"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_promoter_1kup0down_longesttxs.pos",
		"1000u0d_all"=> "/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_promoter_1kup0down_alltxs.pos"
	},
	"Mouse.B38.Ensembl84"=> {
		"1000u0d_longest"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_promoter_1kup0down_longesttxs.pos",
		"1000u0d_all"=> "/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_promoter_1kup0down_alltxs.pos"
	}
);

my %tx2gtf=(
	"Human.B38.Ensembl84"=>"/home/jyin/Data/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc.gtf",
	"Mouse.B38.Ensembl84"=>"/home/jyin/Data/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc.gtf",
);


if(defined $tx2ref{$tx}) {
	print STDERR "Starting analysis using $tx.\n\n" if $verbose;
	print LOG "Starting analysis using $tx.\n\n";
}
else {
	print STDERR "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n" if $verbose;
	print LOG "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n";
}


#Files for merging
#Examples
#JY_315_Isl.genes.results
#JY_315_Isl.isoforms.results


########
#Process
########


#open config files to find samples

my @samples_array;
my %samples_hash; #check no dup samples
my $groupcol;
my %group2samples;

if(defined $samples && length($samples)>0) {
	foreach my $sample (split(",",$samples)) {
		#keep the order of samples
		push @samples_array, $sample;
		unless (defined $samples_hash{$sample}) {
			$samples_hash{$sample}++
		}
		else {
			print STDERR "ERROR:$sample is defined multiple times in $samples.\n";
			print LOG "ERROR:$sample is defined multiple times in $samples.\n";
			exit;
		}
	}
}
elsif(defined $configfile && length($configfile)>0) {
	#another way of getting samples
	
	open(IN,$configfile) || die "Error reading $configfile. $!";
	#first column should be sample, with a header

	my $fileline=0;
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		
		if($fileline==0) {
			#title			
			for(my $num=0;$num<@array;$num++) {
				if($array[$num] eq $group) {
					$groupcol=$num;
					last;
				}
			}
			
			print STDERR "--group $group is identifed at column ",$groupcol+1," in $configfile.\n\n" if $verbose;
			print LOG "--group $group is identifed at column ",$groupcol+1," in $configfile.\n\n";
			
		}
		else {
			my $sample=$array[0];
			push @samples_array, $sample;
			unless (defined $samples_hash{$sample}) {
				$samples_hash{$sample}++;
				$group2samples{$array[$groupcol]}{$sample}++;
			}
			else {
				print STDERR "ERROR:$sample is defined multiple times in $configfile.\n";
				print LOG "ERROR:$sample is defined multiple times in $configfile.\n";
				exit;
			}
		}
		
		$fileline++;
	}
	close IN;
}

print STDERR scalar(@samples_array)," samples identified: ",join(",",@samples_array),"\n\n" if $verbose;
print LOG scalar(@samples_array)," samples identified: ",join(",",@samples_array),"\n\n";

print STDERR scalar(keys %group2samples)," groups of samples identified: ",join(",",sort keys %group2samples),"\n\n" if $verbose;
print LOG scalar(keys %group2samples)," groups of samples identified: ",join(",",sort keys %group2samples),"\n\n";


#----------------
#Find files for merging
#Examples
#JY_315_Isl.genes.results
#JY_315_Isl.isoforms.results

print STDERR "\nReading sample folders.\n" if $verbose;
print LOG "\nReading sample folders.\n";

#my %sample2gene;
#my %sample2tx;
my %sample2folder;

my %sample2tagdir;

foreach my $infolder (split(",",$inputfolders)) {
	#samples in different folders
	#if different folders have the same sample name, use the first found sample
	
	my @infiles=glob("$infolder/*");
	
	foreach my $file (@infiles) {
		if(-d $file) {
			my $samplename=basename($file);
			
			#see whether it is a defined sample
			if(defined $samples_hash{$samplename}) {
				print STDERR $samplename," is found in ",$file,".\n" if $verbose;
				print LOG $samplename," is found in ",$file,".\n" if $verbose;
				
				$sample2folder{$samplename}=abs_path($file);
				
				my @samplefiles=glob("$file/*");
				
				foreach  my $samplefile (@samplefiles) {
					#homer TagDir
					if($samplefile=~/_TagDir$/) {
						$sample2tagdir{$samplename}=abs_path($samplefile);
					}
				}
			}
		}
	}
}



print STDERR scalar(keys %sample2tagdir)," samples identified with Homer TagDir.\n\n" if $verbose;
print LOG scalar(keys %sample2tagdir)," samples identified with Homer TagDir.\n\n";


if( scalar(@samples_array) != scalar(keys %sample2tagdir) ) {
	print STDERR "ERROR:Not all samples have Homer TagDir.\n\n";
	print LOG "ERROR:Not all samples have Homer TagDir.\n\n";
	exit;
}



########
#Print out commands, for local and server run
########

open(S1,">$scriptfile1") || die "Error writing $scriptfile1. $!";


#----------------
#Multiqc with specific folders


#file for multiqc folder search
open(OUT,">$tempfolder/samplefolders.txt") || die $!;
foreach my $sample (@samples_array) {
	print OUT $sample2folder{$sample},"\n";
}
close OUT;

print S1 "$multiqc -l $tempfolder/samplefolders.txt -o $outputfolder/multiqc;\n";

#------------------


#Get peaks for promoters

#all tag dirs
my $tagdirs=join(" ",map {$sample2tagdir{$_}} sort keys %sample2tagdir);

foreach my $pv (sort keys %{$tx2promoter{$tx}}) {
	#promoter version
	
	print STDERR "Printing commands to annotate promoter based on ",$tx2promoter{$tx}{$pv},"\n\n" if $verbose;
	print LOG "Printing commands to annotate promoter based on ",$tx2promoter{$tx}{$pv},"\n\n";
	
	#command to get counting for promoters
	print S1 "$annotatepeaks ",$tx2promoter{$tx}{$pv}," ",$tx2gv{$tx}," -gtf ",$tx2gtf{$tx}," -d $tagdirs > $outputfolder/$pv\_$promotermergedcount_norm 2> $pv\_promoter_annotatepeaks_norm_run.log\n";
	print S1 "$annotatepeaks ",$tx2promoter{$tx}{$pv}," ",$tx2gv{$tx}," -gtf ",$tx2gtf{$tx}," -d $tagdirs -raw > $outputfolder/$pv\_$promotermergedcount_raw 2> $pv\_promoter_annotatepeaks_raw_run.log\n";
}

#Merge called peaks for sample groups
my @groupselectedbeds;

foreach my $group (sort keys %group2samples) {
	#1) merge bed first
	#2) identify reproducible bed
	#3) get counting and annotation
	
	#get peak bed files for each sample
	
	my @peakfiles;
	my @renamed;
	my $cmd;

	my $mergedbed="$outputfolder/$group\_merged.bed";
	my $selectedbed="$outputfolder/$group\_reprod.bed";
	my $selectedbed_sorted="$outputfolder/$group\_reprod_sorted.bed";
	my $selectedbb="$outputfolder/$group\_reprod_sorted.bb";
	
	push @groupselectedbeds,$selectedbed;

	print STDERR "Producing reproducible merged bed for $group, as $selectedbed\n" if $verbose;
	print LOG "Producing reproducible merged bed for $group, as $selectedbed\n";	

	#rename bed
	foreach my $sample (sort keys %{$group2samples{$group}}) {
		my $peakfile=$sample2folder{$sample}."/$sample\_Peaks.bed";
		my $peakfilerenamed="$tempfolder/".basename($peakfile);
		$peakfilerenamed=~s/.bed$/_renamed.bed/;
		
		push @peakfiles, $peakfile;
		push @renamed, $peakfilerenamed;
		
		print S1 "$renamebed -i $peakfile -o $peakfilerenamed -s $sample;"
	}
	
	#merge bed 
	print S1 "cat ",join(" ",@renamed)," | sort -k1,1 -k2,2n | bedtools merge -c 4 -o collapse > $mergedbed;";
	
	#select reproducible bed
	print S1 "$selectbed $mergedbed $selectedbed;";

	#produce bb
	print S1 "cat $selectedbed | grep -v \"#\" | sort -k1,1 -k2,2n > $selectedbed_sorted;$bedtobigbed $selectedbed_sorted ".$tx2ref{$tx}{"chrsize"}." $selectedbb;";
	
	#convert 2 pos? and annotate ...
	

	print S1 "\n";
}

close S1;

#All merge
#Merge called peaks for sample groups

print STDERR "Producing merged bed for all groups, as $allmergedbed\n" if $verbose;
print LOG "Producing merged bed for all groups, as $allmergedbed\n";	


open(S2,">$scriptfile2") || die "Error writing $scriptfile2. $!";

#merge all reprod peaks
print S2 "cat ",join(" ",@groupselectedbeds)," | sort -k1,1 -k2,2n | $mergebed -c 4 -o collapse > $tempfolder/$allmergedbed\_original;";

#rename
print S2 "$renamebed -i $tempfolder/$allmergedbed\_original -o $allmergedbed -c;";

#bb
print S2 "cat $allmergedbed | grep -v \"#\" | sort -k1,1 -k2,2n > $allmergedbed_sorted;$bedtobigbed $allmergedbed_sorted ".$tx2ref{$tx}{"chrsize"}." $allmergedbb;";
	
#convert to pos
print S2 "$bed2pos $allmergedbed -o $allmergedpos;";

#annotate peaks
print S2 "$annotatepeaks ",$allmergedpos," ",$tx2gv{$tx}," -gtf ",$tx2gtf{$tx}," -d $tagdirs > $outputfolder/$allmergedcount_norm 2> allmergedcount_annotatepeaks_norm_run.log;";
print S2 "$annotatepeaks ",$allmergedpos," ",$tx2gv{$tx}," -gtf ",$tx2gtf{$tx}," -d $tagdirs -raw > $outputfolder/$allmergedcount_raw 2> allmergedcount_annotatepeaks_raw_run.log;";

print S2 "\n";

close S2;

#local mode
print STDERR "\nTo run locally, in shell type: sh $scriptfile1;sh $scriptfile2\n\n";
print LOG "\nTo run locally, in shell type: sh $scriptfile1;sh $scriptfile2\n\n";



close LOG;
########
#Functions
########

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


