#!/usr/bin/perl -w
##########################
#
# USE AT YOUR OWN RISK.
#
# 
#
# This script will read a set of bam files
# If a read maps to a contig, add 1 to the counts of contig hash
#
# Map the contigs to the Bins.
# Reconcile the Bins to Samples mapping.
# Print output matrix.
#
# Ambiguous read mapping (reads mapping to more than 1 contig)
# might affect results
#
##########################

use strict;
use warnings;
my %binmap = (); # $binmap{bin}{contig}
my %readmap = (); # $readmap{Sample}{Contig}
my %contig_to_bins=();

# Map the reads to contigs
my @bamfiles = </share/eisenlab/gjospin/Seagrass_JGI/bowtie2_wd/*.bam>;
foreach my $bamfile(@bamfiles){
    $bamfile=~ m/bowtie2_wd\/ZM.(\S+).bam/;
    my $core=$1;
    print STDERR "Processing $bamfile\n";
    next if $core =~ /sorted/;
    open(INBAM,"samtools view $bamfile |");
    while(<INBAM>){
	chomp($_);
	$_ =~ m/^(\S+)\s+(\S+)\s+(\S+)/;
	my $read=$1;
	my $bits =$2;
	my $contig=$3;

	if($contig eq "-"){
	    next;
	}else{
	    $readmap{$core}{$contig}=0 unless defined $readmap{$core}{$contig};
	    $readmap{$core}{$contig}++;
	}
    }
#    last;
}

# Map the contigs to Samples
my @binfiles=</share/eisenlab/gjospin/Seagrass_JGI/SGJGI_dastool_summary.phylum.1.0/bin_by_bin/*/*-contigs.fa>;

foreach my $binfile(@binfiles){

    $binfile=~ m/bin_by_bin\/(\S+)\/\S+-contigs.fa$/;
    my $bincore= $1;
    open(INCONTIG,$binfile);
    while(<INCONTIG>){
	chomp($_);
	if($_=~ m/^>(\S+)/){
	    my $contig=$1;
#	    print STDERR "Contig:$contig\n";
	    $contig_to_bins{$contig}=$bincore;
	}
    }
    close(INCONTIG);
}
my %samples=();
# readmap { sample } { contig }
# merge info and make into a single data frame
my %bin_to_samples=();
foreach my $sample(keys %readmap){

    foreach my $contig(keys %{$readmap{$sample}}){
	$samples{$sample}=1;
	if(exists $contig_to_bins{$contig}){
	    $bin_to_samples{$contig_to_bins{$contig}}{$sample}=0 unless defined $bin_to_samples{$contig_to_bins{$contig}}{$sample};
	    $bin_to_samples{$contig_to_bins{$contig}}{$sample}+=$readmap{$sample}{$contig};
	}
    }
}


# The following block is not needed for Counts but needed if doing Relative Abundance
#
# modify values and print into readable table
# 
# Get read count values from read_count.txt files
# These were the input files for the co-assembly
# Pre-computed so we just read a small text file

my %qc_reads=();
open(INCOUNTQC,"/share/eisenlab/gjospin/Seagrass_JGI/ZM_Plantgenomes_QCfiles/read_counts_QC.txt");
while(<INCOUNTQC>){
    chomp($_);
    $_ =~ m/^ZM\.(\S+)\.filt\.fastq\.gz\s+(\d+)$/;
    my $sample=$1;
    my $qc= $2;
    $qc_reads{$sample}=$qc*2; # the number reported was pairs of reads.
    # We need total number of reads for our purpose.
}
close(INCOUNTQC);
my @s = sort keys %samples;
my $ssize= scalar(@s);
#open(OUT,">SGJGI.Dastool.Bins.RelAbundance.bySamples.txt");

# Printing the output matrix file

open(OUT,">SGJGI.Dastool.Bins.counts.bySamples.txt");
print OUT "Bins/Samples";

for(my $i=0;$i<$ssize;$i++){
    print OUT "\t$s[$i]";
}
print OUT "\n";
foreach my $bin(keys %bin_to_samples){
    print OUT "$bin";

    for(my $i=0;$i<$ssize;$i++){
	unless(exists $bin_to_samples{$bin}{$s[$i]}){
	    print OUT "\t0";
	}else{
	    #my $relab = $bin_to_samples{$bin}{$s[$i]} / $qc_reads{$s[$i]};
	    #print OUT "\t$relab";
	    print OUT "\t$bin_to_samples{$bin}{$s[$i]}";
	}
    }
    print OUT "\n";
}


close(OUT);
