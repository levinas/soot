#! /usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;

use BlastInterface;
use gjoalignment;
use gjoseqlib;
use Proc::ParallelLoop;

my $usage = "Usage: $0 ref.fa reads.fq \n\n";

my ($help, $evalue, $minIden, $minCovQ, $threads, $prefix);

GetOptions("h|help"      => \$help,
           "e|evalue=f"  => \$evalue,
           "i|minIden=f" => \$minIden,
           "c|minCovQ=f" => \$minCovQ,
           "t|threads=i" => \$threads,
           "pre=s"       => \$prefix,
	  ) or die("Error in command line arguments\n");

$help and die $usage;

my $ref   = shift @ARGV;
my $reads = shift @ARGV;

$ref && $reads && -s $ref && -s $reads or die $usage;

$evalue  ||= 10.0;
$minIden ||= 0.8;
$minCovQ ||= 0.8;
$threads ||= 8;

if (!$prefix) {
    $prefix = join('_vs_', $reads, $ref);
    $prefix =~ s/\.(fasta|fa)//g;
}

my $opts = { evalue => $evalue, minIden => $minIden, minCovQ => $minCovQ, threads => $threads, outForm => 'hsp', blastplus => 1 };
my @matches = BlastInterface::blastn($reads, $ref, $opts);
# print STDERR '\@matches = '. Dumper(\@matches);

open(F, ">$prefix.blastn.hsp") or die "Could not open $prefix.blastx.hsp";
for (@matches) { print F join("\t", @$_) . "\n"; }
close(F);

