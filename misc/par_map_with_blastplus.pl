#! /usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;

use BlastInterface;
use gjoalignment;
use gjoseqlib;
use Proc::ParallelLoop;

my $usage = "Usage: $0 ref.fa reads.fq \n\n";

my ($help, $blocksize, $evalue, $minIden, $minCovQ, $threads, $prefix);

GetOptions("h|help"        => \$help,
           "b|blocksize=i" => \$blocksize,
           "e|evalue=f"    => \$evalue,
           "i|minIden=f"   => \$minIden,
           "c|minCovQ=f"   => \$minCovQ,
           "t|threads=i"   => \$threads,
           "pre=s"         => \$prefix,
	  ) or die("Error in command line arguments\n");

$help and die $usage;

my $ref   = shift @ARGV;
my $reads = shift @ARGV;

$ref && $reads && -s $ref && -s $reads or die $usage;

$evalue    ||= 10.0;
$minIden   ||= 0.8;
$minCovQ   ||= 0.8;
$threads   ||= 8;
$blocksize ||= 10000;

if (!$prefix) {
    $prefix = join('_vs_', $reads, $ref);
    $prefix =~ s/\.(fasta|fa)//g;
}

# my @query = gjoseqlib::read_fasta($reads);
# strange error: panic: sv_setpvn called with negative strlen -4294967142 at /space2/fangfang/sas/lib/gjoseqlib.pm line 370.
my @query = gjoseqlib::read_fasta_0($reads);

my $opts = { evalue => $evalue, minIden => $minIden, minCovQ => $minCovQ, threads => 1, outForm => 'hsp', blastplus => 1 };

my $beg = 0;
my $end;

my $blocks = int((@query - 1)/$blocksize) + 1;
pareach [ 1..$blocks ], sub {
    my $i = shift;
    my $beg = ($i-1) * $blocksize;
    my $end = $i * $blocksize - 1;
    # print STDERR join("\t", $beg, $end) . "\n";
    $end = @query - 1 if $end >= @query;
    my @subset = @query[$beg..$end];
    my @matches = BlastInterface::blastn(\@subset, $ref, $opts);
    for (@matches) { print join("\t", @$_) . "\n"; }
}, { Max_Workers => $threads };

# open(F, ">$prefix.blastn.hsp") or die "Could not open $prefix.blastx.hsp";
# for (@matches) { print F join("\t", @$_) . "\n"; }
# close(F);

