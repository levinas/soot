#! /usr/bin/env perl

use strict;
use Data::Dumper;

my $usage = "Usage: $0 contigs.fa\n\n";

my $contig_file = shift @ARGV or die $usage;
my $gene_file = 'gene-detail.txt';

my @genes = map { chomp; $_ } `cut -f1 $gene_file |grep "^YPO"`;
my %lenH = map { chomp; my ($gene, $len) = split /\t/;
                 $gene ? ($gene => $len) : () } `cut -f1,6 $gene_file |grep "^YPO"`;

@genes[1,2] = @genes[2,1];

my @lines = `blastn -query seqs/YPA.f500.fa -db $contig_file -outfmt 6 |cut -f1,3,4,7,8`;

my %seen;
for (@lines) {
    chomp;
    my ($gene, $ident, $bitscore, $q1, $q2) = split /\t/;
    next if $seen{$gene};
    my $len = $lenH{$gene};
    my $cov1 = $q1 < 500 ? (500 - $q1 + 1) : 0; 
    my $cov2 = $q2 > 500 + $len ? ($q2 - 500 - $len) : 0;
    print join("\t", $_, $cov1, $cov2) . "\n";
    $seen{$gene} = 1;
}
