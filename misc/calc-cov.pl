#! /usr/bin/env perl

use strict;
use Data::Dumper;

my $usage = "Usage: $0 bwa-sample-YPO.500/aln.mpileup\n\n";

my $mpileup = shift @ARGV or die $usage;
my $gene_file = 'gene-detail.txt';

my @genes = map { chomp; $_ } `cut -f1 $gene_file |grep "^YPO"`;
my %lenH = map { chomp; my ($gene, $len) = split /\t/;
                 $gene ? ($gene => $len) : () } `cut -f1,6 $gene_file |grep "^YPO"`;

@genes[1,2] = @genes[2,1];
my @folds = (1, 2, 3);

for my $gene (@genes) {
    my $len = $lenH{$gene};
    my $beg = 1 + 500;
    my $end = $len + 500;
    my ($seen, @counts);
    my @lines = `cut -f1-4 $mpileup |grep "^$gene"`;
    my @covs;
    for (@lines) {
        chomp;
        my ($g, $pos, $base, $dp) = split /\t/;
        next if $pos < $beg || $pos > $end;
        push @covs, $dp;
        $seen++;
        for (my $i = 0; $i < @folds; $i++) {
            my $x = $folds[$i];
            $counts[$i]++ if $dp >= $x;
        }
    }
    my $uncov = $len - $seen;
    my $medcov = median(@covs);
    print join("\t", $gene, "len = $len", "med = $medcov", "uncov = $uncov"), ", ";
    for (my $i = 0; $i < @folds; $i++) {
        my $x = $folds[$i];
        my $cnt = $counts[$i];
        my $frac = sprintf "%.1f", ($cnt / $len) * 100;
        # print join("\t", "cov $x"."x = $cnt", "fcov $x"."x = $frac"). ", ";
        print join("\t", "fcov $x"."x = $frac"). ", ";
    }        
    print "\n";
}

sub median {
    my @x = @_;
    return unless @x;
    @x = sort { $a <=> $b } @x;
    return $x[int(@x/2)];
}

