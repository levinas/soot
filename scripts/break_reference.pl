#! /usr/bin/env perl

use strict;
use Data::Dumper;

my $usage = "Usage: $0 ref.gff depth \n\n";

my $ref_file   = shift @ARGV;
my $depth_file = shift @ARGV;

my $gff = read_gff($ref_file);
my @depth_list = map { chomp; [ split /\t/ ] } `cat $depth_file`;
my $median_depth = median(map { $_->[2] } @depth_list);

print "$median_depth\n";


sub read_gff {
    my ($file) = @_;
    my $header = `cat $file | grep "^#"`;
    my @lines = `cat $file | grep -v "^#"`;
    my @features;
    for (@lines) {
        chomp;
        my ($seqname, $source, $feature, $start, $end, $score, $strand, $fname, $attribute) = split /\t/;
        my %hash = map { my ($k,$v) = split /=/; $k => $v } split(/;\s*/, $attribute);
        push @features, [ seqname => $seqname,
                          source => $source,
                          feature => $feature,
                          start => $start,
                          end => $end,
                          score => $end,
                          strand => $strand,
                          fname => $fname,
                          attribute => \%hash ];
    }
    return \@features;
}

sub median {
    return unless @_;
    my @array = sort { $a <=> $b } @_;
    my $index = int($#array / 2);
    return $array[$index];
}
