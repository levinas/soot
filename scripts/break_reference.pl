#! /usr/bin/env perl

use strict;
use Data::Dumper;

my $usage = "Usage: $0 ref.gff depth \n\n";

my $ref_file   = shift @ARGV;
my $depth_file = shift @ARGV;

my $gff = read_gff($ref_file);
my @depth_list = map { chomp; [ split /\t/ ] } `cat $depth_file|head -n100`;
# my $median_depth = median(map { $_->[2] } @depth_list);
my $median_depth = 126;
my $thresh = $median_depth / 10;


my ($db1, $de1, $db2, $de2);
my $p = 1;                      # position in gff
for (@depth_list) {
    my ($contig, $pos, $cov) = @$_;
    my $np = next_pos($gff, $p);
    while ($np && $gff->[$np]->{end} <= $pos) {
        $p = $np;
        $np = next_pos($gff, $p);
    }
    last unless $np;
    # print STDERR '$gff->[$p] = '. Dumper($gff->[$p]);

    $db1 = $gff->[$p]->{start} - $pos;
    $de1 = $gff->[$p]->{end} - $pos;
    $db2 = $gff->[$np]->{start} - $pos;
    $de2 = $gff->[$np]->{end} - $pos;
    print join("\t", $contig, $pos, $cov, $db1, $de1, $db2, $de2) . "\n";
}

sub advance_gff {
    my ($gff, $pos) = @_;
}

sub next_pos {
    my ($gff, $pos) = @_;
    my $d = 1;
    $d++ while $pos+$d < @$gff && $gff->[$pos+$d]->{attribute}->{Parent};
    return undef if $pos+$d < @$gff && $gff->[$pos+$d]->{attribute}->{Parent};
    return $pos+$d;
}

sub read_gff {
    my ($file) = @_;
    my $header = `cat $file | grep "^#"`;
    my @lines = `cat $file | grep -v "^#"`;
    my @features;
    for (@lines) {
        chomp;
        my ($seqname, $source, $feature, $start, $end, $score, $strand, $fname, $attribute) = split /\t/;
        my %hash = map { my ($k,$v) = split /=/; $k => $v } split(/;\s*/, $attribute);
        push @features, { seqname => $seqname,
                          source => $source,
                          feature => $feature,
                          start => $start,
                          end => $end,
                          score => $end,
                          strand => $strand,
                          fname => $fname,
                          attribute => \%hash };
    }
    return \@features;
}

sub median {
    return unless @_;
    my @array = sort { $a <=> $b } @_;
    my $index = int($#array / 2);
    return $array[$index];
}
