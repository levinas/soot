#! /usr/bin/env perl

use strict;

my $usage = "Usage: $0 < reads.gfa > contigs.fa\n\n";

while (<>) {
    next unless /^S\s/;
    chomp;
    my ($type, $segName, $seqSeq) = split /\t/;
    print ">$segName\n$seqSeq\n";
}
