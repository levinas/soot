#! /usr/bin/env perl

use strict;
use Carp;
use Data::Dumper;
use gjoseqlib;

my $usage = "Usage: $0 gene-detail.txt flank \n\n";

my $file  = shift @ARGV or die $usage;
my $flank = shift @ARGV;

my @data = map { [ split/\t/ ] } map { chomp; /^\S/ ? $_ : () } `cat $file`;

run("mkdir -p seqs");
chdir("seqs");

my @seq_all;

my $num_strains = 3;

my %by_group;
my %by_org;

while (@data > 0) {
    my $group;
    my $org;
    for (my $i = 0; $i < $num_strains; $i++) {
        print "$i\n";
        my $strain = shift @data or last;
        my $gene = $strain->[0];
        $group = $gene if $gene =~ /^YPO/;
        $org = substr($gene, 0, 3);
        my $peg = $strain->[6];
        my $loc  = $strain->[1];
        my $seq0 = loc_to_seq($loc, 0, $gene, $peg);
        my $seq  = loc_to_seq($loc, $flank, $gene, $peg);
        push @{$by_group{$group}}, $seq0;
        push @{$by_org{$org}}, $seq;
        # print STDERR '\%by_org = '. Dumper(\%by_org);
        # print STDERR '\%by_group = '. Dumper(\%by_group);
    }
}

for (keys %by_group) {
    my $fname = "$_.fa";
    write_fasta($fname, $by_group{$_});
}

for (keys %by_org) {
    my $fname = $_;
    $fname .= ".f$flank" if $flank > 0;
    $fname .= ".fa";
    write_fasta($fname, $by_org{$_});
}

sub loc_to_seq {
    my ($loc, $flank, $id, $peg) = @_;
    my $func = `echo "$peg" | svr_function_of|cut -f2`; chomp($func);
    my ($gid) = $loc =~ /^(\d+\.\d+)/;
    my ($ctg, $beg, $strand, $len) = $loc =~ /(\S+)_(\d+)([+-])(\d+)/;
    my $new_beg = $strand eq '+' ? $beg - $flank : $beg + $flank;
    my $new_len = $len + 2 * $flank;
    my $new_loc = "$ctg\_$new_beg$strand$new_len";
    my $org_name = `echo $gid|svr_genome_statistics name|cut -f2`; chomp($org_name);
    my $s = `echo "$new_loc" |svr_dna_seq |cut -f2`; chomp($s);
    substr($s, $flank, $len) = uc substr($s, $flank, $len);
    my $def = "$peg [$org_name] [$func]";
    $def = $flank > 0 ? "$new_loc [$loc with $flank flanking bases] $def" : "$loc $def";
    [ $id, $def, $s ];
}

sub run { system(@_) == 0 or confess("FAILED: ". join(" ", @_)); }
