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

for my $gene (@genes) {
    my $len = $lenH{$gene};
    my $beg = 1 + 500;
    my $end = $len + 500;
    my %seen = map { chomp; $_ => 1 } `grep "^$gene" $mpileup |cut -f2`;
    my $cnt1 = 0;               # 5' end contiguous covered bases
    my $cnt2 = 0;               # 3' end contiguous covered bases
    # print STDERR '\%seen = '. Dumper(\%seen);
    
    for (my $i = $beg-1; $i >= 1; $i--) {
        $cnt1++ if $seen{$i};
        last if !$seen{$i};
    }
    for (my $i = $end+1; $i <= $end+500; $i++) {
        $cnt2++ if $seen{$i};
        last if !$seen{$i};
    }
    print join("\t", $gene, $cnt1, $cnt2) . "\n";
}
print "\n";

sub median {
    my @x = @_;
    return unless @x;
    @x = sort { $a <=> $b } @x;
    return $x[int(@x/2)];
}

