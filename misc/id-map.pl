#! /usr/bin/env perl

use strict;

my $usage = "Usage: $0 args\n\n";

while (<>) {
    chomp;
    my ($refseq) = split /\s+/;
    if ($refseq =~ /^YP[AO]/) {
        my ($ctg, $beg, $end) = refseq_to_ncbi_gene($refseq);
        my ($peg) = refseq_to_peg($refseq);
        print join("\t", $refseq, $ctg, $beg, $end, $peg) . "\n";
        # last;
    }
    else {
        print "$_\n";
    }
}

sub refseq_to_peg {
    my ($refseq) = @_;
    my @lines = `echo "$refseq" | svr_aliases_to_pegs`;
    my ($peg) = map { chomp; [split/\t/]->[1] } @lines;
    return $peg;
}

sub refseq_to_ncbi_gene {
    my ($refseq) = @_;
    my $korg = $refseq =~ /^YPA/ ? 'ypa' : 'ype';
    my $url = "http://www.genome.jp/dbget-bin/www_bget?$korg:$refseq";
    my $htm = wget($url);
    ($url) = $htm =~ m|(http://www\.ncbi\.nlm\.nih\.gov/gene/\d+)|;
    # print STDERR "url = $url\n";
    $htm = wget($url);
    my ($ctg, $beg, $end, $complement) = $htm =~ m|"magin_top_m1">(\S+) \((\d+)\.\.(\d+)([,)])|;
    # print STDERR "complement = '$complement'\n";
    if ($complement eq ',') { ($beg, $end) = ($end, $beg); }
    return ($ctg, $beg, $end);
}

sub wget {
    my ($url) = @_;
    # my $htm = `wget --quiet -e robots=off -U firefox --wait 1  -O /dev/stdout '$url'`;
    my $htm = `wget --quiet -e robots=off -U 'Safari 5.1.7 -- Windows' --wait 1  -O /dev/stdout '$url'`;
    sleep 1;
    return $htm;
}
