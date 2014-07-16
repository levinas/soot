#! /usr/bin/env perl

use strict;
use Carp;
use Data::Dumper;
use Getopt::Long;

my $usage = <<"End_of_Usage";

usage: $0 [options] read_dir ref_fasta

       -a algo          - alignment algorithm [bwa_sampe bwa_mem bowtie2 bbmap nasp] (D = bwa_mem)
       -o dir           - output directory (D = .)
       -t int           - number of threads (D = 8)
  
  read_dir structure:
       project_ID_1 / SRR1_1/2.fastq
                     SRR2_1/2.fastq
       project_ID_2 / SRR3_1/2.fastq
                     SRR4_1/2.fastq
       ...

End_of_Usage

my ($help, $algo, $nthread, $outdir);

GetOptions("h|help"        => \$help,
           "a|algo=s"      => \$algo,
           "o|outdir=s"    => \$outdir,
           "t|threads=i"   => \$nthread);

my $read_dir  = shift @ARGV;
my $ref_fasta = shift @ARGV;

$ref_fasta && $read_dir or die $usage;

$nthread ||= 8;
$algo    ||= 'bwa_mem';
$outdir  ||= '.';

my @proj_dirs = map { chomp; $_ } `ls -d $read_dir/*`;

for my $subd (@proj_dirs) {
    my $proj = $subd; $proj =~ s/.*\///;

    print STDERR join(" ", '>', $proj, $subd) . "\n";

    my @pairs = get_read_pairs($subd);

    for (@pairs) {
        my ($f1, $f2) = @$_;
        my $cmd = "snp.pl -t $nthread -a $algo $ref_fasta $f1 $f2";
        print "$cmd\n";
        run($cmd);
    }
}


sub get_read_pairs {
    my ($dir) = @_;
    my @pairs = map { /(\S+)(_R|_)1(_|\.)(\S+)/ ? [ "$dir/$1$2"."1"."$3$4", "$dir/$1$2"."2"."$3$4" ] : () }
                grep { /(fastq|fq)(|\.tar\.gz|\.tgz|\.tar.bz|\.gz|\.bz)$/ } map { chomp; $_ } `ls $dir`;
    wantarray ? @pairs : \@pairs;
}

sub run { system(@_) == 0 or confess("FAILED: ". join(" ", @_)); }
