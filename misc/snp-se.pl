#! /usr/bin/env perl

use strict;
use Carp;

my $usage = "Usage: $0 ref reads output_dir\n\n";

my $ref   = shift @ARGV;
my $reads = shift @ARGV;
my $dir   = shift @ARGV;

$ref && $reads && $dir or die $usage;

run("mkdir -p $dir");

-s "$ref.bwt"         or run("bwa index -a bwtsw $ref");
-s "$dir/aln.sai"     or run("date; time bwa aln -t 8 $ref $reads > $dir/aln.sai; date");
-s "$dir/aln.sam"     or run("bwa samse -t 8 $ref $dir/aln.sai $reads > $dir/aln.sam");
-s "$dir/aln.bam"     or run("samtools view -bS $dir/aln.sam > $dir/aln.bam");
-s "$dir/aln.sorted"  or run("samtools sort $dir/aln.bam $dir/aln.sorted");
-s "$dir/var.raw.vcf" or run("samtools mpileup -uf $ref $dir/aln.sorted.bam | bcftools view -vcg - > $dir/var.raw.vcf");
-s "$dir/var.count"   or run("grep -v '^#' $dir/var.raw.vcf |wc -l > $dir/var.count");

sub run { system($_[0]) == 0 or confess("FAILED: $_[0]"); }
