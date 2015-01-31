#! /usr/bin/env perl

use strict;
use Carp;

my $usage = "Usage: $0 ref reads output_dir mapq\n\n";

my $ref   = shift @ARGV;
my $reads = shift @ARGV;
my $dir   = shift @ARGV;
my $mapq  = shift @ARGV || 30;

$ref && $reads && $dir or die $usage;

run("mkdir -p $dir");

-s "$ref.bwt"           or run("bwa index -a bwtsw $ref");
-s "$dir/aln.sai"       or run("date; time bwa aln -t 8 $ref $reads > $dir/aln.sai; date");
-s "$dir/aln.sam"       or run("bwa samse $ref $dir/aln.sai $reads > $dir/aln.sam");
-s "$dir/aln.bam"       or run("samtools view -bS -q $mapq $dir/aln.sam > $dir/aln.bam");
-s "$dir/aln.sorted"    or run("samtools sort $dir/aln.bam $dir/aln.sorted");
-s "$dir/aln.rmdup"     or run("samtools rmdup -s $dir/aln.sorted.bam $dir/aln.rmdup");
-s "$dir/aln.rmdup.bai" or run("samtools index $dir/aln.rmdup");
-s "$dir/aln.mpileup"   or run("samtools mpileup -f $ref $dir/aln.rmdup > $dir/aln.mpileup");
-s "$dir/var.raw.vcf"   or run("samtools mpileup -uf $ref $dir/aln.rmdup | bcftools view -vcg - > $dir/var.raw.vcf");
-s "$dir/var.AF1.vcf"   or run("grep -v 'AF1=0' $dir/var.raw.vcf > $dir/var.AF1.vcf");
-s "$dir/var.count"     or run("grep -v '^#' $dir/var.AF1.vcf |wc -l > $dir/var.count");

sub run { system($_[0]) == 0 or confess("FAILED: $_[0]"); }
