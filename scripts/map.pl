#! /usr/bin/env perl

use strict;
use Carp;
use Cwd 'abs_path';
use Data::Dumper;
use Getopt::Long;

my $usage = <<"End_of_Usage";

usage: $0 [options] ref.fq reads_1.fq reads_2.fq

       -a algo          - alignment algorithm [bwa_mem bwa_sampe bowtie2 bbmap mosaik] (D = bwa_mem)
       -o dir           - output directory (D = ref_reads_[algo])
       -t int           - number of threads (D = 8)
       -m size          - max memory per thread; suffix K/M/G recognized (D = 2G)
       --vc             - variant calling

End_of_Usage

my ($help, $algo, $memory, $nthread, $outdir, $vc);

GetOptions("h|help"        => \$help,
           "a|algo=s"      => \$algo,
           "m|memory=s"    => \$memory,
           "o|outdir=s"    => \$outdir,
           "t|threads=i"   => \$nthread,
           "vc"            => \$vc);

my $ref   = shift @ARGV;
my $read1 = shift @ARGV;
my $read2 = shift @ARGV;

$ref && $read1 && $read2 or die $usage;

$ref   = abs_path($ref);
$read1 = abs_path($read1);
$read2 = abs_path($read2);

$nthread ||= 8;
$memory  ||= '2G';
$algo    ||= 'bwa_mem';
$outdir  ||= generate_dir_name($algo, $ref, $read1);

if (eval "defined(&map_with_$algo)") {
    print "> $outdir\n";
    run("mkdir -p $outdir");
    chdir($outdir);
    eval "&map_with_$algo";
    print $@ if $@;
    call_variant_with_bcftools() if $vc && !$@;
} else {
    die "Mapping algorithm not defined: $algo\n";
}

sub generate_dir_name {
    my ($algo, $ref, $reads) = @_;
    $ref   =~ s|.*/||; $ref   =~ s/\.(fasta|fna|fa)//;
    $reads =~ s|.*/||; $reads =~ s/\.(fastq|fq).*//; $reads =~ s/_(1|2)//;
    return "$ref\_$reads\_$algo";
}

sub call_variant_with_bcftools {
    -s "mpileup"        or run("samtools mpileup -6 -uf ref.fa aln.bam > mpileup");
    -s "var.raw.vcf"    or run("bcftools view -vcg mpileup > var.raw.vcf");
    -s "var.count.all"  or run("grep -v '^#' var.raw.vcf |wc -l > var.count.all");
    -s "var.count"      or run("grep -v '^#' var.raw.vcf |cut -f4 |grep -v 'N' |wc -l > var.count");
}

sub map_with_bwa_mem {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read_1.fq"      or run("ln -s $read1 read_1.fq");
    -s "read_2.fq"      or run("ln -s $read2 read_2.fq");
    -s "ref.fa.bwt"     or run("bwa index ref.fa");
    -s "aln-pe.sam"     or run("bwa mem -t $nthread ref.fa read_1.fq read_2.fq > aln-pe.sam 2>mem.log");
    -s "aln-pe.bam"     or run("samtools view -@ $nthread -f 0x2 -bS aln-pe.sam > aln-pe.bam"); # keep only properly paired reads
    -s "unmapped.bam"   or run("samtools view -@ $nthread -f 4 -bS aln-pe.sam > unmapped.bam");
    -s "aln.sorted.bam" or run("samtools sort -m $memory -@ $nthread aln-pe.bam aln.sorted");
  # -s "aln.dedup.bam"  or run("samtools rmdup aln.sorted.bam aln.dedup.bam");  # rmdup broken in samtools v1.0 and v1.1
  # -s "aln.bam"        or run("ln -s aln.dedup.bam aln.bam");
    -s "aln.bam"        or run("ln -s aln.sorted.bam aln.bam");
    -s "depth"          or run("samtools depth aln.bam > depth");
    -s "raw.flagstat"   or run("samtools flagstat aln-pe.sam > raw.flagstat");
    -s "flagstat"       or run("samtools flagstat aln.bam > flagstat");
    -s "stats"          or run("samtools stats aln.bam -c 1,8000,1 > stats");
}


sub run { system($_[0]) == 0 or confess("FAILED: $_[0]"); }
