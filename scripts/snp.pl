#! /usr/bin/env perl

use strict;
use Carp;
use Cwd 'abs_path';
use Data::Dumper;
use Getopt::Long;

my $usage = <<"End_of_Usage";

usage: snp.pl [options] ref.fq reads_1.fq reads_2.fq 

       -a algo          - alignment algorithm [bwa_sampe bwa_mem bowtie2 bbmap nasp] (D = bwa_mem)
       -o dir           - output directory (D = ref_reads_[algo])
       -t int           - number of threads (D = 8)

End_of_Usage

my ($help, $algo, $nthread, $outdir);

GetOptions("h|help"        => \$help,
           "a|algo=s"      => \$algo,
           "o|outdir=s"    => \$outdir,
           "t|threads=i"   => \$nthread);

my $ref   = shift @ARGV;
my $read1 = shift @ARGV;
my $read2 = shift @ARGV;

$ref && $read1 && $read2 or die $usage;

$ref   = abs_path($ref);
$read1 = abs_path($read1);
$read2 = abs_path($read2);

$nthread ||= 8;
$algo    ||= 'bwa_mem';
$outdir  ||= generate_dir_name($algo, $ref, $read1);

if (eval "defined(&map_with_$algo)") {
    print "> $outdir\n";
    run("mkdir -p $outdir");
    chdir($outdir);
    eval "&map_with_$algo";
    print $@ if $@;
} else {
    die "Mapping algorithm not defined: $algo\n";
}

sub generate_dir_name {
    my ($algo, $ref, $reads) = @_;
    $ref   =~ s|.*/||; $ref   =~ s/\.(fasta|fna|fa)//;
    $reads =~ s|.*/||; $reads =~ s/\.(fastq|fq).*//; $reads =~ s/_(1|2)//;
    return "$ref\_$reads\_$algo";
}

sub map_with_bwa_mem {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read_1.fq"      or run("ln -s $read1 read_1.fq");
    -s "read_2.fq"      or run("ln -s $read2 read_2.fq");
    -s "ref.fa.bwt"     or run("bwa index ref.fa");
    -s "aln-pe.sam"     or run("bwa mem -t $nthread ref.fa read_1.fq read_2.fq > aln-pe.sam 2>mem.log");
    -s "aln.bam"        or run("samtools view -bS aln-pe.sam > aln.bam");
    -s "aln.sorted.bam" or run("samtools sort aln.bam aln.sorted");
    -s "mpileup"        or run("samtools mpileup -uf ref.fa aln.sorted.bam > mpileup");
    -s "var.raw.vcf"    or run("cat mpileup| bcftools view -vcg - > var.raw.vcf");
    -s "var.count.all"  or run("grep -v '^#' var.raw.vcf |wc -l > var.count.all");
    -s "var.count"      or run("grep -v '^#' var.raw.vcf |cut -f4 |grep -v 'N' |wc -l > var.count");
}

sub map_with_bwa_sampe {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read_1.fq"      or run("ln -s $read1 read_1.fq");
    -s "read_2.fq"      or run("ln -s $read2 read_2.fq");
    -s "ref.fa.bwt"     or run("bwa0.6.2 index ref.fa");
    -s "aln_1.sai"      or run("bwa0.6.2 aln -t $nthread ref.fa read_1.fq > aln_1.sai");
    -s "aln_2.sai"      or run("bwa0.6.2 aln -t $nthread ref.fa read_2.fq > aln_2.sai");
    -s "aln.sam"        or run("bwa0.6.2 sampe ref.fa aln_1.sai aln_2.sai read_1.fq read_2.fq > aln.sam");
    -s "aln.bam"        or run("samtools view -bS aln.sam > aln.bam");
    -s "aln.sorted.bam" or run("samtools sort aln.bam aln.sorted");
    -s "mpileup"        or run("samtools mpileup -uf ref.fa aln.sorted.bam > mpileup");
    -s "var.raw.vcf"    or run("cat mpileup| bcftools view -vcg - > var.raw.vcf");
    -s "var.count.all"  or run("grep -v '^#' var.raw.vcf |wc -l > var.count.all");
    -s "var.count"      or run("grep -v '^#' var.raw.vcf |cut -f4 |grep -v 'N' |wc -l > var.count");
}

sub map_with_bowtie2 {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read_1.fq"      or run("ln -s $read1 read_1.fq");
    -s "read_2.fq"      or run("ln -s $read2 read_2.fq");
    -s "ref.X.bt2"      or run("bowtie2-build ref.fa -f ref");
    -s "aln.sam"        or run("bowtie2 --threads $nthread -x ref -1 read_1.fq -2 read_2.fq -S aln.sam");
    -s "aln.bam"        or run("samtools view -bS aln.sam > aln.bam");
    -s "aln.sorted.bam" or run("samtools sort aln.bam aln.sorted");
    -s "mpileup"        or run("samtools mpileup -uf ref.fa aln.sorted.bam > mpileup");
    -s "var.raw.vcf"    or run("cat mpileup| bcftools view -vcg - > var.raw.vcf");
    -s "var.count.all"  or run("grep -v '^#' var.raw.vcf |wc -l > var.count.all");
    -s "var.count"      or run("grep -v '^#' var.raw.vcf |cut -f4 |grep -v 'N' |wc -l > var.count");
}

sub map_with_bbmap {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read_1.fq"      or run("ln -s $read1 read_1.fq");
    -s "read_2.fq"      or run("ln -s $read2 read_2.fq");
    -s "ref/index"      or run("bbmap.sh ref=ref.fa");
    -s "aln.sam"        or run("bbmap.sh threads=$nthread in=read_1.fq in2=read_2.fq out=aln.sam");
    -s "aln.bam"        or run("samtools view -bS aln.sam > aln.bam");
    -s "aln.sorted.bam" or run("samtools sort aln.bam aln.sorted");
    -s "mpileup"        or run("samtools mpileup -uf ref.fa aln.sorted.bam > mpileup");
    -s "var.raw.vcf"    or run("cat mpileup| bcftools view -vcg - > var.raw.vcf");
    -s "var.count.all"  or run("grep -v '^#' var.raw.vcf |wc -l > var.count.all");
    -s "var.count"      or run("grep -v '^#' var.raw.vcf |cut -f4 |grep -v 'N' |wc -l > var.count");
}


sub map_with_nasp {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read_1.fq"      or run("ln -s $read1 read_1.fq");
    -s "read_2.fq"      or run("ln -s $read2 read_2.fq");
    run("nasp ref.fa");
    # -s "ref/index"      or run("bbmap.sh ref=ref.fa");
    # -s "aln.sam"        or run("bbmap.sh threads=$nthread in=read_1.fq in2=read_2.fq out=aln.sam");
    # -s "aln.bam"        or run("samtools view -bS aln.sam > aln.bam");
    # -s "aln.sorted.bam" or run("samtools sort aln.bam aln.sorted");
    # -s "mpileup"        or run("samtools mpileup -uf ref.fa aln.sorted.bam > mpileup");
    # -s "var.raw.vcf"    or run("cat mpileup| bcftools view -vcg - > var.raw.vcf");
    # -s "var.count.all"  or run("grep -v '^#' var.raw.vcf |wc -l > var.count.all");
    # -s "var.count"      or run("grep -v '^#' var.raw.vcf |cut -f4 |grep -v 'N' |wc -l > var.count");
}


sub run { system($_[0]) == 0 or confess("FAILED: $_[0]"); }
