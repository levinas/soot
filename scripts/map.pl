#! /usr/bin/env perl

use strict;
use Carp;
use Cwd 'abs_path';
use Data::Dumper;
use Getopt::Long;

my $usage = <<"End_of_Usage";

usage: $0 [options] ref.fq reads_1.fq [reads_2.fq]

       -a algo          - alignment algorithm [bwa_mem bwa_sampe bowtie2 bbmap mosaik last] (D = bwa_mem)
       -o dir           - output directory (D = ref_reads_[algo])
       -t int           - number of threads (D = 8)
       -m size          - max memory per thread; suffix K/M/G recognized (D = 2G)
       --vc tool        - variant calling [samtools freebayes] (D = freebayes)

End_of_Usage

my ($help, $algo, $memory, $nthread, $outdir, $paired, $vc);

GetOptions("h|help"        => \$help,
           "a|algo=s"      => \$algo,
           "m|memory=s"    => \$memory,
           "o|outdir=s"    => \$outdir,
           "p|paired"      => \$paired,
           "t|threads=i"   => \$nthread,
           "vc=s"          => \$vc);

my $ref   = shift @ARGV;
my $read1 = shift @ARGV;
my $read2 = shift @ARGV;

$ref && $read1 or die $usage;

$ref   = abs_path($ref);
$read1 = abs_path($read1);
$read2 = abs_path($read2) if $read2;
$read2 = other_read_file_in_pair($read1) if $paired;

$nthread ||= 8;
$memory  ||= '2G'; $memory .= 'G' if $memory =~ /\d$/;
$algo    ||= 'bwa_mem'; $algo .= "_se" if !$read2;
$vc      ||= 'freebayes';
$outdir  ||= generate_dir_name($algo, $ref, $read1);

print "READS = $read1 $read2\n";
print "ALGO = $algo\n";

if (eval "defined(&map_with_$algo)") {
    print "> $outdir\n";
    run("mkdir -p $outdir");
    chdir($outdir);
    eval "&map_with_$algo";
    if ($@) {
        print $@;
    } else {
        compute_stats();
        if ($vc && eval "defined(&call_variant_with_$vc)") {
            eval "&call_variant_with_$vc";
            print $@ if $@;
            compute_consensus() unless $@;
        }
    }
} else {
    die "Mapping algorithm not defined: $algo\n";
}

sub generate_dir_name {
    my ($algo, $ref, $reads) = @_;
    $ref   =~ s|.*/||; $ref   =~ s/\.(fasta|fna|fa)//;
    $reads =~ s|.*/||; $reads =~ s/\.(fastq|fq).*//; $reads =~ s/_(1|2)//;
    return "$ref\_$reads\_$algo";
}

sub other_read_file_in_pair {
    my ($r1) = @_;
    my $r2 = $r1;
    $r2 =~ s/R1\./R2\./;
    return $r2 if -s $r2 && $r2 ne $r1;
}

sub call_variant_with_samtools {
    -s "mpileup"        or run("samtools mpileup -6 -uf ref.fa aln.bam > mpileup");
    -s "var.sam.vcf"    or run("bcftools call -vc mpileup > var.sam.vcf");
    -s "var.sam.count"  or run("grep -v '^#' var.sam.vcf |cut -f4 |grep -v 'N' |wc -l > var.sam.count");
    run("ln -s -f var.sam.vcf var.vcf");
}

sub call_variant_with_freebayes {
    -s "var.fb.vcf"     or run("bash -c 'freebayes-parallel <(fasta_generate_regions.py ref.fa.fai 100000) $nthread -p 1 -f ref.fa aln.bam >var.fb.vcf'");
    -s "var.fb.q10.vcf" or run("vcffilter -f 'QUAL > 10 & DP > 5' var.fb.vcf > var.fb.q10.vcf");
    -s "var.fb.q1.vcf"  or run("vcffilter -f 'QUAL > 1' var.fb.vcf > var.fb.q1.vcf");
    -s "var.fb.count"   or run("grep -v '^#' var.fb.vcf |cut -f4 |grep -v 'N' |wc -l > var.fb.count");
    run("ln -s -f var.fb.q10.vcf var.vcf");
}

sub compute_consensus {
    -s "var.vcf.gz"     or run("bgzip -c var.vcf > var.vcf.gz");
    -s "var.vcf.gz.tbi" or run("tabix var.vcf.gz");
    -s "consensus"      or run("bcftools consensus -c chain -f ref.fa var.vcf.gz >consensus");
}

sub compute_stats {
    -s "ref.fa.fai"     or run("samtools faidx ref.fa");
    -s "raw.flagstat"   or run("samtools flagstat aln.raw.sam > raw.flagstat");
    -s "flagstat"       or run("samtools flagstat aln.bam > flagstat");
    -s "stats"          or run("samtools stats aln.bam -c 1,8000,1 > stats");
  # -s "depth"          or run("samtools depth aln.bam > depth");
    -s "depth"          or run("bedtools genomecov -ibam aln.bam -d > depth");
    -s "depth.hist"     or run("bedtools genomecov -ibam aln.bam > depth.hist");
    -s "uncov.10"       or run("bedtools genomecov -ibam aln.bam -bga | perl -ne 'chomp; \@c=split/\t/; \$ln=\$c[2]-\$c[1]; print join(\"\\t\", \@c, \$ln).\"\\n\" if \$c[3]<10;' > uncov.10" );
    # BED start position 0-based and the end position 1-based (Example: NC_000962,1987085,1987701,0,616; the 0 coverage base really starts at 1987086)
}

sub map_with_bwa_mem {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read_1.fq"      or run("ln -s $read1 read_1.fq");
    -s "read_2.fq"      or run("ln -s $read2 read_2.fq");
    -s "ref.fa.bwt"     or run("bwa index ref.fa");
    -s "aln-pe.sam"     or run("bwa mem -t $nthread ref.fa read_1.fq read_2.fq > aln-pe.sam 2>mem.log");
    -s "aln.raw.sam"    or run("ln -s aln-pe.sam aln.raw.sam");
    -s "aln.keep.bam"   or run("samtools view -@ $nthread -f 0x2 -bS aln.raw.sam > aln.keep.bam"); # keep only properly paired reads
    -s "unmapped.bam"   or run("samtools view -@ $nthread -f 4 -bS aln.raw.sam > unmapped.bam");
    -s "aln.sorted.bam" or run("samtools sort -m $memory -@ $nthread aln.keep.bam aln.sorted");
  # -s "aln.dedup.bam"  or run("samtools rmdup aln.sorted.bam aln.dedup.bam");  # rmdup broken in samtools v1.0 and v1.1
  # -s "aln.bam"        or run("ln -s aln.dedup.bam aln.bam");
    -s "aln.bam"        or run("ln -s aln.sorted.bam aln.bam");
    -s "aln.bam.bai"    or run("samtools index aln.bam");
}

sub map_with_bwa_mem_se {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read.fq"        or run("ln -s $read1 read.fq");
    -s "ref.fa.bwt"     or run("bwa index ref.fa");
    -s "aln-se.sam"     or run("bwa mem -t $nthread ref.fa read.fq > aln-se.sam 2>mem.log");
    -s "aln.raw.sam"    or run("ln -s aln-se.sam aln.raw.sam");
    -s "aln.keep.bam"   or run("samtools view -@ $nthread -bS aln.raw.sam > aln.keep.bam");
    -s "unmapped.bam"   or run("samtools view -@ $nthread -f 4 -bS aln.raw.sam > unmapped.bam");
    -s "aln.sorted.bam" or run("samtools sort -m $memory -@ $nthread aln.keep.bam aln.sorted");
  # -s "aln.dedup.bam"  or run("samtools rmdup aln.sorted.bam aln.dedup.bam");  # rmdup broken in samtools v1.0 and v1.1
  # -s "aln.bam"        or run("ln -s aln.dedup.bam aln.bam");
    -s "aln.bam"        or run("ln -s aln.sorted.bam aln.bam");
    -s "aln.bam.bai"    or run("samtools index aln.bam");
}

sub map_with_bwa_mem_strict {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read_1.fq"      or run("ln -s $read1 read_1.fq");
    -s "read_2.fq"      or run("ln -s $read2 read_2.fq");
    -s "ref.fa.bwt"     or run("bwa index ref.fa");
    -s "aln-pe.sam"     or run("bwa mem -B9 -O16 -E1 -L5 -t $nthread ref.fa read_1.fq read_2.fq > aln-pe.sam 2>mem.log");
    -s "aln.raw.sam"    or run("ln -s aln-pe.sam aln.raw.sam");
    -s "aln.keep.bam"   or run("samtools view -@ $nthread -f 0x2 -q 10 -bS aln.raw.sam > aln.keep.bam"); # keep only properly paired reads
    -s "unmapped.bam"   or run("samtools view -@ $nthread -f 4 -bS aln.raw.sam > unmapped.bam");
    -s "aln.sorted.bam" or run("samtools sort -m $memory -@ $nthread aln.keep.bam aln.sorted");
  # -s "aln.dedup.bam"  or run("samtools rmdup aln.sorted.bam aln.dedup.bam");  # rmdup broken in samtools v1.0 and v1.1
  # -s "aln.bam"        or run("ln -s aln.dedup.bam aln.bam");
    -s "aln.bam"        or run("ln -s aln.sorted.bam aln.bam");
    -s "aln.bam.bai"    or run("samtools index aln.bam");
}

sub map_with_bwa_mem_strict_se {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read.fq"        or run("ln -s $read1 read.fq");
    -s "ref.fa.bwt"     or run("bwa index ref.fa");
    -s "aln-se.sam"     or run("bwa mem -B9 -O16 -E1 -L5 -t $nthread ref.fa read.fq > aln-se.sam 2>mem.log");
    -s "aln.raw.sam"    or run("ln -s aln-se.sam aln.raw.sam");
    -s "aln.keep.bam"   or run("samtools view -@ $nthread -bS -q 10 aln.raw.sam > aln.keep.bam");
    -s "unmapped.bam"   or run("samtools view -@ $nthread -f 4 -bS aln.raw.sam > unmapped.bam");
    -s "aln.sorted.bam" or run("samtools sort -m $memory -@ $nthread aln.keep.bam aln.sorted");
  # -s "aln.dedup.bam"  or run("samtools rmdup aln.sorted.bam aln.dedup.bam");  # rmdup broken in samtools v1.0 and v1.1
  # -s "aln.bam"        or run("ln -s aln.dedup.bam aln.bam");
    -s "aln.bam"        or run("ln -s aln.sorted.bam aln.bam");
    -s "aln.bam.bai"    or run("samtools index aln.bam");
}

sub map_with_last {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read_1.fq"      or run("ln -s $read1 read_1.fq");
    -s "read_2.fq"      or run("ln -s $read2 read_2.fq");
    -s "index.suf"      or run("lastdb -m1111110 index ref.fa");
    -s "out1.maf"       or run("parallel-fastq -j $nthread -k 'lastal -Q1 -d108 -e120 -i1 index' < read_1.fq > out1.maf");
    -s "out2.maf"       or run("parallel-fastq -j $nthread -k 'lastal -Q1 -d108 -e120 -i1 index' < read_2.fq > out2.maf");
  # -s "out1.maf"       or run("lastal -Q1 -d108 -e120 -i1 index read_1.fq > out1.maf"); # sequential
  # -s "out2.maf"       or run("lastal -Q1 -d108 -e120 -i1 index read_2.fq > out2.maf"); # sequential
    -s "aln-pe.maf"     or run("last-pair-probs -m 0.1 out1.maf out2.maf > aln-pe.maf");
    -s "ref.fa.fai"     or run("samtools faidx ref.fa");
    -s "sam.header"     or run("awk '{ print \"\@SQ\\tSN:\"\$1\"\\tLN:\"\$2 }' ref.fa.fai > sam.header");
    -s "aln.raw.sam"    or run("bash -c 'cat sam.header <(maf-convert sam aln-pe.maf) > aln.raw.sam'");
    -s "aln.keep.bam"   or run("samtools view -@ $nthread -bS aln.raw.sam > aln.keep.bam");
    -s "unmapped.bam"   or run("samtools view -@ $nthread -f 4 -bS aln.raw.sam > unmapped.bam");
    -s "aln.sorted.bam" or run("samtools sort -m $memory -@ $nthread aln.keep.bam aln.sorted");
  # -s "aln.dedup.bam"  or run("samtools rmdup aln.sorted.bam aln.dedup.bam");  # rmdup broken in samtools v1.0 and v1.1
  # -s "aln.bam"        or run("ln -s aln.dedup.bam aln.bam");
    -s "aln.bam"        or run("ln -s aln.sorted.bam aln.bam");
    -s "aln.bam.bai"    or run("samtools index aln.bam");
}

sub map_with_last_se {
    -s "ref.fa"         or run("ln -s $ref ref.fa");
    -s "read_1.fq"      or run("ln -s $read1 read.fq");
    -s "index.suf"      or run("lastdb -m1111110 index ref.fa");
    -s "out.maf"        or run("parallel-fastq -j $nthread -k 'lastal -Q1 -d108 -e120 -i1 index' < read.fq > out.maf");
  # -s "out.maf"        or run("lastal -Q1 -d108 -e120 -i1 index read.fq > out.maf"); # sequential
    -s "aln-se.maf"     or run("ln -s out.maf aln-se.maf");
    -s "ref.fa.fai"     or run("samtools faidx ref.fa");
    -s "sam.header"     or run("awk '{ print \"\@SQ\\tSN:\"\$1\"\\tLN:\"\$2 }' ref.fa.fai > sam.header");
    -s "aln.raw.sam"    or run("bash -c 'cat sam.header <(maf-convert sam aln-se.maf) > aln.raw.sam'");
    -s "aln.keep.bam"   or run("samtools view -@ $nthread -bS aln.raw.sam > aln.keep.bam");
    -s "unmapped.bam"   or run("samtools view -@ $nthread -f 4 -bS aln.raw.sam > unmapped.bam");
    -s "aln.sorted.bam" or run("samtools sort -m $memory -@ $nthread aln.keep.bam aln.sorted");
  # -s "aln.dedup.bam"  or run("samtools rmdup aln.sorted.bam aln.dedup.bam");  # rmdup broken in samtools v1.0 and v1.1
  # -s "aln.bam"        or run("ln -s aln.dedup.bam aln.bam");
    -s "aln.bam"        or run("ln -s aln.sorted.bam aln.bam");
    -s "aln.bam.bai"    or run("samtools index aln.bam");
}

sub run { system($_[0]) == 0 or confess("FAILED: $_[0]"); }
