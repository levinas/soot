#! /usr/bin/env perl

use strict;
use Data::Dumper;
use SAPserver;
use SeedUtils;

use Storable qw(nstore retrieve);
use File::Basename;
use lib dirname (__FILE__);
use DT;

my $usage = "Usage: $0 genome-id < var.vcf > var.table\n\n";

my $gid = shift @ARGV or die $usage;
use Getopt::Long;

my ($help, $min_alt_depth, $min_alt_fract, $min_map_quality, $show_html, @show_gids, $known_snp_file);

GetOptions("h|help"     => \$help,
           "d=i"        => \$min_alt_depth,
           "f=f"        => \$min_alt_fract,
           "q=f"        => \$min_map_quality,
           "html"       => \$show_html,
           "gids=s{1,}" => \@show_gids,
           "known=s"    => \$known_snp_file,
	  ) or die("Error in command line arguments\n");

$help and die $usage;

$min_alt_depth   ||= 5;
$min_alt_fract   ||= 0.95;
$min_map_quality ||= 30;


my $vars = read_var();
my $features = get_sorted_features_for_org($gid);

# nstore($vars, 'var.store');
# nstore($features, 'features.store');
# print STDERR '$features = '. Dumper($features);
# exit;

# my $vars = retrieve('var.store');
# my $features = retrieve('features.store');

my @snps;
for my $var (@$vars) {
    my ($ctg, $pos, $ref, $alt, $info) = @{$var}[0, 1, 3, 4, 7];
    my @dp4 = split(/,/, $info->{DP4});
    my $alt_dp = $dp4[2] + $dp4[3];
    my $alt_frac =  sprintf("%.2f", $alt_dp / sum(@dp4));
    my $map_qual = $info->{MQ};

    print STDERR join("\t", $ctg, $pos, $alt_dp, $alt_frac, $map_qual) . "\n";
    
    next unless $alt_frac >= $min_alt_fract && $alt_dp >= $min_alt_depth && $map_qual >= $min_map_quality;

    my $type = length($alt) > length($ref) ? 'Insertion' : length($alt) < length($ref) ? 'Deletion' : undef;
    my $hash = feature_info_for_position($ctg, $pos, $features);
    my ($nt1, $nt2, $aa1, $aa2);
    my ($locus, $gene_name, $func);
    my $gene = $hash->{gene};
    if ($gene) {
        my $dna = lc $gene->[10];
        my $strand = $gene->[5];
        my $p = $hash->{pos_in_gene};
        my $win = SeedUtils::max(length($ref), length($alt));
        # my $win = length($ref);
        my $beg = $p - 1; $beg -= (length($ref)-1) if $strand eq '-';
        my $end = $beg + $win - 1;
        while ($beg % 3) { $beg-- }
        while (($end+1) % 3) { $end++ }
        my $size = $end - $beg + 1;
        my $alt_dna = $alt; $alt_dna = SeedUtils::rev_comp($alt_dna) if $strand eq '-';
        
        my $beg_delta = $strand eq '+' ? ($p-$beg-1) : ($p-$beg-length($ref));
        $nt1 = substr($dna, $beg, $size);
        $nt2 = $nt1; substr($nt2, $beg_delta, length($ref)) = $alt_dna;
        $aa1 = SeedUtils::translate($nt1, undef, 0);
        $aa2 = SeedUtils::translate($nt2, undef, 0);
        if (length($alt) == length($ref)) {
            $type = $aa1 eq $aa2 ? 'Synon' : 'Nonsyn';
        }
        $locus = $gene->[9]->{LocusTag};
        $locus = $locus->[0] if $locus;
        $locus =~ s/\S+://;
        $gene_name = $gene->[9]->{GENE}->[0];
        $gene_name =~ s/\S+://;

        $func = $gene->[8]; $func ||= 'hypothetical protein' if $gene->[0];
        # $func = "$gene_name, $func" if $gene_name;

        # if ($pos == 3459494) {
        #     print STDERR "<<\n";
        #     print join("\t", $gene->[7], $p, $win) . "\n";
        #     print join("\t", $beg, $end, $size) . "\n";
        #     print join("\t", $beg_delta, $p-1) . "\n";
        #     print join("\t", $nt1, $nt2) . "\n";
        #     print join("\t", $aa1, $aa2) . "\n";
        #     print STDERR ">>\n";
        # }
        # print STDERR join("\t", $ctg, $pos, $gene->[0], $func)."\n";
        # print STDERR 'alias = '. Dumper($gene->[9]);
        # print STDERR "\n";
        
    }

    push @snps, [ $ctg, $pos, $ref, $alt, $alt_dp, $alt_frac,
                  $type, $nt1, $nt2, $aa1, $aa2, 
                  $locus, [ $gene->[0], $func ],
                  $hash->{frameshift_region} ? 'frameshift' : undef,
                  [ @{$hash->{left}}[0, 8] ],
                  [ @{$hash->{right}}[0, 8] ]
                ];
    
    # print join("\t", $ctg, $pos, $ref, $alt, $alt_dp, $alt_frac, $nt1, $nt2, $aa1, $aa2, $hash->{frameshift_region} ? 'frameshift' : undef) . "\n";
    # for my $k (qw(gene left right)) {
        # my $v = $hash->{$k} or next;
        # my $func = $v->[8];
        # print join("\t", $k, $v->[0], $v->[6], $func) . "\n";
    # }
    # print '$gene = '. Dumper($hash);
    # print "\n=========\n";
}

my %known_snps;
if ($known_snp_file && -s $known_snp_file) {
    %known_snps = map { my ($ctg, $pos) = split/\t/;
                        $ctg && $pos ? ("$ctg,$pos" => 1) : ()
                    } `cat $known_snp_file`;
}
if ($show_html) {
    my @head = ('Contig', 'Pos', 'Ref', 'Var', 'Var cov', 'Var frac',
                'Type', 'Ref nt', 'Var nt', 'Ref aa', 'Var aa', 
                'Gene ID', 'Gene', 'Ref frame',
                "Neighboring feature (3' end)",
                "Neighboring feature (5' end)" );
    my @rows;
    for (@snps) {
        my $known = $known_snps{"$_->[0],$_->[1]"} and next;
        my @c = map { DT::span_css($_, 'wrap') }
                map { $known ? DT::span_css($_, "opaque") : $_ }
                map { ref $_ eq 'ARRAY' ? $_->[0] ? linked_fid(@$_, \@show_gids) : undef : $_ } @$_;
        push @rows, \@c;
    }
    DT::print_dynamic_table(\@head, \@rows, { title => 'SNPs in BD vs CO92', extra_css => extra_css() });
} else {
    for (@snps) {
        next if $known_snps{"$_->[0],$_->[1]"};
        # my @c = map { ref $_ eq 'ARRAY' ? $_->[0] ? $_->[0] : undef : $_ } @$_;
        my @c = map { ref $_ eq 'ARRAY' ? $_->[0] ? $_->[1] : undef : $_ } @$_;
        print join("\t", @c) . "\n";
    }
}

sub extra_css {
    return <<End_of_CSS;
<style type="text/css">    
  .opaque {
      opacity:0.50;
  }
  .highlight {
    text-decoration:underline;
  }
</style>
End_of_CSS
}

sub linked_fid {
    my ($fid, $txt, $gids) = @_;
    $txt ||= $fid;
    my $url = "http://pubseed.theseed.org/?page=Annotation&feature=$fid";
    if ($gids) { $url .= "&show_genome=$_" for @$gids }
    return add_link($url, $txt);
}

sub add_link {
    my ($url, $txt) = @_;
    $txt ||= $url;
    return "<a href=$url>$txt</a>";
}

sub feature_info_for_position {
    my ($ctg, $pos, $features) = @_;
    
    my ($left, $right, @cover);

    my $index = binary_search_in_sorted_features($ctg, $pos, $features);

    $left = $features->{$ctg}->[$index] if defined($index) && $index >= 0;

    my @cover;
    for (my $i = $index + 1; $i < @{$features->{$ctg}}; $i++) {
        my $fea = $features->{$ctg}->[$i];
        my ($lo, $hi) = @{$fea}[3, 4];
        my $overlap = $lo <= $pos && $pos <= $hi;
        if (! $overlap) { $right = $fea; last; }
        push @cover, $fea;
    }

    my @overlapping_genes = sort { $b->[6] <=> $a->[6] } grep { $_->[0] =~ /peg/ } @cover; # sort genes by length
    my @overlapping_other = grep { $_->[0] !~ /peg/ } @cover;

    my $frameshift = 0;
    my (%func_seen, %func_cnt, %func_lo, %func_hi);
    for (@overlapping_genes) {
        my $len  = $_->[6];
        my $func = $_->[8];
        my $lo   = $_->[3];
        my $hi   = $_->[4];
        $frameshift = 1 if $len % 3;
        if (!hypo_or_mobile($func)) {
            $func_cnt{$func}++;
            $func_lo{$func} = $lo if $lo < $func_lo{$func} || ! defined($func_lo{$func});
            $func_hi{$func} = $hi if $hi > $func_hi{$func};
        }
    }
    for (($left, $right)) {
        my $func = $_->[8];
        my $lo   = $_->[3];
        my $hi   = $_->[4];
        if ($func_cnt{$func}) {
            $func_cnt{$func}++;
            $func_lo{$func} = $lo if $lo < $func_lo{$func} || ! defined($func_lo{$func});
            $func_hi{$func} = $hi if $hi > $func_hi{$func};
        }
    }
    for my $func (keys %func_cnt) {
        my $len = $func_hi{$func} - $func_lo{$func} + 1;
        $frameshift = 1 if $func_cnt{$func} > 1 && $len % 3;
    }
    
    my ($gene) = @overlapping_genes;
    my $pos_in_gene;
    if ($gene) {
        my ($lo, $hi, $strand) = @{$gene}[3, 4, 5];
        $pos_in_gene =  $strand eq '+' ? $pos - $lo + 1 : $hi - $pos + 1;
    }

    my %hash;

    $hash{overlapping_genes} = \@overlapping_genes if @overlapping_genes;
    $hash{overlapping_other} = \@overlapping_other if @overlapping_other;
    $hash{gene}              = $gene               if $gene;
    $hash{pos_in_gene}       = $pos_in_gene        if $pos_in_gene;
    $hash{frameshift_region} = $frameshift         if $frameshift;
    $hash{left}              = $left               if $left;
    $hash{right}             = $right              if $right;
    
    wantarray ? %hash : \%hash;
}

sub hypo_or_mobile {
    my ($func) = @_;
    return !$func || SeedUtils::hypo($func) || $func =~ /mobile/i;
}

sub feature_string {
    my ($feature) = @_;
    $feature ? '['. join(", ", @$feature[0,2,3,4]) .']' : undef;
}

# Assumes the features are in ascending order on left coordinate and then right coordinate.
# Find the index of the rightmost feature who does not have a right neighbor that is to the left of the position
# Return -1 if no such feature can be found.

sub binary_search_in_sorted_features {
    my ($ctg, $pos, $features, $x, $y) = @_;

    return    if !$features || !$features->{$ctg};
    return -1 if $features->{$ctg}->[0]->[4] >= $pos;

    my $feas = $features->{$ctg};
    my $n = @$feas;

    $x = 0      unless defined $x;
    $y = $n - 1 unless defined $y;

    while ($x < $y) {
        my $m = int(($x + $y) / 2);

        # Terminate if:
        #   (1) features[m] is to the left pos, and
        #   (2) features[m+1] covers or is to the right of pos

        my $m2 = $feas->[$m]->[4];
        my $n2 = $feas->[$m+1]->[4];

        return $m if $m2 < $pos && ($n2 >= $pos || !defined($n2));

        if ($m2 < $pos) { $x = $m + 1 } else { $y = $m }
    }
    
    return $x;
}

sub read_var {
    my ($file) = @_;
    my @vars = map  { my $hash = { map { my ($k, $v) = split /=/; $v ||= 1; $k => $v } split(/;/, $_->[7]) }; $_->[7] = $hash; $_ }
               map  { [split /\t/] }
               map  { chomp; $_ }
               grep { !/^#/ } $file ? `cat $file` : <STDIN>;
    wantarray ? @vars : \@vars;
}

sub get_sorted_features_for_org {
    my ($org) = @_;
    my $sap = new SAPserver;

    my $features = $sap->all_features( -ids => [$gid] )->{$gid}; 
    my $locH     = $sap->fid_locations( -ids => $features, -boundaries => 1 );
    my $funcH    = $sap->ids_to_functions( -ids => $features );
    # my $dnaH     = $sap->ids_to_sequences( -ids => $features );
    my $dnaH     = $sap->locs_to_dna( -locations => $locH );
    my $aliasH   = $sap->fids_to_ids( -ids => $features ); # LocusTag, NCBI, RefSeq, GeneID, GENE, Miscellaneous lists
    my @sorted   = sort { $a->[1] cmp $b->[1] || $a->[2] cmp $b->[2] ||
                          $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] }
                   map  { my $loc = $locH->{$_};
                          my ($contig, $beg, $end, $strand) = SeedUtils::parse_location($loc);
                          my ($org, $ctg) = split(/:/, $contig);
                          my $lo = $beg <= $end ? $beg : $end;
                          my $hi = $beg <= $end ? $end : $beg;
                          my $len = $hi - $lo + 1;
                          [ $_, $org, $ctg, $lo, $hi, $strand, $len, $loc, $funcH->{$_}, $aliasH->{$_}, uc $dnaH->{$_} ] 
                        } @$features; 
    my %features;
    for (@sorted) {
        my $ctg = $_->[2];
        push @{$features{$ctg}}, $_;
    }

    wantarray ? %features : \%features;
}

sub sum {
    my @x = @_;
    my $s;
    $s += $_ for @x;
    return $s;
}
