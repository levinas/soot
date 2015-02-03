#! /usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;

my $usage = "Usage: $0 gff 1.vcf [2.vcf ...] > var.table\n\n";

my ($help, $min_alt_depth, $min_alt_fract, $min_map_quality, $show_html);

GetOptions("h|help"     => \$help,
           "d=i"        => \$min_alt_depth,
           "f=f"        => \$min_alt_fract,
           "q=f"        => \$min_map_quality,
           "html"       => \$show_html,
	  ) or die("Error in command line arguments\n");

$help and die $usage;

my $gff_file = shift @ARGV;
my @vcf_files = grep { -s $_ } @ARGV;

$gff_file && @vcf_files or die $usage;

$min_alt_depth   ||= 5;
$min_alt_fract   ||= 0;
$min_map_quality ||= 0;

my $gff = read_gff_tree($gff_file);

print STDERR '$gff = '. Dumper($gff);

sub get_sorted_features_for_org {
    my ($org) = @_;

    my %features;

    # my $features = $sap->all_features( -ids => [$gid] )->{$gid};
    # my $locH     = $sap->fid_locations( -ids => $features, -boundaries => 1 );
    # my $funcH    = $sap->ids_to_functions( -ids => $features );
    # # my $dnaH     = $sap->ids_to_sequences( -ids => $features );
    # my $dnaH     = $sap->locs_to_dna( -locations => $locH );
    # my $aliasH   = $sap->fids_to_ids( -ids => $features ); # LocusTag, NCBI, RefSeq, GeneID, GENE, Miscellaneous lists
    # my @sorted   = sort { $a->[1] cmp $b->[1] || $a->[2] cmp $b->[2] ||
    #                       $a->[3] <=> $b->[3] || $a->[4] <=> $b->[4] }
    #                map  { my $loc = $locH->{$_};
    #                       my ($contig, $beg, $end, $strand) = SeedUtils::parse_location($loc);
    #                       my ($org, $ctg) = split(/:/, $contig);
    #                       my $lo = $beg <= $end ? $beg : $end;
    #                       my $hi = $beg <= $end ? $end : $beg;
    #                       my $len = $hi - $lo + 1;
    #                       [ $_, $org, $ctg, $lo, $hi, $strand, $len, $loc, $funcH->{$_}, $aliasH->{$_}, uc $dnaH->{$_} ]
    #                     } @$features;
    # for (@sorted) {
    #     my $ctg = $_->[2];
    #     push @{$features{$ctg}}, $_;
    # }

    wantarray ? %features : \%features;
}

sub read_gff_tree {
    my ($file) = @_;
    my $header = `cat $file | grep "^#"`;
    my @lines = `cat $file | grep -v "^#" |head`;
    my %id_to_index;
    my %rootH;
    my @features;
    my $index;
    for (@lines) {
        chomp;
        my ($seqname, $source, $feature, $start, $end, $score, $strand, $fname, $attribute) = split /\t/;
        my %hash = map { my ($k,$v) = split /=/; $k => $v } split(/;\s*/, $attribute);
        my $ent = { seqname => $seqname,
                    source => $source,
                    feature => $feature,
                    start => $start,
                    end => $end,
                    score => $end,
                    strand => $strand,
                    fname => $fname,
                    attribute => \%hash };
        my $parent = $attribute->{Parent};
        if (!$parent) {
            push @features, $ent;
            $id_to_index{$seqname} = $index++;
            next;
        }
        while ($parent) {
            $rootH{$seqname} = $parent;
            $parent = $rootH{$parent};
        }
        my $root_index = $id_to_index{$rootH{$seqname}};
        push @{$features[$root_index]->{descendants}}, $ent;
    }
    return \@features;
}
