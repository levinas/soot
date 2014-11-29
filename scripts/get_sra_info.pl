#! /usr/bin/env perl

use strict;

use Carp;
use Data::Dumper;
use Getopt::Long;
use LWP::Simple;
use XML::Simple;
use JSON;

my $usage = <<"End_of_Usage";

usage:  $0 [ options ] < table-with-SRP-IDs > extended-table-with-SRR-info 2> ids.not.found

Options:
        -c int            - 1-based number of column that contains SRP or SRA ID
        -j                - output reads information in JSON
        -p                - output reads information in Assembly-RAST CLI parameters

Examples: 

        $0 SRA068895 SRP020228 SRP020230 > table
        $0 SRA068895 -j > json

End_of_Usage

my ($help, $col, $json_out, $param_out, $download, $match_id);

GetOptions("h|help"  => \$help,
           "c=i"     => \$col,
           "d"       => \$download,
           "j|json"  => \$json_out,
           "m|match" => \$match_id,
           "p|param" => \$param_out,
	  ) or die("Error in command line arguments\n");

$help and die $usage;

my @ids = @ARGV;
my $db = 'SRA';

my @input_rows = @ids ? map { [$_] } @ids : map { chomp; [ split/\t/ ] } <STDIN>;

for (@input_rows) {
    my $line = join("\t", @$_);
    my $id = $col ? $_->[$col-1] : $_->[-1];
    my $data = ncbi_search($db, $id);
    my $opts = { match_id => $match_id ? $id : undef, download => $download };
    if ($json_out) {
        print encode_json($data);
    } elsif ($param_out) {
        my $params = arast_param_for_run_set($data, $opts) or next;
        print $params."\n";
    } else {
        my @rows = key_run_set_stats($data, $opts);
        if (!@rows) {
            print STDERR "No information found for $id\n";
            next;
        } 
        for my $row (@rows) {
            print join("\t", $line, @$row)."\n";
        }
    }
}

sub download_run {
    my ($srr) = @_;
    unless (-s "$srr.fastq" || -s "$srr\_1.fastq" && -s "$srr\_2.fastq") {
        run("fastq-dump --split-3 $srr");
    }
}

sub arast_param_for_run_set {
    my ($hash, $opts) = @_;
    my $match_id = $opts->{match_id};
    my $download = $opts->{download};
    my $pkg = $hash->{EXPERIMENT_PACKAGE};
    my @set = ref($pkg) eq 'ARRAY' ? @$pkg : $pkg ? ($pkg) : ();
    my @params;
    my $desc;
    for (@set) {
        my $title = $_->{SAMPLE}->{TITLE};
        my $submission = $_->{SUBMISSION}->{accession};
        my $study = $_->{EXPERIMENT}->{STUDY_REF}->{IDENTIFIERS}->{PRIMARY_ID};
        my $study2 = $_->{STUDY}->{IDENTIFIERS}->{PRIMARY_ID};
        my $rs = $_->{RUN_SET}->{RUN};
        my $acc = $rs->{accession};
        my $bases = $rs->{total_bases};
        my $nreads = $rs->{Statistics}->{nreads};
        my $reads = $rs->{Statistics}->{Read};
        my @lens = map { $_->{average} } @$reads;
        my $rlen = mean(\@lens);

        # Usually: SRP = $study, SRA = $submission
        next if $match_id && $study ne $match_id && $submission ne $match_id;

        download_run($acc) if $download;

        my $param = ($nreads == 1) ? "--single $acc.fastq" :
                    ($nreads == 2) ? "--pair $acc\_1.fastq $acc\_2.fastq" : undef or next;

        $desc ||= "-m '$study: $title'";
        push @params, $param;
    }
    return unless @params;
    join(" ", $desc, @params);
}
    
sub key_run_set_stats {
    my ($hash, $opts) = @_;
    my $match_id = $opts->{match_id};
    my $pkg = $hash->{EXPERIMENT_PACKAGE};
    my @set = ref($pkg) eq 'ARRAY' ? @$pkg : $pkg ? ($pkg) : ();
    my @rows;
    for (@set) {
        my $title = $_->{SAMPLE}->{TITLE};
        my $study = $_->{EXPERIMENT}->{STUDY_REF}->{IDENTIFIERS}->{PRIMARY_ID};
        my $study2 = $_->{STUDY}->{IDENTIFIERS}->{PRIMARY_ID};
        my $submission = $_->{SUBMISSION}->{accession};
        my $rs = $_->{RUN_SET}->{RUN};
        my $acc = $rs->{accession};
        my $bases = $rs->{total_bases};
        my $nreads = $rs->{Statistics}->{nreads};
        my $reads = $rs->{Statistics}->{Read};
        my @lens = map { $_->{average} } @$reads;
        my $rlen = mean(\@lens);

        # Usually: SRP = $study, SRA = $submission
        # print join("\t", $match_id, $study, $submission) . "\n";
        
        next if $match_id && $study ne $match_id && $submission ne $match_id;
        push @rows, [$acc, $nreads, $rlen, $bases, $submission, $study, $title];
    }
    wantarray ? @rows : \@rows;
}

sub mean {
    my ($array) = @_;
    my $sum;
    $sum += $_ for @$array;
    @$array ? $sum / @$array : undef;
}

sub ncbi_search {
    my ($db, $query) = @_;

    my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    my $url = $base . "esearch.fcgi?db=$db&term=$query&usehistory=y";
    my $output = get($url);
    # print $output;

    my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
    my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
    $url = $base . "esummary.fcgi?db=$db&query_key=$key&WebEnv=$web";
    my $docsums = get($url);
    # print "$docsums";

    $url = $base . "efetch.fcgi?db=$db&query_key=$key&WebEnv=$web";
    $url .= "&rettype=abstract&retmode=text";
    my $data = get($url);
    return XMLin($data);
}

sub run { system(@_) == 0 or confess("FAILED: ". join(" ", @_)); }
