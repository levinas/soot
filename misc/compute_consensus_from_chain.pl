#! /usr/bin/env perl

use strict;

my $usage = "$0 ref.fa chain > consensus\n";

my $ref_file = shift @ARGV;
my $chain_file = shift @ARGV;

my @seqs = read_fasta($ref_file);
my $s = $seqs[0]->[2];

my @new_seqs;
my @lines = `cat $chain_file | grep chain`;
for (@lines) {
    my @cols = split /\s+/;
    my ($beg, $end_, $id) = @cols[5, 6, 7];
    my $len = $end_ - $beg;
    my $subseq = substr($s, $beg, $len);
    push @new_seqs, [$id, undef, $subseq];
}

write_fasta(\@new_seqs);

sub read_fasta
{
    my $dataR = ( $_[0] && ref $_[0] eq 'SCALAR' ) ?  $_[0] : slurp( @_ );
    $dataR && $$dataR or return wantarray ? () : [];

    my $is_fasta = $$dataR =~ m/^[\s\r]*>/;
    my @seqs = map { $_->[2] =~ tr/ \n\r\t//d; $_ }
               map { /^(\S+)([ \t]+([^\n\r]+)?)?[\n\r]+(.*)$/s ? [ $1, $3 || '', $4 || '' ] : () }
               split /[\n\r]+>[ \t]*/m, $$dataR;

    #  Fix the first sequence, if necessary
    if ( @seqs )
    {
        if ( $is_fasta )
        {
            $seqs[0]->[0] =~ s/^>//;  # remove > if present
        }
        elsif ( @seqs == 1 )
        {
            $seqs[0]->[1] =~ s/\s+//g;
            @{ $seqs[0] } = ( 'raw_seq', '', join( '', @{$seqs[0]} ) );
        }
        else  #  First sequence is not fasta, but others are!  Throw it away.
        {
            shift @seqs;
        }
    }

    wantarray() ? @seqs : \@seqs;
}

sub slurp
{
    my ( $fh, $close );
    if ( $_[0] && ref $_[0] eq 'GLOB' )
    {
        $fh = shift;
    }
    elsif ( $_[0] && ! ref $_[0] )
    {
        my $file = shift;
        if    ( -f $file                       ) { }
        elsif (    $file =~ /^<(.*)$/ && -f $1 ) { $file = $1 }  # Explicit read
        else                                     { return undef }
        open( $fh, '<', $file ) or return undef;
        $close = 1;
    }
    else
    {
        $fh = \*STDIN;
        $close = 0;
    }

    my $out = '';
    my $inc = 1048576;
    my $end =       0;
    my $read;
    while ( $read = read( $fh, $out, $inc, $end ) ) { $end += $read }
    close( $fh ) if $close;

    \$out;
}

sub write_fasta
{
    my ( $fh, $close, $unused ) = output_filehandle( shift );
    ( unshift @_, $unused ) if $unused;

    ( ref( $_[0] ) eq "ARRAY" ) or confess "Bad sequence entry passed to print_alignment_as_fasta\n";

    #  Expand the sequence entry list if necessary:

    if ( ref( $_[0]->[0] ) eq "ARRAY" ) { @_ = @{ $_[0] } }

    foreach my $seq_ptr ( @_ ) { print_seq_as_fasta( $fh, @$seq_ptr ) }

    close( $fh ) if $close;
}

sub output_filehandle
{
    my $file = shift;

    #  Null string or undef

    return ( \*STDOUT, 0 ) if ( ! defined( $file ) || ( $file eq "" ) );

    #  FILEHANDLE

    return ( $file, 0 ) if ( ref( $file ) eq "GLOB" );

    #  Some other kind of reference; return the unused value

    return ( \*STDOUT, 0, $file ) if ref( $file );

    #  File name

    my $fh;
    open( $fh, '>', $file ) || die "Could not open output $file\n";
    return ( $fh, 1 );
}

sub print_seq_as_fasta
{
    my $fh = ( ref $_[0] eq 'GLOB' ) ? shift : \*STDOUT;
    return if ( @_ < 2 ) || ( @_ > 3 ) || ! ( defined $_[0] && defined $_[-1] );
    #  Print header line
    print $fh  ( @_ == 3 && defined $_[1] && $_[1] =~ /\S/ ) ? ">$_[0] $_[1]\n" : ">$_[0]\n";
    #  Print sequence, 60 chars per line
    print $fh  join( "\n", $_[-1] =~ m/.{1,60}/g ), "\n";
}
