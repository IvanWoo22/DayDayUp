#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Term::ProgressBar;

my %list;
my @table;

open my $WHITELIST, "<", $ARGV[0];
my $sum = 0;
while (<$WHITELIST>) {
    chomp;
    my ( $sample, undef ) = split /\t/;
    foreach my $j ( 0 .. 21 ) {
        $table[$sum][$j] = 0;
    }
    $list{$sample} = $sum;
    $sum++;
}
close($WHITELIST);

open my $SAM, "<", $ARGV[1];
my $line     = `wc -l < $ARGV[1]`;
my $progress = Term::ProgressBar->new($line);
$line = 0;
while (<$SAM>) {
    my @tmp = split /\t/;
    my ( undef, $name, undef ) = split( /_/, $tmp[0] );
    $line++;
    if ( ( defined($name) ) and ( exists( $list{$name} ) ) ) {
        my $seq = $tmp[9];

        if ( $seq =~
/(?:(?:A(?:G(?:T(?:T(?:[C,T](?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.CTT[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GCTT[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][0]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:T(?:T(?:[C,T](?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TTT[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GTTT[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][1]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.C[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TC[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TTC[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TTC[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TTC[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][2]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TTT[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TTT[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TTT[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][3]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.C[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]C[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T]C[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T]C[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T]C[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T]C[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][4]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]T[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T]T[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T]T[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T]T[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T]T[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][5]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.CGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]CGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T][C,T]CGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T][C,T]CGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T][C,T]CGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T][C,T]CGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T][C,T]CGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][6]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T][C,T]TGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T][C,T]TGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T][C,T]TGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T][C,T]TGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T][C,T]TGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][7]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:[C,T](?:G(?:G(?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.CAA[C,T]AATTAATAGA[C,T]TGGAT)|.GCAA[C,T]AATTAATAGA[C,T]TGGAT)|.GGCAA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]GGCAA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T][C,T]GGCAA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T][C,T][C,T]GGCAA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T][C,T][C,T]GGCAA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T][C,T][C,T]GGCAA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T][C,T][C,T]GGCAA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T][C,T][C,T]GGCAA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][8]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:[C,T](?:G(?:G(?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.TAA[C,T]AATTAATAGA[C,T]TGGAT)|.GTAA[C,T]AATTAATAGA[C,T]TGGAT)|.GGTAA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]GGTAA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T][C,T]GGTAA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T][C,T][C,T]GGTAA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T][C,T][C,T]GGTAA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T][C,T][C,T]GGTAA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T][C,T][C,T]GGTAA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T][C,T][C,T]GGTAA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][9]++;
        }

        if ( $seq =~
/(?^:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.CAATTAATAGA[C,T]TGGAT)|.ACAATTAATAGA[C,T]TGGAT)|.AACAATTAATAGA[C,T]TGGAT)|.[C,T]AACAATTAATAGA[C,T]TGGAT)|.G[C,T]AACAATTAATAGA[C,T]TGGAT)|.GG[C,T]AACAATTAATAGA[C,T]TGGAT)|.[C,T]GG[C,T]AACAATTAATAGA[C,T]TGGAT)|.[C,T][C,T]GG[C,T]AACAATTAATAGA[C,T]TGGAT)|.[C,T][C,T][C,T]GG[C,T]AACAATTAATAGA[C,T]TGGAT)|.T[C,T][C,T][C,T]GG[C,T]AACAATTAATAGA[C,T]TGGAT)|.TT[C,T][C,T][C,T]GG[C,T]AACAATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T][C,T][C,T]GG[C,T]AACAATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T][C,T][C,T]GG[C,T]AACAATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][10]++;
        }

        if ( $seq =~
/(?^:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.TAATTAATAGA[C,T]TGGAT)|.ATAATTAATAGA[C,T]TGGAT)|.AATAATTAATAGA[C,T]TGGAT)|.[C,T]AATAATTAATAGA[C,T]TGGAT)|.G[C,T]AATAATTAATAGA[C,T]TGGAT)|.GG[C,T]AATAATTAATAGA[C,T]TGGAT)|.[C,T]GG[C,T]AATAATTAATAGA[C,T]TGGAT)|.[C,T][C,T]GG[C,T]AATAATTAATAGA[C,T]TGGAT)|.[C,T][C,T][C,T]GG[C,T]AATAATTAATAGA[C,T]TGGAT)|.T[C,T][C,T][C,T]GG[C,T]AATAATTAATAGA[C,T]TGGAT)|.TT[C,T][C,T][C,T]GG[C,T]AATAATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T][C,T][C,T]GG[C,T]AATAATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T][C,T][C,T]GG[C,T]AATAATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][11]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.CTGGAT)|.ACTGGAT)|.GACTGGAT)|.AGACTGGAT)|.TAGACTGGAT)|.ATAGACTGGAT)|.AATAGACTGGAT)|.TAATAGACTGGAT)|.TTAATAGACTGGAT)|.ATTAATAGACTGGAT)|.AATTAATAGACTGGAT)|.[C,T]AATTAATAGACTGGAT)|.A[C,T]AATTAATAGACTGGAT)|.AA[C,T]AATTAATAGACTGGAT)|.[C,T]AA[C,T]AATTAATAGACTGGAT)|.G[C,T]AA[C,T]AATTAATAGACTGGAT)|.GG[C,T]AA[C,T]AATTAATAGACTGGAT)|.[C,T]GG[C,T]AA[C,T]AATTAATAGACTGGAT)|.[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGACTGGAT)|.[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGACTGGAT)|.T[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGACTGGAT)|.TT[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGACTGGAT)|.[C,T]TT[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGACTGGAT)|.G[C,T]TT[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGACTGGAT))/
          )
        {
            $table[ $list{$name} ][12]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:[C,T](?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TTGGAT)|.ATTGGAT)|.GATTGGAT)|.AGATTGGAT)|.TAGATTGGAT)|.ATAGATTGGAT)|.AATAGATTGGAT)|.TAATAGATTGGAT)|.TTAATAGATTGGAT)|.ATTAATAGATTGGAT)|.AATTAATAGATTGGAT)|.[C,T]AATTAATAGATTGGAT)|.A[C,T]AATTAATAGATTGGAT)|.AA[C,T]AATTAATAGATTGGAT)|.[C,T]AA[C,T]AATTAATAGATTGGAT)|.G[C,T]AA[C,T]AATTAATAGATTGGAT)|.GG[C,T]AA[C,T]AATTAATAGATTGGAT)|.[C,T]GG[C,T]AA[C,T]AATTAATAGATTGGAT)|.[C,T][C,T]GG[C,T]AA[C,T]AATTAATAGATTGGAT)|.[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGATTGGAT)|.T[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGATTGGAT)|.TT[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGATTGGAT)|.[C,T]TT[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGATTGGAT)|.G[C,T]TT[C,T][C,T][C,T]GG[C,T]AA[C,T]AATTAATAGATTGGAT))/
          )
        {
            $table[ $list{$name} ][13]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.CC[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TCC[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TTCC[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TTCC[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TTCC[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][14]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.CT[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TCT[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TTCT[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TTCT[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TTCT[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][15]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TC[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TTC[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TTTC[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TTTC[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TTTC[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][16]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TTT[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TTTT[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TTTT[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TTTT[C,T]GG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][17]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.CCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]CCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T]CCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T]CCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T]CCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T]CCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][18]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.CTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]CTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T]CTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T]CTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T]CTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T]CTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][19]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T]TCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T]TCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T]TCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T]TCGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][20]++;
        }

        if ( $seq =~
/(?:(?:A(?:G(?:[C,T](?:T(?:T(?:[C,T](?:G(?:G(?:[C,T](?:A(?:A(?:[C,T](?:A(?:A(?:T(?:T(?:A(?:A(?:T(?:A(?:G(?:A(?:[C,T](?:T(?:G(?:G(?:A.|.T)|.AT)|.GAT)|.GGAT)|.TGGAT)|.[C,T]TGGAT)|.A[C,T]TGGAT)|.GA[C,T]TGGAT)|.AGA[C,T]TGGAT)|.TAGA[C,T]TGGAT)|.ATAGA[C,T]TGGAT)|.AATAGA[C,T]TGGAT)|.TAATAGA[C,T]TGGAT)|.TTAATAGA[C,T]TGGAT)|.ATTAATAGA[C,T]TGGAT)|.AATTAATAGA[C,T]TGGAT)|.[C,T]AATTAATAGA[C,T]TGGAT)|.A[C,T]AATTAATAGA[C,T]TGGAT)|.AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.T[C,T]TTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.TT[C,T]TTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.[C,T]TT[C,T]TTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT)|.G[C,T]TT[C,T]TTGG[C,T]AA[C,T]AATTAATAGA[C,T]TGGAT))/
          )
        {
            $table[ $list{$name} ][21]++;
        }
    }
    $progress->update($line);
}
close($SAM);

foreach my $i ( 0 .. $sum - 1 ) {
    foreach my $j ( 0 .. 6 ) {
        print( join( "\t", @{ $table[$i] }[ $j * 2 .. $j * 2 + 1 ] ) );
        my $cov = 0;
        $cov =
          $table[$i][ $j * 2 + 1 ] /
          ( $table[$i][ $j * 2 + 1 ] + $table[$i][ $j * 2 ] )
          if ( $table[$i][ $j * 2 + 1 ] + $table[$i][ $j * 2 ] > 0 );
        print("\t$cov\t");
    }
    print( join( "\t", @{ $table[$i] }[ 14 .. 21 ] ) );
    print "\n";
}

__END__
