#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my %list;

open my $LS, "<", $ARGV[0];
while (<$LS>) {
    chomp;
    $list{$_} = 1;
}
close($LS);

open my $FA, "<", $ARGV[1];
while (<$FA>) {
    chomp;
    my $line = $_;
    my @tmp  = split ",", $line;
    my ( $gene, undef ) = split /\./, $tmp[0];
    $tmp[0] = $gene;
    if ( exists( $list{$gene} ) ) {
        print( join( "\t", @tmp ) );
        print("\tT\n");
    }
    else {
        print( join( "\t", @tmp ) );
        print("\t\n");
    }
}
close($FA);
