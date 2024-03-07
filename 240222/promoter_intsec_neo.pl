#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;

my @chr_set;
foreach ( 1 .. 5 ) {
    $chr_set[$_] = AlignDB::IntSpan->new;
}

open GFF, "<", $ARGV[0];
while (<GFF>) {
    chomp;
    my @tmp   = split( "\t", $_ );
    my $chr   = $tmp[0];
    my $start = $tmp[1];
    my $end   = $tmp[2];
    $chr_set[$chr]->AlignDB::IntSpan::add_pair( $start, $end );
}
close GFF;

open SEG, "<", $ARGV[1];
while (<SEG>) {
    chomp;
    my @tmp  = split( " ", $_ );
    my $set1 = AlignDB::IntSpan->new;
    my $set2 = AlignDB::IntSpan->new;
    $set1->add_pair( $tmp[1], $tmp[2] );
    $set2->add_pair( $tmp[4], $tmp[5] );
    my $intsec1 = $set1->intersect( $chr_set[ $tmp[0] ] );
    my $intsec2 = $set2->intersect( $chr_set[ $tmp[3] ] );

    if ( ( $intsec1->cardinality >= 5 ) and ( $intsec2->cardinality >= 5 ) ) {
        my $rl1 = $intsec1->as_string;
        my $rl2 = $intsec2->as_string;
        print "$tmp[0] $tmp[1] $tmp[2] $rl1 $tmp[3] $tmp[4] $tmp[5] $rl2 \n";
    }
}
close SEG;
