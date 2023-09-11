#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;

my @chr_set;
foreach ( 1 .. 5 ) {
    $chr_set[$_] = AlignDB::IntSpan->new;
}
my %cov = ( "I" => 1, "II" => 2, "III" => 3, "IV" => 4, "V" => 5 );

open SEG, "<", $ARGV[0];
while (<SEG>) {
    chomp;
    my @tmp = split( " ", $_ );
    $chr_set[ $cov{ $tmp[0] } ]->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
    $chr_set[ $cov{ $tmp[3] } ]->AlignDB::IntSpan::add_pair( $tmp[4], $tmp[5] );
}
close SEG;

open GFF, "<", $ARGV[1];
while (<GFF>) {
    chomp;
    my @tmp = split( "\t", $_ );
    my $set = AlignDB::IntSpan->new;
    $set->add_pair( $tmp[1], $tmp[2] );
    my $intsec = $set->intersect( $chr_set[ $tmp[0] ] );
    if ( $intsec->cardinality >= 0.8 * ( $tmp[2] - $tmp[1] + 1 ) ) {
        print "$_\n";
    }
}
close GFF;
