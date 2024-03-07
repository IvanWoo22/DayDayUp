#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;

open IN, "<", $ARGV[0];
open LK, "<", $ARGV[1];
open OT, ">", $ARGV[2];

my @sit_chr;
my %sit_meth;
foreach ( 1 .. 5 ) {
    $sit_chr[$_] = AlignDB::IntSpan->new;
}

while (<IN>) {
    chomp;
    my @tmp = split( /\t/, $_ );
    $sit_chr[ $tmp[0] ]->AlignDB::IntSpan::add( $tmp[1] );
    my $tmp = $tmp[0] . "S" . $tmp[1];
    $sit_meth{$tmp} = $tmp[3];
}
close IN;

my %cov = ( "1" => 1, "2" => 2, "3" => 3, "4" => 4, "5" => 5 );

while (<LK>) {
    chomp;
    my @lin  = split( " ", $_ );
    my $set1 = AlignDB::IntSpan->new;
    my $set2 = AlignDB::IntSpan->new;
    $set1->add_pair( $lin[1], $lin[2] );
    $set2->add_pair( $lin[4], $lin[5] );
    my $sum1 = $sit_chr[ $cov{ $lin[0] } ]->AlignDB::IntSpan::overlap($set1);
    my $sum2 = $sit_chr[ $cov{ $lin[3] } ]->AlignDB::IntSpan::overlap($set2);
    my $arg_meth1 = "NULL";
    my $arg_meth2 = "NULL";

    if ( $sum1 > 0 ) {
        my @intsec1 =
          $sit_chr[ $cov{ $lin[0] } ]->AlignDB::IntSpan::intersect($set1)
          ->as_array;
        my $meth1 = 0;
        foreach ( 0 .. $#intsec1 ) {
            my $tmp = $cov{ $lin[0] } . "S" . $intsec1[$_];
            $meth1 += $sit_meth{$tmp};
        }
        $arg_meth1 = $meth1 / ( $#intsec1 + 1 );
    }
    if ( $sum2 > 0 ) {
        my @intsec2 =
          $sit_chr[ $cov{ $lin[3] } ]->AlignDB::IntSpan::intersect($set2)
          ->as_array;
        my $meth2 = 0;
        foreach ( 0 .. $#intsec2 ) {
            my $tmp = $cov{ $lin[3] } . "S" . $intsec2[$_];
            $meth2 += $sit_meth{$tmp};
        }
        $arg_meth2 = $meth2 / ( $#intsec2 + 1 );
    }
    print OT "$_\t$sum1\t$sum2\t$arg_meth1\t$arg_meth2\n";
}
close LK;
close OT;
