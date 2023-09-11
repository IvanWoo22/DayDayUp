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

my %cov = ( "I" => 1, "II" => 2, "III" => 3, "IV" => 4, "V" => 5 );

while (<LK>) {
    chomp;
    my @lin  = split( "\t", $_ );
    my $set1 = AlignDB::IntSpan->new;
    my $set2 = AlignDB::IntSpan->new;
    $set1->AlignDB::IntSpan::add( $lin[7] );
    $set2->AlignDB::IntSpan::add( $lin[9] );
    my $sum1   = $sit_chr[ $cov{ $lin[0] } ]->AlignDB::IntSpan::overlap($set1);
    my $sum2   = $sit_chr[ $cov{ $lin[3] } ]->AlignDB::IntSpan::overlap($set2);
    my $count1 = () = $lin[6] =~ /,/g;
    my $count2 = () = $lin[8] =~ /,/g;

    if ( ( $sum1 > 3 * $count1 ) and ( $sum2 > 3 * $count2 ) ) {
        my @intsec1 =
          $sit_chr[ $cov{ $lin[0] } ]->AlignDB::IntSpan::intersect($set1)
          ->as_array;
        my @intsec2 =
          $sit_chr[ $cov{ $lin[3] } ]->AlignDB::IntSpan::intersect($set2)
          ->as_array;

        my $meth1 = 0;
        foreach ( 0 .. $#intsec1 ) {
            my $tmp = $cov{ $lin[0] } . "S" . $intsec1[$_];
            $meth1 += $sit_meth{$tmp};
        }
        my $arg_meth1 = $meth1 / ( $#intsec1 + 1 );

        my $meth2 = 0;
        foreach ( 0 .. $#intsec2 ) {
            my $tmp = $cov{ $lin[3] } . "S" . $intsec2[$_];
            $meth2 += $sit_meth{$tmp};
        }
        my $arg_meth2 = $meth2 / ( $#intsec2 + 1 );

        my $delta = $arg_meth2 - $arg_meth1;
        print OT join( "\t", @lin )
          . "\t$delta\t$sum1\t$sum2\t$arg_meth1\t$arg_meth2\n";
    }
}
close LK;
close OT;
