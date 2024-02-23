#!/usr/bin/perl
use strict;
use warnings;
use AlignDB::IntSpan;

open CP,  "<", $ARGV[0];
open CR,  "<", $ARGV[1];
open BR,  "<", $ARGV[2];
open TC,  "<", $ARGV[3];
open ACR, ">", $ARGV[4];
open ABR, ">", $ARGV[5];
open ATC, ">", $ARGV[6];
open BTC, ">", $ARGV[7];
my $partial = $ARGV[8];

my @cr = <CR>;
my @br = <BR>;
my @tc = <TC>;
close CR;
close BR;
close TC;

my @mcx_cr;
my @mcx_br;
my @mcx_tc;

foreach ( 1 .. 5 ) {
    chomp( $cr[$_] );
    chomp( $br[$_] );
    chomp( $tc[$_] );
    my @tmp_cr = split( ": ", $cr[$_] );
    my @tmp_br = split( ": ", $br[$_] );
    my @tmp_tc = split( ": ", $tc[$_] );

    $mcx_cr[$_] = AlignDB::IntSpan->new;
    $mcx_br[$_] = AlignDB::IntSpan->new;
    $mcx_tc[$_] = AlignDB::IntSpan->new;

    $mcx_cr[$_]->add_runlist( $tmp_cr[1] );
    $mcx_br[$_]->add_runlist( $tmp_br[1] );
    $mcx_tc[$_]->add_runlist( $tmp_tc[1] );
}

while (<CP>) {
    chomp;
    my @tmp  = split( " ", $_ );
    my $set1 = AlignDB::IntSpan->new;
    my $set2 = AlignDB::IntSpan->new;
    my $out1;
    my $out2;
    if ( $tmp[1] < $tmp[2] ) {
        $out1 = $tmp[0] . " " . $tmp[1] . " " . $tmp[2];
        $set1->add_pair( $tmp[1], $tmp[2] );
    }
    else {
        $out1 = $tmp[0] . " " . $tmp[2] . " " . $tmp[1];
        $set1->add_pair( $tmp[2], $tmp[1] );
    }

    if ( $tmp[4] < $tmp[5] ) {
        $out2 = $tmp[3] . " " . $tmp[4] . " " . $tmp[5];
        $set2->add_pair( $tmp[4], $tmp[5] );
    }
    else {
        $out2 = $tmp[3] . " " . $tmp[5] . " " . $tmp[4];
        $set2->add_pair( $tmp[5], $tmp[4] );
    }

    my $cr_comp;
    if (
        (
            $set1->overlap( $mcx_cr[ $cov{ $tmp[0] } ] ) >=
            $partial * $set1->size
        )
        and (
            $set2->overlap( $mcx_cr[ $cov{ $tmp[3] } ] ) >=
            $partial * $set2->size )
      )
    {
        $cr_comp = 3;
    }
    elsif (
        $set2->overlap( $mcx_cr[ $cov{ $tmp[3] } ] ) >=
        $partial * $set2->size )
    {
        $cr_comp = 2;
    }
    elsif (
        $set1->overlap( $mcx_cr[ $cov{ $tmp[0] } ] ) >=
        $partial * $set1->size )
    {
        $cr_comp = 1;
    }
    else {
        $cr_comp = 0;
    }
    my $br_comp;
    if (
        (
            $set1->overlap( $mcx_br[ $cov{ $tmp[0] } ] ) >=
            $partial * $set1->size
        )
        and (
            $set2->overlap( $mcx_br[ $cov{ $tmp[3] } ] ) >=
            $partial * $set2->size )
      )
    {
        $br_comp = 3;
    }
    elsif (
        $set2->overlap( $mcx_br[ $cov{ $tmp[3] } ] ) >=
        $partial * $set2->size )
    {
        $br_comp = 2;
    }
    elsif (
        $set1->overlap( $mcx_br[ $cov{ $tmp[0] } ] ) >=
        $partial * $set1->size )
    {
        $br_comp = 1;
    }
    else {
        $br_comp = 0;
    }
    my $tc_comp;
    if (
        (
            $set1->overlap( $mcx_tc[ $cov{ $tmp[0] } ] ) >=
            $partial * $set1->size
        )
        and (
            $set2->overlap( $mcx_tc[ $cov{ $tmp[3] } ] ) >=
            $partial * $set2->size )
      )
    {
        $tc_comp = 3;
    }
    elsif (
        $set2->overlap( $mcx_tc[ $cov{ $tmp[3] } ] ) >=
        $partial * $set2->size )
    {
        $tc_comp = 2;
    }
    elsif (
        $set1->overlap( $mcx_tc[ $cov{ $tmp[0] } ] ) >=
        $partial * $set1->size )
    {
        $tc_comp = 1;
    }
    else {
        $tc_comp = 0;
    }

    if ( ( $cr_comp == 1 ) and ( $br_comp == 1 ) and ( $tc_comp == 1 ) ) {
        print ACR "$out1 $out2\n";
    }
    elsif ( ( $cr_comp == 2 ) and ( $br_comp == 2 ) and ( $tc_comp == 2 ) ) {
        print ACR "$out2 $out1\n";
    }
    elsif ( ( $cr_comp == 3 ) and ( ( $br_comp == 1 ) and ( $tc_comp == 1 ) ) )
    {
        print ABR "$out1 $out2\n";
    }
    elsif ( ( $cr_comp == 3 ) and ( ( $br_comp == 2 ) and ( $tc_comp == 2 ) ) )
    {
        print ABR "$out2 $out1\n";
    }
    elsif ( ( $cr_comp == 3 ) and ( $br_comp == 3 ) and ( $tc_comp == 1 ) ) {
        print ATC "$out1 $out2\n";
    }
    elsif ( ( $cr_comp == 3 ) and ( $br_comp == 3 ) and ( $tc_comp == 2 ) ) {
        print ATC "$out2 $out1\n";
    }
    elsif ( ( $cr_comp == 3 ) and ( $br_comp == 3 ) and ( $tc_comp == 3 ) ) {
        print BTC "$out1 $out2\n";
    }
}
close CP;
close ACR;
close ABR;
close ATC;
close BTC;