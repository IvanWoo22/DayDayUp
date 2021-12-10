#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

while (<STDIN>) {
    chomp;
    my @tmp = split "\t";
    if ( -e "temp/$tmp[1].csv" ) {
        my %count;
        foreach ( -20 .. 20 ) {
            $count{$_} = 0;
        }
        open( my $IN, "<", "temp/$tmp[1].csv" );
        readline($IN);
        while (<$IN>) {
            chomp;
            my @tmmp = split ",";
            if ( abs( $tmmp[1] - $tmp[0] ) <= 20 ) {
                if ( $tmp[2] eq "+" ) {
                    $count{ $tmmp[1] - $tmp[0] }++;
                }
                else {
                    $count{ $tmp[0] - $tmmp[1] }++;
                }
            }
            if ( $tmmp[1] == $tmp[0] ) {
                warn("$tmp[0]\t$tmp[1]\n");
            }
        }
        if ( $count{-20} > 0 ) {
            print("1");
        }
        else {
            print("0");
        }
        foreach ( -19 .. 20 ) {
            if ( $count{$_} > 0 ) {
                print("\t1");
            }
            else {
                print("\t0");
            }
        }
        print("\n");
    }
}
