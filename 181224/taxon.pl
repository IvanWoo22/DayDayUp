#!/usr/bin/perl -w
use strict;

open IN, "<", $ARGV[0];
open OUT1, ">", $ARGV[1];
open OUT2, ">", $ARGV[2];

while (my $r1 = <IN>) {
    chomp($r1);
    chomp(my $r2 = <IN>);
    if (($r1 =~ m/M\t/) and ($r2 =~ m/M\t/)) {
        print OUT1 "$r1\n$r2\n";
    }
    else {
        print OUT2 "$r1\n$r2\n";
    }
}

close IN;
close OUT1;
close OUT2;
