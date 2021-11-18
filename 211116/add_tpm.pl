#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my %list;

open my $LS, "<", $ARGV[0];
while (<$LS>) {
    chomp;
    my $line = $_;
    my @tmp  = split /\t/, $line;
    $list{ $tmp[0] } = $line;
}
close($LS);

open my $FA, "<", $ARGV[1];
my $no = $ARGV[2];
while (<$FA>) {
    chomp;
    my $line = $_;
    my @tmp  = split /\t/, $line;
    my @out  = split /\//, $tmp[ $no - 1 ];
    print("$line");
    foreach (@out) {
        print("\t$list{$_}");
    }
    print("\n");
}
close($FA);
