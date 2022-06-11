#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open(my $TAXON, "<", $ARGV[0]);
my %info_tbl;
while (<$TAXON>) {
    chomp;
    my @tbl = split /\t/;
    $info_tbl{$tbl[0]} = $tbl[4] . "\t" . $tbl[5] . "\t" . $tbl[6];
}
close($TAXON);

while (<STDIN>) {
    chomp;
    my $name = $_;
    foreach my $KEY (keys(%info_tbl)) {
        if ($name =~ /$KEY/) {
            print("$name\t$info_tbl{$KEY}\n");
        }
    }
}