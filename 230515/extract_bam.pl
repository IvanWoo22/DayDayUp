#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

while (<>) {
    next if s/^#//;
    chomp;
    my $record = 0;
    my @info   = split;
    if ( $info[5] ne "*" ) {
        while ( $info[5] =~ m/([0-9]+)I/g ) {
            $record = $record + $1;
        }
        while ( $info[5] =~ m/([0-9]+)S/g ) {
            $record = $record + $1;
        }
        if ( $record < 30 ) {
            print "$info[0]\n";
        }
    }
}

__END__