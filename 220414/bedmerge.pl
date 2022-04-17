#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my ( $chr, $start, $end ) = ( 0, -1, -1 );

while (<STDIN>) {
    chomp;
    my @tmp = split;
    if ( $tmp[1] != $end ) {
        print "$chr\t$start\t$end\n";
        $chr   = $tmp[0];
        $start = $tmp[1];
        $end   = $tmp[2];
    }
    else {
        $end = $tmp[2];
    }
}

print "$chr\t$start\t$end\n";

__END__
