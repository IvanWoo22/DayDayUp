#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my ( $position, $near_site ) = ( 0, 0 );
while (<STDIN>) {
    chomp;
    my @tmp = split;
    unless ( ( $tmp[1] eq $near_site ) and ( $tmp[3] eq $position ) ) {
        print( join( "\t", @tmp ) );
    }
    $near_site = $tmp[3];
    $position  = $tmp[1];
}

__END__
