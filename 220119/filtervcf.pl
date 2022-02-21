#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

while (<STDIN>) {
    chomp;
    my @tmp      = split;
    my $base_var = 0;
    my ( undef, $a_count, undef ) = split( ":", $tmp[5] );
    my ( undef, $c_count, undef ) = split( ":", $tmp[6] );
    my ( undef, $g_count, undef ) = split( ":", $tmp[7] );
    my ( undef, $t_count, undef ) = split( ":", $tmp[8] );
    if ( $a_count >= 5 ) {
        $base_var++;
    }
    if ( $c_count >= 5 ) {
        $base_var++;
    }
    if ( $g_count >= 5 ) {
        $base_var++;
    }
    if ( $t_count >= 5 ) {
        $base_var++;
    }
    if ( ( $tmp[3] >= 10 ) and ( $base_var > 1 ) ) {
        print( join( "\t", @tmp ) );
        print("\t$base_var\n");
    }
}
__END__
