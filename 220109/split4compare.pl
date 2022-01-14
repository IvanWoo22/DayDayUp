#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub INNER_DIF {
    my $RAW_ARRAY = shift;
    my $LIMIT     = shift;
    my $COUNT     = 0;
    foreach ( 0 .. $#{$RAW_ARRAY} - 1 ) {
        if ( abs( ${$RAW_ARRAY}[$_] - ${$RAW_ARRAY}[ $_ + 1 ] ) > $LIMIT ) {
            $COUNT++;
        }
    }
    if ( $COUNT > 0 ) {
        return (0);
    }
    else {
        return (1);
    }
}

open( my $TSV, "<", $ARGV[0] );

readline($TSV);
my ( $count1, $count2 ) = ( 0, 0 );
while (<$TSV>) {
    chomp;
    my @tmp = split /\t/;
    my @val_ch1;
    my @val_ch2;
    foreach my $i ( 1 .. $ARGV[1] ) {
        if ( $tmp[$i] ne "NA" ) {
            push( @val_ch1, $tmp[$i] );
        }
    }
    foreach my $i ( $ARGV[1] + 1 .. $#tmp ) {
        if ( $tmp[$i] ne "NA" ) {
            push( @val_ch2, $tmp[$i] );
        }
    }
    if ( ( $#val_ch1 > $ARGV[2] - 2 ) && ( $#val_ch2 > $ARGV[3] - 2 ) ) {
        if (   ( INNER_DIF( \@val_ch1, $ARGV[4] ) )
            && ( INNER_DIF( \@val_ch2, $ARGV[4] ) ) )
        {
            if (   ( $tmp[ $ARGV[1] ] ne "NA" )
                && ( $tmp[ $ARGV[1] + 1 ] ne "NA" )
                && ( abs( $tmp[ $ARGV[1] ] - $tmp[ $ARGV[1] + 1 ] ) ) >
                $ARGV[5] )
            {
                print( join( "\t", @tmp ) );
                print "\n";
            }
            $count2++;
        }
        $count1++;
    }
}

close($TSV);

warn("$count1\t$count2\n");

__END__
