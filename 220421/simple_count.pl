#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $REGION, "<", $ARGV[0] );
my %regions;
my @order;
while (<$REGION>) {
    chomp;
    my @tmp = split;
    if ( exists( $regions{ $tmp[0] } ) ) {
        warn("There has been a [$tmp[0]] range!\n");
    }
    else {
        push( @order, $tmp[0] );
        foreach ( 0 .. $tmp[1] ) {
            ${ $regions{ $tmp[0] } }[$_] = 0;
        }
    }
}
close($REGION);

open( my $SITE, "<", $ARGV[1] );
while (<$SITE>) {
    chomp;
    my @tmp = split;
    ${ $regions{ $tmp[0] } }[ $tmp[1] ] = 1;
}
close($SITE);

foreach my $i ( 0 .. $#order ) {
    foreach my $j ( 0 .. $#{ $regions{ $order[$i] } } ) {
        print "${$regions{$order[$i]}}[$j]\n";
    }
}

__END__
