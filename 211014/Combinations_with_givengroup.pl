#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Algorithm::Combinatorics qw(combinations);

open( my $POOLIN,  "<", $ARGV[0] );
open( my $GROUPIN, "<", $ARGV[1] );

my $NUMBER = $ARGV[2];

my @POOL;

while (<$POOLIN>) {
    chomp;
    push( @POOL, $_ );
}

my %all;
while (<$GROUPIN>) {
    chomp;
    my @list = split "\t";
    my %POOL;
    for my $i (@list) {
        $POOL{$i} = 1;
    }
    my @temp_pool = grep { !exists( $POOL{$_} ); } @POOL;
    my $iter      = combinations( \@temp_pool, $NUMBER );
    while ( my $c = $iter->next ) {
        unless ( exists( $all{ join( "", sort( @list, @{$c} ) ) } ) ) {
            print( join( "\t", sort( @list, @{$c} ) ) );
            print("\n");
            $all{ join( "", sort( @list, @{$c} ) ) } = 1;
        }
    }
}

__END__
