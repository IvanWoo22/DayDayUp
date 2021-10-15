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
        print( join( "\t", @list ) );
        print("\t");
        print( join( "\t", @{$c} ) );
        print("\n");
    }
}

__END__
