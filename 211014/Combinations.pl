#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Algorithm::Combinatorics qw(combinations);

open( my $POOLIN, "<", $ARGV[0] );
my $NUMBER = $ARGV[1];

my @POOL;

while (<$POOLIN>) {
    chomp;
    push( @POOL, $_ );
}

my $iter = combinations( \@POOL, $NUMBER );

while ( my $c = $iter->next ) {
    print( join( "\t", @{$c} ) );
    print("\n");
}

__END__
