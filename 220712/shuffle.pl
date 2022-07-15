#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ( $i = @$array ; --$i ; ) {
        my $j = int rand( $i + 1 );
        next if $i == $j;
        @$array[ $i, $j ] = @$array[ $j, $i ];
    }
}

my ( @info, @vector );
while (<STDIN>) {
    chomp;
    my @tmp = split;
    push( @info,   join( "\t", @tmp[ 0 .. 3 ] ) );
    push( @vector, join( "\t", @tmp[ 4 .. 45 ] ) );
}

foreach my $id ( 0 .. $#info ) {
    print("$info[$id]\t");
    my @tmp = split "\t", $vector[$id];
    fisher_yates_shuffle( \@tmp );
    print( join( "\t", @tmp ) . "\n" );
}

__END__
