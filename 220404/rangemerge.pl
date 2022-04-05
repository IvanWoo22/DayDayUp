#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my @chr;
my @range;
my @depth;
while (<STDIN>) {
    chomp;
    my @tmp = split;
    push( @chr,   $tmp[0] );
    push( @range, $tmp[1] );
    push( @depth, $tmp[3] );
}

my ( $chr, $start, $end, $count, $dep_count ) = ( 0, -1, -1, 0, 0 );

foreach my $id ( 0 .. $#chr ) {
    if ( $range[$id] != $end ) {
        if ( $count > 0 ) {
            my $ave_dep = $dep_count / $count;
            my $length  = $end - $start;
            print "$chr\t$start\t$end\t$length\t$ave_dep\n";
        }
        $chr       = $chr[$id];
        $start     = $range[$id];
        $end       = $range[$id] + 1000;
        $count     = 1;
        $dep_count = $depth[$id];
    }
    else {
        $end = $range[$id] + 1000;
        $count++;
        $dep_count = $dep_count + $depth[$id];
    }
}

my $ave_dep = $dep_count / $count;
my $length  = $end - $start;
print "$chr\t$start\t$end\t$length\t$ave_dep\n";

__END__
