#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my ( @chr, @range, @depth );
while (<STDIN>) {
    chomp;
    my @tmp = split;
    push( @chr,   $tmp[0] );
    push( @range, $tmp[1] );
    push( @depth, $tmp[2] );
}

my ( $chr, $start, $end, $count, $gap ) = ( 0, -1, -1, 0, $ARGV[0] );
my @dep_count;
foreach my $id ( 0 .. $#chr ) {
    if ( ( $range[$id] - $end > $gap ) || ( $chr[$id] ne $chr ) ) {
        if ( $count > 0 ) {
            my $length = $end - $start;
            print "$chr\t$start\t$end\t$length\t$count\t";
            print( join( "\t", @dep_count[ 1 .. 4 ] ) );
            print "\n";
        }
        $chr                 = $chr[$id];
        $start               = $range[$id];
        $count               = 1;
        @dep_count[ 1 .. 4 ] = ( 0, 0, 0, 0 );
    }
    else {
        $count++;
    }
    $dep_count[ $depth[$id] ]++;
    $end = $range[$id] + 1;
}

my $length = $end - $start;
print("$chr\t$start\t$end\t$length\t$count\t");
print( join( "\t", @dep_count[ 1 .. 4 ] ) );
print("\n");

__END__
