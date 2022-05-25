#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub euclidean {
    my @TMP1 = split( "\t", $_[0] );
    my @TMP2 = split( "\t", $_[1] );
    if ( $#TMP1 == $#TMP2 ) {
        my $SQUARE_SUM = 0;
        foreach my $CID ( 0 .. $#TMP1 ) {
            $SQUARE_SUM += ( $TMP1[$CID] - $TMP2[$CID] )**2;
        }
        my $EUL_DISTANCE = sqrt($SQUARE_SUM);
        return ($EUL_DISTANCE);
    }
    else {
        warn("The two vectors were not in equal length!\n");
    }
}

my ( @chr, @position, @pv, @vector );
push( @chr,      "-1" );
push( @position, "-1" );
push( @pv,       "-1" );
push( @vector,   "-1" );
while (<STDIN>) {
    chomp;
    my @tmp = split;
    push( @chr,      $tmp[0] );
    push( @position, $tmp[1] );
    push( @pv,       $tmp[2] );
    push( @vector,   join( "\t", @tmp[ 4 .. 45 ] ) );
}
push( @chr,      "-1" );
push( @position, "-1" );
push( @pv,       "-1" );
push( @vector,   "-1" );

foreach my $id ( 1 .. ( $#chr - 1 ) ) {
    my $length    = 1000;
    my $distance  = "NA";
    my $near_site = "NA";
    if (   ( abs( $position[$id] - $position[ $id - 1 ] ) <= 100 )
        && ( $chr[$id] eq $chr[ $id - 1 ] ) )
    {
        $near_site = $position[ $id - 1 ];
        $length    = abs( $position[$id] - $position[ $id - 1 ] );
        $distance  = euclidean( $vector[$id], $vector[ $id - 1 ] );
    }
    if (   ( abs( $position[$id] - $position[ $id + 1 ] ) <= 100 )
        && ( $chr[$id] eq $chr[ $id + 1 ] ) )
    {
        if ( $length > abs( $position[$id] - $position[ $id + 1 ] ) ) {
            $near_site = $position[ $id + 1 ];
            $length    = abs( $position[$id] - $position[ $id + 1 ] );
            $distance  = euclidean( $vector[$id], $vector[ $id + 1 ] );
        }
    }
    print(
"$chr[$id]\t$position[$id]\t$pv[$id]\t$near_site\t$distance\t$length\t$vector[$id]\n"
    );
}

__END__
