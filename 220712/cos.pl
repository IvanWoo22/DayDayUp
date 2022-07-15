#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub cosine {
    my @TMP1 = split( "\t", $_[0] );
    my @TMP2 = split( "\t", $_[1] );
    if ( $#TMP1 == $#TMP2 ) {
        my $MULTI_SUM = 0;
        my $TMP1_DIS  = 0;
        my $TMP2_DIS  = 0;
        foreach my $CID ( 0 .. $#TMP1 ) {
            $MULTI_SUM += ( $TMP1[$CID] * $TMP2[$CID] );
            $TMP1_DIS  += $TMP1[$CID]**2;
            $TMP2_DIS  += $TMP2[$CID]**2;
        }
        my $COS_DISTANCE;
        if ( sqrt($TMP1_DIS) * sqrt($TMP2_DIS) > 0 ) {
            $COS_DISTANCE =
              1 - $MULTI_SUM / ( sqrt($TMP1_DIS) * sqrt($TMP2_DIS) );
        }
        else {
            $COS_DISTANCE = 1;
        }
        return ($COS_DISTANCE);
    }
    else {
        warn("The two vectors were not in equal length!\n");
    }
}

my ( @chr, @position, @pv, @vector );
while (<STDIN>) {
    chomp;
    my @tmp = split;
    push( @chr,      $tmp[0] );
    push( @position, $tmp[1] );
    push( @pv,       $tmp[2] );
    push( @vector,   join( "\t", @tmp[ 4 .. 45 ] ) );
}

foreach my $id1 ( 0 .. $#chr ) {
    print("\t$chr[$id1]_$position[$id1]");
}
print("\n");

foreach my $id1 ( 0 .. $#chr ) {
    print("$chr[$id1]_$position[$id1]");
    foreach my $id2 ( 0 .. $#chr ) {
        my $distance = cosine( $vector[$id1], $vector[$id2] );
        if ( ( $distance == 1 ) or ( $distance <= 0 ) ) {
            printf( "\t%.0f", $distance );
        }
        else {
            printf( "\t%.3f", $distance );
        }
    }
    print("\n");
}

__END__
