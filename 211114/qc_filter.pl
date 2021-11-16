#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $TSV,  "<", $ARGV[0] );
open( my $TSV1, ">", $ARGV[1] );
open( my $TSV2, ">", $ARGV[2] );

sub IS_NA {
    my @TMP = @_;
    my $I   = 0;
    foreach my $VAR (@TMP) {
        if ( $VAR ne "NA" ) {
            $I++;
        }
    }
    return ($I);
}

while (<$TSV>) {
    chomp;
    my $line = $_;
    my @tmp  = split( /,/, $line );
    $tmp[0] =~ s/"//g;
    my $sums = 0;
    my $alls = 0;
    foreach my $j ( 1 .. 15 ) {
        my $k = ( $j - 1 ) * 7 + 1;
        if (   ( IS_NA( @tmp[ $k, $k + 1, $k + 6 ] ) == 3 )
            && ( IS_NA( @tmp[ $k + 2 .. $k + 5 ] ) > 1 ) )
        {
            $sums++;
        }
        if ( IS_NA( @tmp[ $k .. $k + 6 ] ) == 7 ) {
            $alls++;
        }
    }
    if ( $sums > 6 ) {
        print $TSV1 "$line\n";
    }
    if ( $alls > 3 ) {
        print $TSV2 "$line\n";
    }
}
close($TSV);

__END__