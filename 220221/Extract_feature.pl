#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $LIST, "<", $ARGV[0] );
open( my $TSV,  "<", $ARGV[1] );

my %list;
while (<$LIST>) {
    chomp;
    my @tmp = split;
    $list{ $tmp[0] } = join( "\t", @tmp[ 1 .. $#tmp ] );
}
close($LIST);

while (<$TSV>) {
    chomp;
    my @tmp = split;
    if ( exists( $list{ $tmp[0] } ) ) {
        warn("$tmp[0]\t$list{ $tmp[0] }\n");
        print( join( "\t", @tmp ) );
        print("\n");
    }
}
close($TSV);

__END__
