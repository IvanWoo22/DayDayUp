#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my %dict_motif1;
my @element = ( "A", "T", "C", "G" );
for ( my $i = 0 ; $i < 10 ; $i++ ) {
    if ( $i == 0 ) {
        foreach my $base (@element) {
            $dict_motif1{$base} = 0;
        }
    }
    else {
        foreach my $old ( keys %dict_motif1 ) {
            foreach my $base (@element) {
                my $new = $old . $base;
                $dict_motif1{$new} = 0;
            }
            delete $dict_motif1{$old};
        }
    }
}

my %dict_motif2 = %dict_motif1;

open( my $NC, "<", $ARGV[0] );
while (<$NC>) {
    chomp;
    my @tmp = split( "", $_ );
    $dict_motif1{ join( "", @tmp[ 6 .. 15 ] ) }++;
}
close($NC);

open( my $TR, "<", $ARGV[1] );
while (<$TR>) {
    chomp;
    my @tmp = split( "", $_ );
    $dict_motif2{ join( "", @tmp[ 6 .. 15 ] ) }++;
}
close($TR);

foreach my $motif ( sort keys %dict_motif1 ) {
    print "$motif\t$dict_motif1{$motif}\t$dict_motif2{$motif}\n";
}

__END__
