#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

my %chr_ceil;
my %chr_floor;

open( my $CTR, "<", $ARGV[0] );
while (<$CTR>) {
    chomp;
    my @tmp = split;
    if ( exists( $chr_ceil{ $tmp[0] } ) ) {
        if ( $tmp[1] < $chr_floor{ $tmp[0] } ) {
            $chr_floor{ $tmp[0] } = $tmp[1];
        }
        if ( $tmp[2] > $chr_ceil{ $tmp[0] } ) {
            $chr_ceil{ $tmp[0] } = $tmp[2];
        }
    }
    else {
        $chr_floor{ $tmp[0] } = $tmp[1];
        $chr_ceil{ $tmp[0] }  = $tmp[2];
    }
}
close($CTR);

my %ctr_range;
foreach my $chr ( keys(%chr_ceil) ) {
    $ctr_range{$chr} = AlignDB::IntSpan->new;
    $ctr_range{$chr}
      ->AlignDB::IntSpan::add_range( $chr_floor{$chr}, $chr_ceil{$chr} );
}

open( my $LIST, "<", $ARGV[1] );
while (<$LIST>) {
    chomp;
    my @tmp = split;
    unless ( $ctr_range{ $tmp[0] }
        ->AlignDB::IntSpan::contains_any( $tmp[1], $tmp[2] ) )
    {
        print( join( "\t", @tmp ) );
        print("\n");
    }
}
close($LIST);

__END__
