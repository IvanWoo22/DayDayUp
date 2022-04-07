#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

my %chr_N;
open( my $CTR, "<", $ARGV[0] );
while (<$CTR>) {
    chomp;
    my @tmp = split;
    if ( exists( $chr_N{ $tmp[0] } ) ) {
        $chr_N{ $tmp[0] }->AlignDB::IntSpan::add_range( $tmp[1], $tmp[2] );
    }
    else {
        $chr_N{ $tmp[0] } = AlignDB::IntSpan->new;
        $chr_N{ $tmp[0] }->AlignDB::IntSpan::add_range( $tmp[1], $tmp[2] );
    }
}
close($CTR);

open( my $LIST, "<", $ARGV[1] );
while (<$LIST>) {
    chomp;
    my @tmp = split;
    unless (
        $chr_N{ $tmp[0] }->AlignDB::IntSpan::contains_any( $tmp[1], $tmp[2] ) )
    {
        print( join( "\t", @tmp ) );
        print("\n");
    }
}
close($LIST);

__END__
