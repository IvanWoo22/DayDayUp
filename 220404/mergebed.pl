#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

my %chr_bed;
foreach my $BED_FILE (@ARGV) {
    open( my $CTR, "<", $BED_FILE );
    while (<$CTR>) {
        chomp;
        my @tmp = split;
        if ( exists( $chr_bed{ $tmp[0] } ) ) {
            $chr_bed{ $tmp[0] }
              ->AlignDB::IntSpan::add_range( $tmp[1], $tmp[2] );
        }
        else {
            $chr_bed{ $tmp[0] } = AlignDB::IntSpan->new;
            $chr_bed{ $tmp[0] }
              ->AlignDB::IntSpan::add_range( $tmp[1], $tmp[2] );
        }
    }
    close($CTR);
}
foreach my $key ( keys(%chr_bed) ) {
    my $out = $chr_bed{ $key }->AlignDB::IntSpan::ranges;
    print("$key\t${$out}[0]\n");
}
