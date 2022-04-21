#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

open( my $SITES, "<", $ARGV[0] );
my %site_list;
while (<$SITES>) {
    chomp;
    my @tmp = split;
    if ( exists( $site_list{ $tmp[0] } ) ) {
        $site_list{ $tmp[0] }->AlignDB::IntSpan::add( $tmp[1] );
    }
    else {
        $site_list{ $tmp[0] } = AlignDB::IntSpan->new;
        $site_list{ $tmp[0] }->AlignDB::IntSpan::add( $tmp[1] );
    }
}
close($SITES);

open( my $REGIONS, "<", $ARGV[1] );
while (<$REGIONS>) {
    chomp;
    my ( $header, $start, $end ) = split;
    my $set = AlignDB::IntSpan->new;
    $set->AlignDB::IntSpan::add_pair( $start, $end );
    if ( $site_list{$header}->AlignDB::IntSpan::overlap($set) > 0 ) {
        my $count = $site_list{$header}->AlignDB::IntSpan::overlap($set);
        print("$header\t$start\t$end\t$count\n");
    }
    else {
        print("$header\t$start\t$end\t0\n");
    }
}
close($REGIONS);

__END__
