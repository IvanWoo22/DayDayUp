#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

sub cosine {
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

my ( %pv, %vector, @site_id );
my $site = AlignDB::IntSpan->new;
while (<STDIN>) {
    chomp;
    my @tmp = split;
    push( @site_id, $tmp[0] . "_" . $tmp[1] );
    $site->AlignDB::IntSpan::add( $tmp[1] );
    $pv{ $tmp[0] . "_" . $tmp[1] }     = $tmp[2];
    $vector{ $tmp[0] . "_" . $tmp[1] } = join( "\t", @tmp[ 4 .. 45 ] );
}

foreach my $id ( 0 .. $#site_id ) {
    my ( $chr, $position ) = split( "_", $site_id[$id] );
    my $left_set  = AlignDB::IntSpan->new;
    my $right_set = AlignDB::IntSpan->new;
    $left_set->AlignDB::IntSpan::add_range( $position - 100, $position - 1 );
    $right_set->AlignDB::IntSpan::add_range( $position + 1, $position + 100 );
    my @near_site =
      $site->AlignDB::IntSpan::nearest_island($position)->as_array;
    my $left_count = $site->AlignDB::IntSpan::intersect($left_set)->cardinality;
    my $right_count =
      $site->AlignDB::IntSpan::intersect($right_set)->cardinality;
    my $length       = abs( $position - $near_site[0] );
    my $eul_distance = cosine( $vector{ $site_id[$id] },
        $vector{ $chr . "_" . $near_site[0] } );
    print(
"$chr\t$position\t$pv{$site_id[$id]}\t$near_site[0]\t$eul_distance\t$length\t$left_count\t$right_count\t$vector{$site_id[$id]}\n"
    );
}

__END__
