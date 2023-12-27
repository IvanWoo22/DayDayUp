#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;

my @chr_set;
foreach ( 1 .. 5 ) {
    $chr_set[$_] = AlignDB::IntSpan->new;
}

open my $BED, "<", $ARGV[0];
while (<$BED>) {
    chomp;
    my @tmp = split( "\t", $_ );
    if ( $tmp[0] =~ /^\d+$/ ) {
        my $chr   = $tmp[0];
        my $start = $tmp[1];
        my $end   = $tmp[2];
        $chr_set[$chr]->AlignDB::IntSpan::add_pair( $start, $end );
    }
}
close($BED);

my %cov = ( "I" => 1, "II" => 2, "III" => 3, "IV" => 4, "V" => 5 );

open my $SEG, "<", $ARGV[1];
while (<$SEG>) {
    chomp;
    my @tmp  = split( " ", $_ );
    my $set1 = AlignDB::IntSpan->new;
    my $set2 = AlignDB::IntSpan->new;
    $set1->add_pair( $tmp[1], $tmp[2] );
    $set2->add_pair( $tmp[4], $tmp[5] );
    my $intsec1 = $set1->intersect( $chr_set[ $cov{ $tmp[0] } ] );
    my $intsec2 = $set2->intersect( $chr_set[ $cov{ $tmp[3] } ] );
    my ( $e1, $e2, $p1, $p2 ) = (0,0,0,0);

    if ( $intsec1->cardinality >= 5 ) {
        $e1 = 1;
    }
    if ( $intsec2->cardinality >= 5 ) {
        $e2 = 1;
    }
    $p1 = $intsec1->cardinality / $set1->cardinality;
    $p2 = $intsec2->cardinality / $set2->cardinality;
    print "$_\t$e1\t$e2\t$p1\t$p2\n";
}
close($SEG);

__END__
