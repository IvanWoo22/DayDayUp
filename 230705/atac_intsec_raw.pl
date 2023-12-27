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
    my @tmp = split( " ", $_ );
    my $set = AlignDB::IntSpan->new;
    $set->add_pair( $tmp[1], $tmp[2] );
    my $intsec = $set->intersect( $chr_set[ $cov{ $tmp[0] } ] );
    my ( $e1, $e2 ) = ( $intsec->cardinality, $set->cardinality );
    my $p = $intsec->cardinality / $set->cardinality;
    print "$_\t$e1\t$e2\t$p\n";
}
close($SEG);

__END__
