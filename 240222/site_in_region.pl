#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;

my %chr_set;
open my $R_BED, "<", $ARGV[0];
while (<$R_BED>) {
    chomp;
    my @tmp = split( "\t", $_ );
    my $chr = $tmp[0];
    unless ( exists( $chr_set{$chr} ) ) {
        $chr_set{$chr} = AlignDB::IntSpan->new;
    }
    my $start = $tmp[1];
    my $end   = $tmp[2];
    $chr_set{$chr}->AlignDB::IntSpan::add_pair( $start, $end );
}
close($R_BED);

open my $S_BED, "<", $ARGV[1];
while (<$S_BED>) {
    chomp;
    my @tmp = split( "\t", $_ );
    if ( exists( $chr_set{ $tmp[0] } )
        and ( $chr_set{ $tmp[0] }->AlignDB::IntSpan::contains_any( $tmp[1] ) ) )
    {
        print("$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\n");
    }
}
close($S_BED);

__END__