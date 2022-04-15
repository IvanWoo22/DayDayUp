#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

my %chr;
open( my $CHR_RANGE, "<", $ARGV[0] );
while (<$CHR_RANGE>) {
    chomp;
    my @tmp = split;
    $chr{ $tmp[0] } = AlignDB::IntSpan->new;
    $chr{ $tmp[0] }->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
}
close($CHR_RANGE);
my %het_chr;
open( my $HET, "<", $ARGV[1] );
while (<$HET>) {
    chomp;
    my @tmp = split;
    if ( exists( $het_chr{ $tmp[0] } ) ) {
        $het_chr{ $tmp[0] }->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
    }
    else {
        $het_chr{ $tmp[0] } = AlignDB::IntSpan->new;
        $het_chr{ $tmp[0] }->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
    }
}
close($HET);

my (%eu_chr);
foreach my $chr ( keys(%het_chr) ) {
    $eu_chr{$chr} = AlignDB::IntSpan->new;
    $eu_chr{$chr} = $chr{$chr}->AlignDB::IntSpan::diff( $het_chr{$chr} );
}

open( my $SITES, "<", $ARGV[2] );
my ( $het, $eu ) = ( 0, 0 );
while (<$SITES>) {
    chomp;
    my ( $header, $pos ) = split "\t", $_;
    if ( $het_chr{$header}->AlignDB::IntSpan::contains($pos) ) {
        $het++;
    }
    elsif ( $eu_chr{$header}->AlignDB::IntSpan::contains($pos) ) {
        $eu++;
    }
}
close($SITES);
print "$het\n$eu\n";

__END__
