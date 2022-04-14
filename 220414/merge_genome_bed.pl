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

my %unatac_chr;
open( my $ATAC, "<", $ARGV[2] );
while (<$ATAC>) {
    chomp;
    my @tmp = split;
    if ( exists( $unatac_chr{ $tmp[0] } ) ) {
        $unatac_chr{ $tmp[0] }->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
    }
    else {
        $unatac_chr{ $tmp[0] } = AlignDB::IntSpan->new;
        $unatac_chr{ $tmp[0] }->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
    }
}
close($ATAC);

my ( %eu_atac_chr, %het_atac_chr, %eu_unatac_chr, %het_unatac_chr );
foreach my $chr ( keys(%het_chr) ) {
    $eu_atac_chr{$chr}    = AlignDB::IntSpan->new;
    $het_atac_chr{$chr}   = AlignDB::IntSpan->new;
    $eu_unatac_chr{$chr}  = AlignDB::IntSpan->new;
    $het_unatac_chr{$chr} = AlignDB::IntSpan->new;
    $eu_atac_chr{$chr}    = $chr{$chr}->AlignDB::IntSpan::diff( $het_chr{$chr} )
      ->AlignDB::IntSpan::diff( $unatac_chr{$chr} );
    $het_atac_chr{$chr} =
      $chr{$chr}->AlignDB::IntSpan::diff( $unatac_chr{$chr} )
      ->AlignDB::IntSpan::intersect( $het_chr{$chr} );
    $eu_unatac_chr{$chr} = $chr{$chr}->AlignDB::IntSpan::diff( $het_chr{$chr} )
      ->AlignDB::IntSpan::intersect( $unatac_chr{$chr} );
    $het_unatac_chr{$chr} =
      $chr{$chr}->AlignDB::IntSpan::intersect( $het_chr{$chr} )
      ->AlignDB::IntSpan::intersect( $unatac_chr{$chr} );
    my ( $size1, $size2, $size3, $size4 ) = (
        $eu_atac_chr{$chr}->AlignDB::IntSpan::size,
        $het_atac_chr{$chr}->AlignDB::IntSpan::size,
        $eu_unatac_chr{$chr}->AlignDB::IntSpan::size,
        $het_unatac_chr{$chr}->AlignDB::IntSpan::size
    );
    print(
"$chr\n$size1: $eu_atac_chr{$chr}\n$size2: $het_atac_chr{$chr}\n$size3: $eu_unatac_chr{$chr}\n$size4: $het_unatac_chr{$chr}\n\n"
    );
}

__END__
