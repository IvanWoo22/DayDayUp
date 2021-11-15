#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my @info;
my %list;
open my $LS, "<", $ARGV[0];
while (<$LS>) {
    chomp;
    $list{$_} = 1;
}
close($LS);
open my $FA, "<", $ARGV[1];
@info = <$FA>;
close($FA);

my $i = 0;
while ( $i <= $#info ) {
    my ( undef, $gene_raw, undef ) = split /\|/, $info[$i];
    my ( $gene, undef ) = split /\./, $gene_raw;
    if ( exists( $list{$gene} ) ) {
        print( $info[$i] );
        my $j = 1;
        until ( ( $i + $j > $#info ) or ( $info[ $i + $j ] =~ /^>/ ) ) {
            print( $info[ $i + $j ] );
            $j++;
        }
        $i += $j;
    }
    else {
        my $j = 1;
        until ( ( $i + $j > $#info ) or ( $info[ $i + $j ] =~ /^>/ ) ) {
            $j++;
        }
        $i += $j;
    }
}

__END__
