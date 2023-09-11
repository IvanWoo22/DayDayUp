#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;

my %chr_gene_range;
open my $GENE_RANGE, "<", $ARGV[0];
while (<$GENE_RANGE>) {
    chomp;
    my @tmp = split( "\t", $_ );
    ${ $chr_gene_range{ $tmp[0] } }{ $tmp[3] } = AlignDB::IntSpan->new;
    ${ $chr_gene_range{ $tmp[0] } }{ $tmp[3] }
      ->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
}
close $GENE_RANGE;

my %gene_promoter;
open my $PROMOTER_RANGE, "<", $ARGV[1];
while (<$PROMOTER_RANGE>) {
    chomp;
    my @tmp = split( "\t", $_ );
    $gene_promoter{ $tmp[3] } = AlignDB::IntSpan->new;
    $gene_promoter{ $tmp[3] }->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
}
close $PROMOTER_RANGE;

my %cov = ( "I" => 1, "II" => 2, "III" => 3, "IV" => 4, "V" => 5 );
open my $RANGE_LIST, "<", $ARGV[2];
while (<$RANGE_LIST>) {
    chomp;
    my @tmp = split( " ", $_ );
    my ( $gene_list1, $gene_list2 ) = ( "", "" );
    my $set1 = AlignDB::IntSpan->new;
    $set1->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
    my $oset1 = AlignDB::IntSpan->new;
    my $set2  = AlignDB::IntSpan->new;
    $set2->AlignDB::IntSpan::add_pair( $tmp[4], $tmp[5] );
    my $oset2 = AlignDB::IntSpan->new;

    foreach my $gene ( keys( %{ $chr_gene_range{ $cov{ $tmp[0] } } } ) ) {
        if (
            (
                $set1->AlignDB::IntSpan::overlap(
                    ${ $chr_gene_range{ $cov{ $tmp[0] } } }{$gene}
                ) >= 0.3 * ${ $chr_gene_range{ $cov{ $tmp[0] } } }{$gene}
                ->AlignDB::IntSpan::size
            )
            or (
                $set1->AlignDB::IntSpan::overlap(
                    ${ $chr_gene_range{ $cov{ $tmp[0] } } }{$gene}
                ) >= 0.8 * $set1->AlignDB::IntSpan::size
            )
          )
        {
            $gene_list1 = $gene_list1 . $gene . ",";
            $oset1->AlignDB::IntSpan::merge( $gene_promoter{$gene} );
        }
    }
    foreach my $gene ( keys( %{ $chr_gene_range{ $cov{ $tmp[3] } } } ) ) {
        if (
            (
                $set2->AlignDB::IntSpan::overlap(
                    ${ $chr_gene_range{ $cov{ $tmp[3] } } }{$gene}
                ) >= 0.3 * ${ $chr_gene_range{ $cov{ $tmp[3] } } }{$gene}
                ->AlignDB::IntSpan::size
            )
            or (
                $set2->AlignDB::IntSpan::overlap(
                    ${ $chr_gene_range{ $cov{ $tmp[3] } } }{$gene}
                ) >= 0.8 * $set2->AlignDB::IntSpan::size
            )
          )
        {
            $gene_list2 = $gene_list2 . $gene . ",";
            $oset2->AlignDB::IntSpan::merge( $gene_promoter{$gene} );
        }
    }
    print( join( "\t", @tmp ) );
    print "\t$gene_list1\t$oset1\t$gene_list2\t$oset2\n";
}
close $RANGE_LIST;
