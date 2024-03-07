#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;

open my $GFF, "<", $ARGV[0];
my $BED_SUFFIX = $ARGV[1];
open my $BED_PMT,        ">", $BED_SUFFIX . "_promoter.bed";
open my $BED_5UTR,       ">", $BED_SUFFIX . "_5UTR.bed";
open my $BED_3UTR,       ">", $BED_SUFFIX . "_3UTR.bed";
open my $BED_CDS,        ">", $BED_SUFFIX . "_CDS.bed";
open my $BED_GENE,       ">", $BED_SUFFIX . "_gene.bed";
open my $YML_INTERGENIC, ">", $BED_SUFFIX . "_intergenic.yml";
open my $YML_INTRON,     ">", $BED_SUFFIX . "_intron.yml";
my $longest = 0;
my $chr;
my $gene_range;
my %chr_intron_range = (
    1 => AlignDB::IntSpan->new,
    2 => AlignDB::IntSpan->new,
    3 => AlignDB::IntSpan->new,
    4 => AlignDB::IntSpan->new,
    5 => AlignDB::IntSpan->new
);
my %chr_intergenic_range = (
    1 => AlignDB::IntSpan->new,
    2 => AlignDB::IntSpan->new,
    3 => AlignDB::IntSpan->new,
    4 => AlignDB::IntSpan->new,
    5 => AlignDB::IntSpan->new
);
$chr_intergenic_range{1}->AlignDB::IntSpan::add_pair( 1, 30427671 );
$chr_intergenic_range{2}->AlignDB::IntSpan::add_pair( 1, 19698289 );
$chr_intergenic_range{3}->AlignDB::IntSpan::add_pair( 1, 23459830 );
$chr_intergenic_range{4}->AlignDB::IntSpan::add_pair( 1, 18585056 );
$chr_intergenic_range{5}->AlignDB::IntSpan::add_pair( 1, 26975502 );

while (<$GFF>) {
    if ( ( !/^##/ ) and ( !/^Chr[MC]/ ) ) {
        chomp;
        my @tmp = split "\t";
        if ( $tmp[2] eq "mRNA" ) {
            if ( $tmp[8] =~ /longest=1/ ) {
                $longest = 1;
                if ( $tmp[8] =~ /ID=([^;]+);/ ) {
                    my $cid = $1;
                    $tmp[0] =~ s/Chr//;
                    $chr = $tmp[0];
                    my ( $start, $end );
                    if ( $tmp[6] eq "+" ) {
                        $start = $tmp[3] - 2000;
                        $start = 1 if $start < 1;
                        $end   = $tmp[3] - 1;
                        print $BED_GENE "$chr\t$start\t$tmp[4]\t$cid\n";
                    }
                    else {
                        $start = $tmp[4] + 1;
                        $end   = $tmp[4] + 2000;
                        print $BED_GENE "$chr\t$tmp[3]\t$end\t$cid\n";
                    }
                    print $BED_PMT "$chr\t$start\t$end\t$cid\n";
                    $gene_range = AlignDB::IntSpan->new;
                    $gene_range->AlignDB::IntSpan::add_pair( $tmp[3], $tmp[4] );
                    if ( exists( $chr_intergenic_range{$chr} ) ) {
                        $chr_intergenic_range{$chr}
                          ->AlignDB::IntSpan::remove_range( $tmp[3], $tmp[4] );
                    }
                }
            }
            else {
                $longest = 0;
                $chr_intron_range{$chr}->AlignDB::IntSpan::merge($gene_range);
            }
        }
        elsif ( $longest == 1 ) {
            if ( $tmp[2] eq "five_prime_UTR" ) {
                if ( $tmp[8] =~ /Parent=([^;]+);/ ) {
                    my $cid = $1;
                    $tmp[0] =~ s/Chr//;
                    my ( $start, $end );
                    $start = $tmp[3];
                    $end   = $tmp[4];
                    print $BED_5UTR "$tmp[0]\t$start\t$end\t$cid\n";
                }
            }
            elsif ( $tmp[2] eq "three_prime_UTR" ) {
                if ( $tmp[8] =~ /Parent=([^;]+);/ ) {
                    my $cid = $1;
                    $tmp[0] =~ s/Chr//;
                    my ( $start, $end );
                    $start = $tmp[3];
                    $end   = $tmp[4];
                    print $BED_3UTR "$tmp[0]\t$start\t$end\t$cid\n";
                }
            }
            elsif ( $tmp[2] eq "CDS" ) {
                if ( $tmp[8] =~ /Parent=([^;]+);/ ) {
                    my $cid = $1;
                    $tmp[0] =~ s/Chr//;
                    my ( $start, $end );
                    $start = $tmp[3];
                    $end   = $tmp[4];
                    print $BED_CDS "$tmp[0]\t$start\t$end\t$cid\n";
                }
            }
            elsif ( $tmp[2] eq "exon" ) {
                $gene_range->AlignDB::IntSpan::remove_range( $tmp[3], $tmp[4] );
            }
        }
    }
}
close($GFF);
close($BED_PMT);
close($BED_5UTR);
close($BED_3UTR);
close($BED_CDS);
print $YML_INTERGENIC ("---\n");
foreach ( 1 .. 5 ) {
    print $YML_INTERGENIC ( $_ . ": "
          . $chr_intergenic_range{$_}->AlignDB::IntSpan::as_string
          . "\n" );
}
print $YML_INTRON ("---\n");
foreach ( 1 .. 5 ) {
    print $YML_INTRON (
        $_ . ": " . $chr_intron_range{$_}->AlignDB::IntSpan::as_string . "\n" );
}
close($YML_INTERGENIC);
close($YML_INTRON);

__END__
