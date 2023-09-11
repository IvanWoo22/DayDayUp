#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;

open my $LIST,    "<", $ARGV[0];
open my $GFF,     "<", $ARGV[1];
open my $BED_CDS, ">", $ARGV[2];
open my $BED_PMT, ">", $ARGV[3];

my ( %id, %range );
while (<$LIST>) {
    chomp;
    $id{$_}    = 1;
    $range{$_} = AlignDB::IntSpan->new;
}
close($LIST);

while (<$GFF>) {
    if ( ( !/^##/ ) and ( !/^Chr[MC]/ ) ) {
        chomp;
        my @tmp = split "\t";
        if ( $tmp[2] eq "mRNA" ) {
            if ( $tmp[8] =~ /ID=(\w+)\.(\d+)\.Araport11\.447;/ ) {
                my $cid = $1 . "." . $2;
                if ( exists( $id{$cid} ) ) {
                    $tmp[0] =~ s/Chr//;
                    print $BED_CDS "$tmp[0]\t$tmp[3]\t$tmp[4]\t$cid\n";
                }
            }
        }
        elsif ( $tmp[2] eq "CDS" ) {
            if ( $tmp[8] =~ /ID=(\w+)\.(\d+)\.Araport11\.447\.CDS\.1;/ ) {
                my $cid = $1 . "." . $2;
                if ( exists( $id{$cid} ) ) {
                    $tmp[0] =~ s/Chr//;
                    my ( $start, $end );
                    if ( $tmp[6] eq "+" ) {
                        $start = $tmp[3] - 2000;
                        $end   = $tmp[3] - 1;
                    }
                    else {
                        $start = $tmp[4] + 1;
                        $end   = $tmp[4] + 2000;
                    }
                    print $BED_PMT "$tmp[0]\t$start\t$end\t$cid\n";
                }
            }
        }
    }
}
close($GFF);
close($BED_CDS);
close($BED_PMT);
