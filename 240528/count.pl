#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub SEQ_A_PERCENT {
    my $seq   = shift;
    my $count = () = $seq =~ /A/g;
    my $out   = $count / length($seq);
    return ($out);
}

my ( %set_count, %pir_count );
open( my $COUNT, "<", $ARGV[0] );

while (<$COUNT>) {
    chomp;
    my ( undef, $count, $set, $pir, undef, $tail ) = split( /\t/, $_ );
    if ( $tail eq "" ) {
        if ( exists( $set_count{$set} ) ) {
            $set_count{$set}->{"nona"} += $count;
        }
        else {
            $set_count{$set}->{"nona"} = $count;
            $set_count{$set}->{"a"}    = 0;
        }
        if ( exists( $pir_count{$pir} ) ) {
            $pir_count{$pir}->{"nona"} += $count;
        }
        else {
            $pir_count{$pir}->{"nona"} = $count;
            $pir_count{$pir}->{"a"}    = 0;
        }
    }
    elsif (( ( length($tail) > 4 ) and ( SEQ_A_PERCENT($tail) >= 0.8 ) )
        or ( $tail eq "AAA" )
        or ( $tail eq "AAAA" ) )
    {
        if ( exists( $set_count{$set} ) ) {
            $set_count{$set}->{"a"} += $count;
        }
        else {
            $set_count{$set}->{"nona"} = 0;
            $set_count{$set}->{"a"}    = $count;
        }
        if ( exists( $pir_count{$pir} ) ) {
            $pir_count{$pir}->{"a"} += $count;
        }
        else {
            $pir_count{$pir}->{"nona"} = 0;
            $pir_count{$pir}->{"a"}    = $count;
        }
    }

}
close($COUNT);

open( my $SET, ">", "$ARGV[1].set.tsv" );
foreach my $id ( sort keys(%set_count) ) {
    print $SET "$id\t$set_count{$id}->{nona}\t$set_count{$id}->{a}\n";
}
close($SET);

open( my $PIR, ">", "$ARGV[1].pir.tsv" );
foreach my $id ( sort keys(%pir_count) ) {
    print $PIR "$id\t$pir_count{$id}->{nona}\t$pir_count{$id}->{a}\n";
}
close($PIR);

__END__
