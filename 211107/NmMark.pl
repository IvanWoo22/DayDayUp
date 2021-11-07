#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $TSV1, "<", $ARGV[0] );
open( my $TSV2, "<", $ARGV[1] );

my ( %table1, %table2, %table3, %table4, %table5, %table6 );
while (<$TSV1>) {
    chomp;
    my @tbl = split /\s+/;
    $table1{ $tbl[0] } = 1;
}
close($TSV1);

while (<$TSV2>) {
    chomp;
    my @tbl  = split /\s+/;
    my $name = $tbl[0];
    $table2{$name} = join( "\t", @tbl );
    $table3{$name} = $tbl[-1];
    $table4{$name} = $tbl[-4];
    $table5{$name} = $tbl[-3];
    $table6{$name} = $tbl[-2];
}
close($TSV2);

my ( %ranks, %rank1, %rank2, %rank3 );
my $rank0 = 1;
foreach my $t ( sort { $table3{$b} <=> $table3{$a} } keys %table3 ) {
    $ranks{$t} = $rank0;
    $rank0++;
}

$rank0 = 1;
foreach my $t ( sort { $table4{$b} <=> $table4{$a} } keys %table4 ) {
    $rank1{$t} = $rank0;
    $rank0++;
}

$rank0 = 1;
foreach my $t ( sort { $table5{$b} <=> $table5{$a} } keys %table5 ) {
    $rank2{$t} = $rank0;
    $rank0++;
}

$rank0 = 1;
foreach my $t ( sort { $table6{$b} <=> $table6{$a} } keys %table6 ) {
    $rank3{$t} = $rank0;
    $rank0++;
}

foreach my $t ( sort { $a <=> $b } keys %table3 ) {
    my @tbl  = split( /\s+/, $table2{$t} );
    my $name = $tbl[0];
    if ( exists( $table1{$name} ) ) {
        print( join( "\t", @tbl ) );
        print("\t$rank1{$t}\t$rank2{$t}\t$rank3{$t}\t$ranks{$t}\t1\n");
    }
    else {
        print( join( "\t", @tbl ) );
        print("\t$rank1{$t}\t$rank2{$t}\t$rank3{$t}\t$ranks{$t}\n");
    }
}

__END__
