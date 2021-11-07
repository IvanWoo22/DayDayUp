#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $TSV1, "<", $ARGV[0] );
open( my $TSV2, "<", $ARGV[1] );

my %table1;
my %table2;
my %table3;
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
}
close($TSV2);

my %rank;
my $rank0 = 1;
foreach my $t ( sort { $table3{$b} <=> $table3{$a} } keys %table3 ) {
    $rank{$t} = $rank0;
    $rank0++;
}

foreach my $t ( sort { $a <=> $b } keys %table3 ) {
    my @tbl  = split( /\s+/, $table2{$t} );
    my $name = $tbl[0];
    if ( exists( $table1{$name} ) ) {
        print( join( "\t", @tbl ) );
        print("\t$rank{$t}\t1\n");
    }
    else {
        print( join( "\t", @tbl ) );
        print("\t$rank{$t}\n");
    }
}

__END__
