#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my %table1;
my $col_num1;

open( my $TSV1, "<", $ARGV[0] );
foreach (<$TSV1>) {
    chomp;
    my @tbl  = split /\s+/;
    my $name = "$tbl[0]\t$tbl[1]";
    $col_num1 = $#tbl;
    $table1{$name} = join( "\t", @tbl[ 2 .. $col_num1 ] );
}
close($TSV1);

foreach my $FILE_INDEX ( 1 .. $#ARGV ) {
    my %table2;
    my $col_num2;
    open( my $TSV2, "<", $ARGV[$FILE_INDEX] );
    while (<$TSV2>) {
        chomp;
        my @tbl  = split /\s+/;
        my $name = "$tbl[0]\t$tbl[1]";
        $col_num2 = $#tbl;
        $table2{$name} = join( "\t", @tbl[ 2 .. $col_num2 ] );
    }
    close($TSV2);
    my %count;
    my @uniq = grep { ++$count{$_} < 2 } ( keys(%table1), keys(%table2) );
    foreach my $key (@uniq) {
        if ( ( exists( $table1{$key} ) ) and ( exists( $table2{$key} ) ) ) {
            $table1{$key} .= "\t" . $table2{$key};
        }
        elsif ( exists( $table2{$key} ) ) {
            $table1{$key} =
              join( "\t", ("NA") x ( $col_num1 - 1 ) ) . "\t" . $table2{$key};
        }
        elsif ( exists( $table1{$key} ) ) {
            $table1{$key} .= "\t" . join( "\t", ("NA") x ( $col_num2 - 1 ) );
        }
    }
}

for my $key ( keys %table1 ) {
    print("$key\t$table1{$key}\n");
}

__END__
