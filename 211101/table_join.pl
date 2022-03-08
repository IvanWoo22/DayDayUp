#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $TSV1, "<", $ARGV[0] );
open( my $TSV2, "<", $ARGV[1] );

my %table1;
my %table2;
my $col_num1;
my $col_num2;

sub ADD_COL {
    my ( $TBL1, $TBL2, $WIDTH1, $WIDTH2, $NAME ) = @_;
    if ( ( exists( ${$TBL1}{$NAME} ) ) and ( exists( ${$TBL2}{$NAME} ) ) ) {
        push( @{ ${$TBL1}{$NAME} }, @{ ${$TBL2}{$NAME} } );
    }
    elsif ( exists( ${$TBL2}{$NAME} ) ) {
        @{ ${$TBL1}{$NAME} }[ 0 .. $WIDTH1 - 1 ] = ("N/A") x $WIDTH1;
        push( @{ ${$TBL1}{$NAME} }, @{ ${$TBL2}{$NAME} } );
    }
    elsif ( exists( ${$TBL1}{$NAME} ) ) {
        @{ ${$TBL1}{$NAME} }[ $WIDTH1 .. $WIDTH1 + $WIDTH2 - 1 ] =
          ("N/A") x $WIDTH2;
    }
    return ($TBL1);
}

while (<$TSV1>) {
    chomp;
    my @tbl  = split /\s+/;
    my $name = $tbl[0];
    $col_num1 = $#tbl;
    @{ $table1{$name} } = @tbl[ 1 .. $col_num1 ];
}
close($TSV1);

while (<$TSV2>) {
    chomp;
    my @tbl  = split /\s+/;
    my $name = $tbl[0];
    $col_num2 = $#tbl;
    @{ $table2{$name} } = @tbl[ 1 .. $col_num2 ];
}
close($TSV2);

my %count;
foreach my $e ( keys(%table1), keys(%table2) ) {
    $count{$e}++;
}

$col_num1 = $col_num1-0;
$col_num2 = $col_num2-0;

for my $key ( keys %count ) {
    if ( ( exists( $table1{$key} ) ) and ( exists( $table2{$key} ) ) ) {
        print("$key\t");
        print( join( "\t", @{ $table1{$key} } ) );
        print("\t");
        print( join( "\t", @{ $table2{$key} } ) );
        print("\n");
    }
    elsif ( exists( $table2{$key} ) ) {
        print("$key\t");
        print( join( "\t", ("N/A") x $col_num1 ) );
        print("\t");
        print( join( "\t", @{ $table2{$key} } ) );
        print("\n");
    }
    elsif ( exists( $table1{$key} ) ) {
        print("$key\t");
        print( join( "\t", @{ $table1{$key} } ) );
        print("\t");
        print( join( "\t", ("N/A") x $col_num2 ) );
        print("\n");
    }
}

__END__
