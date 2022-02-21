#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $LIST, "<", $ARGV[0] );
open( my $TSV,  "<", $ARGV[1] );
readline($LIST);
my %list;
while (<$LIST>) {
    chomp;
    my @tmp = split;
    $list{ $tmp[1] . "\t" . $tmp[2] } = $tmp[0];
}
close($LIST);

while (<$TSV>) {
    chomp;
    my @tmp = split;
    if ( exists( $list{ $tmp[0] . "\t" . $tmp[1] } ) ) {
        my $out = $tmp[0] . "\t" . $tmp[1];
        print("$list{$out}\t");
        print( join( "\t", @tmp ) );
        print("\n");
    }
}
close($TSV);

__END__
