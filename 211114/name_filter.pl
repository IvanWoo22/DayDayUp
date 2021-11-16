#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open( my $TSV1, "<", $ARGV[0] );
open( my $TSV2, "<", $ARGV[1] );

my %table;

while (<$TSV1>) {
    chomp;
    $table{$_} = 1;
}
close($TSV1);

while (<$TSV2>) {
    chomp;
    my $line = $_;
    my @tmp  = split /,/, $line;
    $tmp[0] =~ s/"//g;
    if ( exists( $table{ $tmp[0] } ) ) {
        $table{$tmp[0]} ++;
        print("$line\n");
    }
}
close($TSV2);

foreach(keys(%table)){
    if ($table{$_}<2) {
        warn("$_\n");
    }
}

__END__
