#!/usr/bin/env perl
use strict;
use warnings;

my %feature;
open( my $TIME, "<", $ARGV[0] );
while (<$TIME>) {
    chomp;
    my @tmp = split "\t";
    $feature{ $tmp[0] } = 1;
}
close($TIME);

open( my $LINE, "<", $ARGV[1] );
while (<$LINE>) {
    chomp;
    my @tmp  = split "\t";
    my $time = $tmp[2] . $tmp[3] . $tmp[4] . $tmp[5] . $tmp[6] . $tmp[7];
    if ( exists( $feature{$time} ) ) {
        print "$tmp[0]\t$tmp[1]\n";
    }
}
close($LINE);

__END__
