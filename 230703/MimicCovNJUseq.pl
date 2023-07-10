#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my %cov_count;
while (<>) {
    chomp;
    my ( $contig, $pos, $dir, undef, undef, undef, undef, undef, $end_count ) =
      split "\t";
    if ( $dir eq "+" ) {
        foreach ( $pos - 18 .. $pos ) {
            if ( exists( $cov_count{ $contig . $dir . $_ } ) ) {
                $cov_count{ $contig . $dir . $_ } += $end_count;
            }
            else {
                $cov_count{ $contig . $dir . $_ } = $end_count;
            }
        }
    }
    else {
        foreach ( $pos .. $pos + 18 ) {
            if ( exists( $cov_count{ $contig . $dir . $_ } ) ) {
                $cov_count{ $contig . $dir . $_ } += $end_count;
            }
            else {
                $cov_count{ $contig . $dir . $_ } = $end_count;
            }
        }
    }
}

foreach my $key ( sort { $cov_count{$b} <=> $cov_count{$a} } keys %cov_count ) {
    print("$key\t$cov_count{$key}\n");
}

__END__
