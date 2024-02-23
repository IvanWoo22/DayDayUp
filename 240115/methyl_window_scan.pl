#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;
use POSIX;

open my $BED, "<", $ARGV[0];

my %chr_count;
my %chr_methyl;
while (<$BED>) {
    chomp;
    my @line = split "\t";
    if ( exists( $chr_count{ $line[0] } ) ) {
        ${ $chr_count{ $line[0] } }[ ceil( ( $line[1] + 0.5 ) / 500 ) ]++;
        ${ $chr_count{ $line[0] } }[ floor( ( $line[1] + 0.5 ) / 500 ) ]++;
        ${ $chr_methyl{ $line[0] } }[ ceil( ( $line[1] + 0.5 ) / 500 ) ] +=
          $line[3];
        ${ $chr_methyl{ $line[0] } }[ floor( ( $line[1] + 0.5 ) / 500 ) ] +=
          $line[3];
    }
    else {
        ${ $chr_count{ $line[0] } }[ ceil( ( $line[1] + 0.5 ) / 500 ) ]  = 1;
        ${ $chr_count{ $line[0] } }[ floor( ( $line[1] + 0.5 ) / 500 ) ] = 1;
        ${ $chr_methyl{ $line[0] } }[ ceil( ( $line[1] + 0.5 ) / 500 ) ] =
          $line[3];
        ${ $chr_methyl{ $line[0] } }[ floor( ( $line[1] + 0.5 ) / 500 ) ] =
          $line[3];
    }
}

foreach my $chr ( 1 .. 5 ) {
    foreach my $window ( 0 .. $#{ $chr_count{$chr} } ) {
        if ( defined( ${ $chr_count{$chr} }[$window] )
            and ( ${ $chr_count{$chr} }[$window] > 3 ) )
        {
            my $start          = $window * 500;
            my $end            = $start + 1000;
            my $count          = ${ $chr_count{$chr} }[$window];
            my $average_methyl = ${ $chr_methyl{$chr} }[$window] / $count;
            print("$chr\t$start\t$end\t$count\t$average_methyl\n");
        }
    }
}

__END__
