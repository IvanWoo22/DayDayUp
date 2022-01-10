#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use List::Util qw"max min";
use Statistics::Basic;

open( my $TSV, "<", $ARGV[0] );

readline($TSV);
my @dif;
while(<$TSV>){
    chomp;
    my @tmp  = split /\t/;
    my @val_ch;
    foreach my $i (1..$#tmp){
        if ($tmp[$i] ne "NA") {
            push(@val_ch,$tmp[$i]);
        }
    }
    if ($#val_ch > 2) {
        push(@dif,max(@val_ch)-min(@val_ch));
    }
}

close($TSV);

my $median = Statistics::Basic::median(@dif);
my $std = Statistics::Basic::stddev(@dif);
print "$#dif\t$median\t$std\n";

__END__