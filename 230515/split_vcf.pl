#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use List::Util;

while (<>) {
    next if s/^#//;
    chomp;
    my %record;
    my @info = split;
    if (($info[9] ne ".") && ($info[10] ne ".")) {
        for my $field (split /;/, $info[7]) {
            my ($key, $raw_value) = split /=/, $field;
            @{$record{$key}} = split /,/, $raw_value;
        }
        my (undef, undef, $N_AO, $N_DP, undef, undef, undef, undef, undef) =
            split /:/, $info[9];
        my (undef, undef, $T_AO, $T_DP, undef, undef, undef, undef, undef) =
            split /:/, $info[10];
        my @N_AO = split /,/, $N_AO;
        my @T_AO = split /,/, $T_AO;
        if (grep(/ins|del/, @{$record{"TYPE"}})
            && ($T_DP > 10)
            && ($N_DP > 10)) {
            foreach my $no (0 .. $#{$record{"TYPE"}}) {
                if (${$record{"TYPE"}}[$no] eq "ins"
                    || ${$record{"TYPE"}}[$no] eq "del") {
                    if (($N_AO[$no] == 0) && ($T_AO[$no] > 0)) {
                        print(join("\t", @info[ 0, 1, 3, 4, 7, 9, 10 ]));
                        print("\n");
                    }
                }
            }
        }
    }
}
__END__
