#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use open ':std', ':encoding(UTF-8)';

while (<>) {
    s/[\x{ff10}-\x{ff19}]/convert_to_halfwidth($&)/ge;
    print;
}

sub convert_to_halfwidth {
    my $char = shift;
    my $ord = ord($char);
    my $halfwidth = chr($ord - 0xfee0);
    return $halfwidth;
}

__END__