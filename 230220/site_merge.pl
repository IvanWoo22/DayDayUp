#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my @chr;
my @range;
my @depth;

my $gap = $ARGV[0];
my $maxlength = $ARGV[1];
my $cpgcount = $ARGV[2];

while (<STDIN>) {
    chomp;
    my @tmp = split;
    push(@chr, $tmp[0]);
    push(@range, $tmp[1]);
    push(@depth, $tmp[2]);
}

my ($chr, $start_point, $end_point, $site_count) = (0, 0, 0, 0);
my @dep_count;
foreach my $id (0 .. $#chr) {
    if (($range[$id] <= $gap + $end_point) && ($start_point + $maxlength >= $range[$id]) && ($chr eq $chr[$id])) {
        $end_point = $range[$id];
        $site_count++;
        $dep_count[$depth[$id]]++;
    }
    else {
        if ($site_count > $cpgcount) {
            print "$chr\t$start_point\t$end_point\t$site_count\t";
            print(join("\t", @dep_count[1 .. 2]));
            print "\n";
        }
        $site_count = 1;
        $start_point = $range[$id];
        $end_point = $range[$id];
        $chr = $chr[$id];
        @dep_count[1, 2] = (0, 0);
        $dep_count[$depth[$id]]++;
    }
}

if ($site_count > $cpgcount) {
    print "$chr\t$start_point\t$end_point\t$site_count\t";
    print(join("\t", @dep_count[1 .. 2]));
    print "\n";
}

__END__