#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my %LIST;
open my $LIST, "<", $ARGV[0];
while ( my $line = <$LIST> ) {
    chomp $line;
    my @sample = split("\t",$line);
    $LIST{$sample[1]}="$sample[0].$sample[1]";
}
close($LIST);
chomp(my $sample_raw=$ARGV[1]);
my @sample = split(",",$sample_raw);
foreach(1..$#sample){
    print(",$LIST{$sample[$_]}");
}
print("\n");

__END__