#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open my $CSV, "<", $ARGV[0];
my $title = <$CSV>;
chomp $title;
my @title = split(",",$title);
my %id;
foreach my $id (1..$#title){
    my ($kind, undef) = split(/\./,$title[$id]);
    if ($kind eq $ARGV[1]) {
        $id{$id} = 1;
        print(",$title[$id]");
    }
}
print("\n");

while ( my $line = <$CSV> ) {
    chomp $line;
    my @mV = split(",",$line);
    print("$mV[0]");
    foreach my $id (1..$#mV){
        if (exists($id{$id})) {
            print(",$mV[$id]");
        }
    }
    print("\n");
}
close($CSV);

__END__