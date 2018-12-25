#!/usr/bin/perl -w
use strict;

open IN , "<", $ARGV[0];
chomp(my $tit1 = <IN>);
chomp(my $tit2 = <IN>);
chomp(my $tit3 = <IN>);

while (<IN>){
	chomp;
	my @tmp = split /\t/, $_;
	if ($tmp[6] eq "="){
		my @list = ($tmp[1], $tmp[3], $tmp[5], @tmp[7..9]);
		my $output = join("\t", @list);
		print "$output\n";
	}
}

close IN;
