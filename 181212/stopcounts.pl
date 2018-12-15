#!/usr/bin/perl -w
use strict;

open IN1,"<", $ARGV[0];
open IN2,"<", $ARGV[1];
open OUT,">", $ARGV[2];

my $map_len = 18;
my %bas;
chomp(my $raw = <IN1>);

while(<IN2>){
	chomp;
	my @tmp = split /\t/, $_;
	my @ls = split //, $tmp[2];
	my $base = join('',@ls[-$map_len..-1]);
	if ($raw =~ /$tmp[2]/){
		if (exists $bas{$base}){
			$bas{$base} += $tmp[1];
		}else{
			$bas{$base} = $tmp[1];
		}
	}
}


my @refbas = split //, $raw;
my $len = $#refbas;
foreach my $i (0..$len-$map_len+1) {
	my @tmp = @refbas[$i..$i+$map_len-1];
	my $barc = join('',@tmp);
	if (exists $bas{$barc}){
                print OUT "$barc\t$bas{$barc}\r\n";
        }else{
                print OUT "$barc\tNo Match.\r\n";
        }
}


close IN1;
close IN2;
close OUT;
