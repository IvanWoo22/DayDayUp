#!/usr/bin/perl -w
use strict;

open IN1,"<", $ARGV[0];
open IN2,"<", $ARGV[1];
my $OUT = "nor.$ARGV[0]";
open OUT1,">", $OUT;
open OUT2,">", $ARGV[2];

my $map_len = $ARGV[3];

my %bas;
while(<IN1>){
	chomp;
	my @tmp = split /\t/, $_;
	my @ls = split //, $tmp[2];
	my $base = join('',@ls[-$map_len..-1]);
	if (exists $bas{$base}){
		$bas{$base} += $tmp[1];
	}else{
		$bas{$base} = $tmp[1];
	}
}


foreach my $char (sort { $bas{$b} <=> $bas{$a} } keys %bas)
{
    print OUT1 "$char\t$bas{$char}\n";
}


my $raw = <IN2>;
chomp($raw);
my @refbas = split //, $raw;
my $len = $#refbas;
foreach my $i (0..$len-$map_len+1) {
	my @tmp = @refbas[$i..$i+$map_len-1];
	my $barc = join('',@tmp);
	if (exists $bas{$barc}){
                print OUT2 "$barc\t$bas{$barc}\r\n";
        }else{
                print OUT2 "$barc\tNo Match.\r\n";
        }
}


close IN1;
close IN2;
close OUT1;
close OUT2;
