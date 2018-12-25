#!/usr/bin/perl -w
use strict;

open IN1 ,"<", $ARGV[0];
open IN2 ,"<", $ARGV[1];
open OUT ,">", $ARGV[2];


my @sitecount;


while (my $r1 = <IN2>){
	chomp($r1);
	chomp(my $r2 = <IN2>);
	my @tmp1 = split /\t/, $r1;
	my @tmp2 = split /\t/, $r2;
	if ( $tmp1[4] < 0 ){
		my $sit = $tmp1[3]-$tmp1[4]-2;
		if (exists $sitecount[$sit]){
			$sitecount[$sit]++;
		}else{
			$sitecount[$sit]=1;
		}
	}elsif( $tmp2[4] < 0 ){
		my $sit = $tmp2[3]-$tmp2[4]-1;
		if (exists $sitecount[$sit]){
			$sitecount[$sit]++;
		}else{
			$sitecount[$sit]=1;
		}
	}
}

close IN2;

chomp(my $ref = <IN1>);
my @refbas = split //, $ref;
my $len = $#refbas;
foreach my $i (0..$len-14) {
	my @tmp = @refbas[$i..$i+14];
	my $barc = join('',@tmp);
	if (exists $sitecount[$i+14]){
                print OUT "$barc\t$sitecount[$i+14]\n";
        }else{
                print OUT "$barc\tNo Match.\n";
        }
}



close IN1;
close OUT;
