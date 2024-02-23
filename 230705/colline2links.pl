#!/usr/bin/env perl
use strict;
use warnings;
use AlignDB::IntSpan;

my %genes;
my %genee;
open GFF,"<",$ARGV[0];
while(<GFF>){
        chomp;
        my @tmp = split("\t",$_);
        $genes{$tmp[1]} = $tmp[2];
        $genee{$tmp[1]} = $tmp[3];
}
close GFF;

open COLL,"<",$ARGV[1];
my @collin = <COLL>;
close COLL;

print"---\n";

my %chr = ( 1=>"I", 2=>"II", 3=>"III", 4=>"IV", 5=>"V" );

my @seg;
foreach(1..5){
        $seg[$_] = AlignDB::IntSpan->new;
}

my $i = 0;
while($i <= $#collin){
        if($collin[$i] =~ /## Alignment (\w+)/){
                $i++;
                my @at;
                until(($i > $#collin) or ($collin[$i] =~ /## Alignment/)){
                        my @tmp = split("\t",$collin[$i]);
                        push(@at,$tmp[1]);
                        $i++;
                }
                $i--;
                if($at[0] =~ /AT([1-5])G/){
                        my $chr = $1;
                        my $s = $genes{$at[0]} - 2000;
                        my $e = $genee{$at[-1]} + 2000;
                        $seg[$chr]->add_pair($s,$e);
                }
        }
        $i++;
}

foreach(1..5){
    print"$chr{$_}: $seg[$_]\n";
}