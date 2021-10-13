#!/usr/bin/perl -w
use strict;
use AlignDB::IntSpan;

open IN1, "<", $ARGV[0];
open IN2, "<", $ARGV[1];
open OUT, ">", $ARGV[2];

my $line = 0;
my %loci;
my @start;
my @end;
my @info;
$start[0] = 0;
while (my $gffline = <IN1>) {
    $line++;
    chomp($gffline);
    my @gfftmp = split /\t/, $gffline;
    $info[$line] = "chr$gfftmp[0]\t$gfftmp[3]-$gfftmp[4]\t$gfftmp[8]";
    $start[$line] = $gfftmp[3];
    $end[$line] = $gfftmp[4];
    $loci{$gfftmp[3]} = $line;
}

while (my $cpgline = <IN2>) {
    chomp($cpgline);
    my @cpgtmp = split /\t/, $cpgline;
    my $cpg = AlignDB::IntSpan->new;
    $cpg->add_range($cpgtmp[0], $cpgtmp[1]);

    my $score = 0;
    my $mark = 0;
    foreach (1 .. $line) {
        my $gff = AlignDB::IntSpan->new;
        $gff->add_range($start[$_], $end[$_]);
        my $merge = $gff->intersect($cpg);
        if ($merge->cardinality > $score) {
            $mark = $_;
            $score = $merge->cardinality;
        }
    }
    if ($mark == 0) {
        my $gene = AlignDB::IntSpan->new(@start);
        my $island = $gene->nearest_island($cpgtmp[0]);
        $mark = $loci{$island};
    }
    print OUT "$info[$mark]\t$cpgtmp[0]-$cpgtmp[1]\t$cpgtmp[3]/$cpgtmp[2]\t$cpgtmp[4]\t$cpgtmp[5]\t$cpgtmp[6]\t$cpgtmp[7]\n" || die "Error in $ARGV[0] in $cpgline: $!";
}
close IN1;
close IN2;
close OUT;
