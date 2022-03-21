#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

open(my $TSV, "<", $ARGV[0]);
open(my $LIST, "<", $ARGV[1]);

my %cgid;
my %chrlist;
while(<$LIST>){
    chomp;
    my @tmp = split;
    $cgid{$tmp[1]."_".$tmp[2]} = $tmp[0];
    push(@{$chrlist{$tmp[1]}},$tmp[2]);
}
close($LIST);

my (%prom,%body,%prom_name,%body_name);
while(<$TSV>){
    chomp;
    my @tmp = split;
    foreach(@{$chrlist{$tmp[0]}}){
        if (($_ > $tmp[1]) and ($_ < $tmp[2])) {
            if(exists($prom{$cgid{$tmp[0]."_".$_}})){
                $prom{$cgid{$tmp[0]."_".$_}} .= ";$tmp[6]";
                $prom_name{$cgid{$tmp[0]."_".$_}} .= ";$tmp[7]";
            }else{
                $prom{$cgid{$tmp[0]."_".$_}} = $tmp[6];
                $prom_name{$cgid{$tmp[0]."_".$_}} = $tmp[7];
            }
        }elsif(($_ > $tmp[3]) and ($_ < $tmp[4])) {
            if(exists($body{$cgid{$tmp[0]."_".$_}})){
                $body{$cgid{$tmp[0]."_".$_}} .= ";$tmp[6]";
                $body_name{$cgid{$tmp[0]."_".$_}} .= ";$tmp[7]";
            }else{
                $body{$cgid{$tmp[0]."_".$_}} = $tmp[6];
                $body_name{$cgid{$tmp[0]."_".$_}} = $tmp[7];
            }
        }
    }
}
close($TSV);

foreach my $chr (keys(%chrlist)){
    foreach my $pos (@{$chrlist{$chr}}){
        my $cg = $cgid{$chr."_".$pos};
        $prom{$cg} = "" if (!exists($prom{$cg}));
        $prom_name{$cg} = "" if (!exists($prom_name{$cg}));
        $body{$cg} = "" if (!exists($body{$cg}));
        $body_name{$cg} = "" if (!exists($body_name{$cg}));
        print("$cg\t$chr\t$pos\t$prom{$cg}\t$prom_name{$cg}\t$body{$cg}\t$body_name{$cg}\n");
    }
}

__END__