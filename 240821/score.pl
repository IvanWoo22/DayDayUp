#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Text::CSV;

my $csv = Text::CSV->new( { binary => 1 } );

my $beta_file = $ARGV[0];
open( my $BETA, "<", $beta_file ) or die "Could not open '$beta_file': $!";

my $header = $csv->getline($BETA);
shift @$header;
my @sample = @$header;

my %beta;
while ( my $row = $csv->getline($BETA) ) {
    my $row_name = shift @$row;
    for my $i ( 0 .. $#sample ) {
        $beta{$row_name}{ $sample[$i] } = $row->[$i];
    }
}
close $BETA;

my $score_file = $ARGV[1];
my %score;
for my $i ( 0 .. $#sample ) {
    $score{ $sample[$i] } = [];
}
open( my $SCORE, "<", $score_file ) or die "Could not open '$score_file': $!";
while (<$SCORE>) {
    chomp;
    my @fields = split( /\t/, $_ );
    my @marker = @fields[ 0 .. 4 ];
    my @coef   = @fields[ 5 .. 9 ];
    for my $i ( 0 .. $#sample ) {
        if (    ( $beta{ $marker[0] }{ $sample[$i] } ne "NA" )
            and ( $beta{ $marker[1] }{ $sample[$i] } ne "NA" )
            and ( $beta{ $marker[2] }{ $sample[$i] } ne "NA" )
            and ( $beta{ $marker[3] }{ $sample[$i] } ne "NA" )
            and ( $beta{ $marker[4] }{ $sample[$i] } ne "NA" ) )
        {
            my $score =
              $beta{ $marker[0] }{ $sample[$i] } * $coef[0] +
              $beta{ $marker[1] }{ $sample[$i] } * $coef[1] +
              $beta{ $marker[2] }{ $sample[$i] } * $coef[2] +
              $beta{ $marker[3] }{ $sample[$i] } * $coef[3] +
              $beta{ $marker[4] }{ $sample[$i] } * $coef[4];
            push( @{ $score{ $sample[$i] } }, $score );
        }
        else {
            push( @{ $score{ $sample[$i] } }, "NA" );
        }
    }
}
close $SCORE;

for my $i ( 0 .. $#sample ) {
    my $out = join( "\t", @{ $score{ $sample[$i] } } );
    print "$sample[$i]\t$out\n";
}

__END__
