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
    my @marker = split( /\+/, $fields[0] );
    my @coef   = split( /,/,  $fields[1] );

    for my $i ( 0 .. $#sample ) {
        my $all_valid = 1;
        for my $m (@marker) {
            if ( !defined $beta{$m}{ $sample[$i] }
                || $beta{$m}{ $sample[$i] } eq "NA" )
            {
                $all_valid = 0;
                last;
            }
        }
        if ( $all_valid != 0 ) {
            my $score = 0;
            for my $j ( 0 .. $#marker ) {
                $score += $beta{ $marker[$j] }{ $sample[$i] } * $coef[$j];
            }
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
