#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my $raw_pattern = $ARGV[0];
my $mis_count   = $ARGV[1];
my $one_match   = fuzzy_pattern( $raw_pattern, $mis_count );
print "$one_match\n";

sub fuzzy_pattern {
    my ( $original_pattern, $mismatches_allowed ) = @_;
    $mismatches_allowed >= 0
      or die "Number of mismatches must be greater than or equal to zero\n";
    my $new_pattern =
      make_approximate( $original_pattern, $mismatches_allowed );
    return qr/$new_pattern/;
}

sub make_approximate {
    my ( $pattern, $mismatches_allowed ) = @_;
    if    ( $mismatches_allowed == 0 ) { return $pattern }
    elsif ( length($pattern) <= $mismatches_allowed ) {
        $pattern =~ tr/ACTG/./;
        return $pattern;
    }
    else {
        my ( $first, $rest ) = $pattern =~ /^(.)(.*)/;
        my $after_match = make_approximate( $rest, $mismatches_allowed );
        if ( $first =~ /[ACGTY]/ ) {
            my $after_miss = make_approximate( $rest, $mismatches_allowed - 1 );
            return "(?:$first$after_match|.$after_miss)";
        }
        else { return "$first$after_match" }
    }
}

__END__