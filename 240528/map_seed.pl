#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub READ_COMPONENTS {
    my ($FILE) = @_;
    my %COMPONENTS;
    open( my $FH, '<', $FILE ) or die "Could not open file '$FILE' $!";
    while (<$FH>) {
        chomp;
        my ( $SET, $INTERSECTION, $NAME, $SEQ, undef ) = split /\t/;
        my @NAMES  = split( /\|/, $NAME );
        my @SEQS   = split( /\|/, $SEQ );
        my $prefix = substr( $INTERSECTION, 0, 15 );
        push @{ $COMPONENTS{$prefix} },
          {
            set          => $SET,
            intersection => $INTERSECTION,
            names        => \@NAMES,
            sequences    => \@SEQS
          };
    }
    close $FH;
    return \%COMPONENTS;
}

sub READ_INPUT {
    my ($FILE) = @_;
    my @INPUT;
    open( my $FH, '<', $FILE ) or die "Could not open file '$FILE' $!";
    while (<$FH>) {
        chomp;
        my ( $SEQ, $VAL ) = split /\t/;
        push @INPUT, { sequence => $SEQ, value => $VAL } if $VAL >= 5;
    }
    close $FH;
    return \@INPUT;
}

sub FIND_BEST {
    my ( $INPUT, $COMPONENTS ) = @_;
    my $BEST_SET;
    my $BEST_SEQ_NAME;
    my $BEST_SEQ;
    my $REMAINING_SEQ;
    my $HIGHEST_SIMILARITY = 0;
    foreach my $COMPONENT (@$COMPONENTS) {
        if ( index( $INPUT, $COMPONENT->{intersection} ) == 0 ) {
            for ( my $i = 0 ; $i < @{ $COMPONENT->{sequences} } ; $i++ ) {
                if ( index( $INPUT, $COMPONENT->{sequences}[$i] ) == 0 ) {
                    my $SIMILARITY = length( $COMPONENT->{sequences}[$i] );
                    if ( $SIMILARITY > $HIGHEST_SIMILARITY ) {
                        $HIGHEST_SIMILARITY = $SIMILARITY;
                        $BEST_SET           = $COMPONENT->{set};
                        $BEST_SEQ_NAME      = $COMPONENT->{names}[$i];
                        $BEST_SEQ           = $COMPONENT->{sequences}[$i];
                        $REMAINING_SEQ      = substr( $INPUT, $SIMILARITY );
                    }
                }
            }
            last;
        }
    }
    return ( $BEST_SET, $BEST_SEQ_NAME, $BEST_SEQ, $REMAINING_SEQ );
}

my $components_file = $ARGV[0];
my $input_file      = $ARGV[1];
my $components      = READ_COMPONENTS($components_file);
my $input_sequences = READ_INPUT($input_file);

foreach my $input (@$input_sequences) {
    my $prefix = substr( $input->{sequence}, 0, 15 );
    if ( exists $components->{$prefix} ) {
        my ( $best_set, $best_seq_name, $best_seq, $remaining_seq ) =
          FIND_BEST( $input->{sequence}, $components->{$prefix} );
        if ( defined $best_set ) {
            print
"$input->{sequence}\t$input->{value}\t$best_set\t$best_seq_name\t$best_seq\t$remaining_seq\n";
        }
        else {
            print
"$input->{sequence}\t$input->{value}\tNA\tNA\tNA\t$input->{sequence}\n";
        }
    }
    else {
        print
"$input->{sequence}\t$input->{value}\tNA\tNA\tNA\t$input->{sequence}\n";
    }
}

__END__
