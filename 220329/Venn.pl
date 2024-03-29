#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Algorithm::Combinatorics qw(combinations);

my ( %point, %point_info, @all_set );
foreach my $file_name (@ARGV) {
    open( my $FH, "<", $file_name );
    while (<$FH>) {
        chomp;
        my ( $id, $chr, $pos ) = split( "\t", $_ );
        push( @{ $point{$file_name} }, $id );
        $point_info{$id} = $chr . "\t" . $pos;
    }
    push( @all_set, @{ $point{$file_name} } );
}

my %all_point;
foreach my $ELEMENT (@all_set) {
    $all_point{$ELEMENT}++;
}
my @all_point = keys %all_point;

foreach my $number ( 1 .. $#ARGV + 1 ) {
    my $iter = combinations( \@ARGV, $number );
    while ( my $choose = $iter->next ) {
        my ( @in, @out );
        my %choose;
        foreach my $element ( @{$choose} ) {
            $choose{$element} = 1;
        }
        foreach my $element (@ARGV) {
            push( @in,  @{ $point{$element} } ) if exists( $choose{$element} );
            push( @out, @{ $point{$element} } )
              unless exists( $choose{$element} );
        }
        my ( @in_point, @out_point, %in_point, %out_point );
        foreach my $ELEMENT (@in) {
            $in_point{$ELEMENT}++;
        }
        @in_point = grep { $in_point{$_} == $number; } ( keys %in_point );
        foreach my $ELEMENT (@out) {
            $out_point{$ELEMENT}++;
        }
        @out_point = keys %out_point;
        my @point =
          grep {
                  ( exists( $in_point{$_} ) )
              and ( $in_point{$_} == $number )
              and ( not( exists( $out_point{$_} ) ) )
          } @all_point;
        if ( @point > 0 ) {
            my %no;
            my @list;
            foreach my $i ( 0 .. $#ARGV ) {
                $no{ $ARGV[$i] } = $i;
                $list[$i] = 0;
            }
            foreach my $j (@$choose) {
                $list[ $no{$j} ] = 1;
            }
            foreach my $id (@point) {
                print "$id\t$point_info{$id}\t" . join( "\t", @list ) . "\n";
            }
        }
    }
}

__END__
