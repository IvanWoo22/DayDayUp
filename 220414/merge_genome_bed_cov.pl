#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use AlignDB::IntSpan;

my @chr;
foreach my $FILE ( 0 .. ( $#ARGV - 1 ) ) {
    open( my $CHR_RANGE, "<", $ARGV[$FILE] );
    while (<$CHR_RANGE>) {
        chomp;
        my @tmp = split;
        if ( exists( ${ $chr[$FILE] }{ $tmp[0] } ) ) {
            ${ $chr[$FILE] }{ $tmp[0] }
              ->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
        }
        else {
            ${ $chr[$FILE] }{ $tmp[0] } = AlignDB::IntSpan->new;
            ${ $chr[$FILE] }{ $tmp[0] }
              ->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
        }
    }
    close($CHR_RANGE);
}

open( my $WINDOW, "<", $ARGV[$#ARGV] );
while (<$WINDOW>) {
    chomp;
    my @tmp = split;
    my $win = AlignDB::IntSpan->new;
    $win->AlignDB::IntSpan::add_pair( $tmp[1], $tmp[2] );
    my $count = 0;
    foreach my $sample ( 0 .. ( $#ARGV - 1 ) ) {
        if (
            ${ $chr[$sample] }{ $tmp[0] }->AlignDB::IntSpan::overlap($win) > 0 )
        {
            $count++;
        }
    }
    print("$tmp[0]\t$tmp[1]\t$tmp[2]\t$count\n");
}
close($WINDOW);

__END__
