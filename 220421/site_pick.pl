#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use AlignDB::IntSpan;

open( my $SITES, "<", $ARGV[0] );
my %site_list;
while (<$SITES>) {
    chomp;
    my $line = $_;
    my @tmp  = split( ":", $line );
    if ( exists( $site_list{ $tmp[0] } ) ) {
        $site_list{ $tmp[0] }->AlignDB::IntSpan::add( $tmp[1] );
    }
    else {
        $site_list{ $tmp[0] } = AlignDB::IntSpan->new;
        $site_list{ $tmp[0] }->AlignDB::IntSpan::add( $tmp[1] );
    }
}
close($SITES);

while (<STDIN>) {
    chomp;
    my $line = $_;
    my ( $header, $position, undef ) = split( "\t", $line );
    if (    ( exists( $site_list{$header} ) )
        and ( $site_list{$header}->AlignDB::IntSpan::contains_any($position) ) )
    {
        print("$line\n");
    }
}

__END__
