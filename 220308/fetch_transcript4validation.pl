#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

my %sequences;

sub SEQ_TR_TU {
    my $SEQ = shift;
    return ( $SEQ =~ tr/Uu/Tt/r );
}

open my $FA, "<", $ARGV[0];
my $id;
while ( my $line = <$FA> ) {
    chomp $line;
    if ( $line =~ /^(>.*)$/ ) {
        $id = $1;
    }
    elsif ( $line !~ /^\s*$/ ) {
        $sequences{$id} .= $line;
    }
}
close($FA);

sub FETCH {
    my $KEY = shift;
    my $FETCH_TITLE;
    foreach my $TITLE ( keys(%sequences) ) {
        if ( $TITLE =~ /$KEY/ ) {
            $FETCH_TITLE = $TITLE;
            last;
        }
    }
    return ( $FETCH_TITLE, $sequences{$FETCH_TITLE} );
}

open my $LS, "<", $ARGV[1];
my $site = 0;
while ( my $line = <$LS> ) {
    chomp $line;
    $site++;
    my @tmp = split( "\t", $line );
    $tmp[1] = SEQ_TR_TU( $tmp[1] );
    if ( split( /\//, $tmp[0] ) > 1 ) {
        my $opt = 0;
        foreach my $key ( split( /\//, $tmp[0] ) ) {
            $opt++;
            my ( $title, $seq ) = FETCH($key);
            open my $OUT, ">", $ARGV[2] . "/site_" . $site . "_" . $opt . ".fa";
            print $OUT ("$title$tmp[1]|$tmp[2]$tmp[3]$tmp[4]\n$seq\n");
            close($OUT);
        }
    }
    else {
        my ( $title, $seq ) = FETCH( $tmp[0] );
        open my $OUT, ">", $ARGV[2] . "/site_" . $site . ".fa";
        print $OUT ("$title$tmp[1]|$tmp[2]$tmp[3]$tmp[4]\n$seq\n");
        close($OUT);
    }
}
close($LS);

__END__
