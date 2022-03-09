#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

sub SEQ_TR_TU {
    my $SEQ = shift;
    return ( $SEQ =~ tr/Uu/Tt/r );
}

open my $FA, "<", $ARGV[0];
my $id;
my %sequences;
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

open my $TPM, "<", $ARGV[1];
my %tpm;
while ( my $line = <$TPM> ) {
    chomp $line;
    my ( $raw_id, $tpm ) = split( "\t", $line );
    my @tmp = split( /\./, $raw_id );
    $tpm{ $tmp[0] } = $tpm;
}
close($TPM);

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

open my $LS, "<", $ARGV[2];
my $site = 0;
while ( my $line = <$LS> ) {
    chomp $line;
    $site++;
    my @tmp = split( "\t", $line );
    $tmp[1] = SEQ_TR_TU( $tmp[1] );
    $tmp[2] = SEQ_TR_TU( $tmp[2] );
    if ( split( /\//, $tmp[0] ) > 1 ) {
        my $opt = 0;
        foreach my $key ( split( /\//, $tmp[0] ) ) {
            $opt++;
            my ( $title, $seq ) = FETCH($key);
            open my $OUT, ">", $ARGV[3] . "/site_" . $site . "_" . $opt . ".fa";
            unless ( exists( $tpm{ $tmp[0] } ) ) {
                $tpm{ $tmp[0] } = "NA";
            }
            print $OUT (
                "$title$tmp[3]$tmp[4]$tmp[5]|$tmp[1]|$tmp[2]|$tpm{$tmp[0]}\n$seq\n"
            );
            close($OUT);
        }
    }
    else {
        my ( $title, $seq ) = FETCH( $tmp[0] );
        open my $OUT, ">", $ARGV[3] . "/site_" . $site . ".fa";
        unless ( exists( $tpm{ $tmp[0] } ) ) {
            $tpm{ $tmp[0] } = "NA";
        }
        print $OUT (
            "$title$tmp[3]$tmp[4]$tmp[5]|$tmp[1]|$tmp[2]|$tpm{$tmp[0]}\n$seq\n"
        );
        close($OUT);
    }
}
close($LS);

__END__
