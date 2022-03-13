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
    my @FETCH_TITLE;
    my @FETCH_SEQ;
    foreach my $TITLE ( keys(%sequences) ) {
        if ( $TITLE =~ /$KEY/ ) {
            push(@FETCH_TITLE,$TITLE);
            push(@FETCH_SEQ,$sequences{$TITLE});
        }
    }
return ( \@FETCH_TITLE, \@FETCH_SEQ, $#FETCH_TITLE );
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
            my ( $title, $seq, $count ) = FETCH($key);
	    if($count==0){
            open my $OUT, ">", $ARGV[3] . "/site_" . $site . "_" . $opt . ".fa";
            unless ( exists( $tpm{ $tmp[0] } ) ) {
                $tpm{ $tmp[0] } = "N/A";
            }
		print $OUT (
"${$title}[0]$tmp[3]$tmp[4]$tmp[5]|$tmp[1]|$tmp[2]|$tpm{$tmp[0]}\n${$seq}[0]\n"
            	);
            close($OUT);
	    }else {
            foreach my $optt (1 .. ($count + 1)) {
                open my $OUT, ">", $ARGV[3] . "/site_" . $site . "_" . $opt . "_" . $optt . ".fa";
                unless (exists($tpm{ $tmp[0] })) {
                    $tpm{ $tmp[0] } = "N/A";
                }
                print $OUT (
                    "${$title}[$optt - 1]$tmp[3]$tmp[4]$tmp[5]|$tmp[1]|$tmp[2]|$tpm{$tmp[0]}\n${$seq}[$optt - 1]\n"
                );
                close($OUT);
            }
        }
        }
    }
    else {
        my ( $title, $seq, $count ) = FETCH( $tmp[0] );
        if($count==0) {
            open my $OUT, ">", $ARGV[3] . "/site_" . $site . ".fa";
            unless (exists($tpm{ $tmp[0] })) {
                $tpm{ $tmp[0] } = "N/A";
            }
            print $OUT (
                "${$title}[0]$tmp[3]$tmp[4]$tmp[5]|$tmp[1]|$tmp[2]|$tpm{$tmp[0]}\n${$seq}[0]\n"
            );
            close($OUT);
        }else{
            foreach my $optt (1 .. ($count + 1)) {
                open my $OUT, ">", $ARGV[3] . "/site_" . $site . "_0_" . $optt . ".fa";
                unless (exists($tpm{ $tmp[0] })) {
                    $tpm{ $tmp[0] } = "N/A";
                }
                print $OUT (
                    "${$title}[$optt - 1]$tmp[3]$tmp[4]$tmp[5]|$tmp[1]|$tmp[2]|$tpm{$tmp[0]}\n${$seq}[$optt - 1]\n"
                );
                close($OUT);
            }
        }
    }
}
close($LS);

__END__

