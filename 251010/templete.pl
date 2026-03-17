#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use AlignDB::IntSpan;

my $geneid  = "";
my $transid = "";

Getopt::Long::GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'refstr=s'  => \my $refstr,
    'geneid:s'  => \$geneid,
    'transid:s' => \$transid,
    'in|i=s'    => \my $in_sam,
    'out|o=s'   => \my $out_fq,
    'stdout'    => \my $stdout,
) or Getopt::Long::HelpMessage(1);

if ( ( $out_fq && $stdout ) || ( !$out_fq && !$stdout ) ) {
    die "ERROR: You must specify exactly one of --out|-o or --stdout\n";
}

my %trans_range;
my %trans_chr;
my %trans_dir;
while (<>) {
    chomp;
    my ( $chr, $start, $end, $dir, $info ) = split /\t/;
    $chr  =~ s/chr//;
    $info =~ /$refstr($transid\w+\.[0-9]+)/;
    if ( exists( $trans_chr{$1} ) ) {
        $trans_range{$1}->AlignDB::IntSpan::add_range( $start, $end );
    }
    else {
        $trans_chr{$1}   = $chr;
        $trans_dir{$1}   = $dir;
        $trans_range{$1} = AlignDB::IntSpan->new();
        $trans_range{$1}->AlignDB::IntSpan::add_range( $start, $end );
    }
}

sub COORDINATE_POS {
    my $INDEX = $_[0];
    my $SITE  = $_[1];
    my $ISLAND;
    if ( $trans_dir{$INDEX} eq "+" ) {
        $ISLAND = $trans_range{$INDEX}->AlignDB::IntSpan::at($SITE);
    }
    else {
        $ISLAND = $trans_range{$INDEX}->AlignDB::IntSpan::at( -$SITE );
    }
    my $ABS_SITE =
      $trans_chr{$INDEX} . "\t" . $ISLAND . "\t" . $trans_dir{$INDEX};
    return ($ABS_SITE);
}

open( my $IN_FH, "<", $in_sam );
my ( %gene_id, @qname_out, @qname_abs_site );
while (<$IN_FH>) {
    chomp;
    my ( $qname, $rname, $site, $cigar, $seq ) = split /\t/;
    $rname =~ /($transid\w+\.[0-9]+)/;
    my $trans_id = $1;
    $rname =~ /($geneid\w+)/;
    my $gene_id  = $1;
    my $abs_site = COORDINATE_POS( $trans_id, $site );

    if ( exists( $gene_id{$abs_site} ) ) {
        unless ( $gene_id{$abs_site} =~ /$gene_id/ ) {
            $gene_id{$abs_site} .= "/" . $gene_id;
        }
    }
    else {
        $gene_id{$abs_site} = $gene_id;
    }
    push( @qname_out,      "$qname\t$rname\t$site\t$cigar\t$seq" );
    push( @qname_abs_site, $abs_site );
}
close($IN_FH);

my $OUT_FH;
if ( defined $stdout ) {
    $OUT_FH = *STDOUT;
}
else {
    open( $OUT_FH, ">", $out_fq );
}
foreach ( 0 .. $#qname_out ) {
    print $OUT_FH (
        "$qname_out[$_]\t$qname_abs_site[$_]\t$gene_id{$qname_abs_site[$_]}\n");
}
close($OUT_FH) unless defined $stdout;
__END__
