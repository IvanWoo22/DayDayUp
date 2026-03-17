#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use Getopt::Long;
use AlignDB::IntSpan;

=head1 NAME
dedup.pl -- Deduplication by finding out these two sites of transcripts located on the same position.
=head1 SYNOPSIS
    perl dedup.pl --refstr "Parent=" --transid "ENST" --info data/hsa_exon.info
=cut

my $geneid  = "";
my $transid = "";
Getopt::Long::GetOptions(
    'help|h'    => sub { Getopt::Long::HelpMessage(0) },
    'refstr=s'  => \my $refstr,
    'transid:s' => \$transid,
    'geneid:s'  => \$geneid,
    'info=s'    => \my $info_fh,
    'out|o=s'   => \my $out_fq,
    'stdout'    => \my $stdout,
) or Getopt::Long::HelpMessage(1);

if ( ( $out_fq && $stdout ) || ( !$out_fq && !$stdout ) ) {
    die "ERROR: You must specify exactly one of --out|-o or --stdout\n";
}

sub COUNT_ALIGNMENT_LENGTH {
    my $CIGAR  = shift;
    my $LENGTH = 0;
    while ( $CIGAR =~ /([0-9]+)=/g ) {
        $LENGTH += $1;
    }
    while ( $CIGAR =~ /([0-9]+)X/g ) {
        $LENGTH += $1;
    }
    while ( $CIGAR =~ /([0-9]+)D/g ) {
        $LENGTH += $1;
    }
    return $LENGTH;
}

open( my $INFO_FH, "<", $info_fh );
my %trans_range;
my %trans_chr;
my %trans_dir;
while (<$INFO_FH>) {
    chomp;
    my ( $chr, $start, $end, $dir, $info ) = split( /\t/, $_ );
    if ( $info =~ /\b$refstr([^;]+)/ ) {
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
    else {
        warn "Warning: Cannot extract transcript ID from: $info\n";
        next;
    }
}
close($INFO_FH);

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

my $trans_id_regex;
if ( defined $transid && length $transid ) {
    $trans_id_regex = qr/^($transid\w+\.\d+)/;
}
else {
    $trans_id_regex = qr/^(\w+\.\d+)/;
}

my $gene_id_regex;
if ( defined $geneid && length $geneid ) {
    $gene_id_regex = qr/^($geneid\w+)/;
}
else {
    $gene_id_regex = qr/^(\w+)/;
}

my ( %read_site_exist, %gene_id, @read_out, @read_seq_info, @read_abs_site,
    %read_exist );
while (<>) {
    chomp;
    my ( $read_name, $map_name, $site, $cigar, $seq ) = split /\t/;
    my $length  = COUNT_ALIGNMENT_LENGTH($cigar);
    my $end_pos = $site + $length - 2;
    if ( $map_name =~ $trans_id_regex ) {
        my $trans_id = $1;
        if ( $map_name =~ $gene_id_regex ) {
            my $gene_id  = $1;
            my $abs_site = COORDINATE_POS( $trans_id, $end_pos );
            if ( exists( $gene_id{$abs_site} ) ) {
                unless ( $gene_id{$abs_site} =~ /$gene_id/ ) {
                    $gene_id{$abs_site} .= "/" . $gene_id;
                }
            }
            else {
                $gene_id{$abs_site} = $gene_id;
            }
            my $read_abs_site = join( "\t", $read_name, $abs_site );
            unless ( exists $read_site_exist{$read_abs_site} ) {
                $read_site_exist{$read_abs_site} = 1;
                push( @read_out,      "$read_name" );
                push( @read_seq_info, "$length\t$seq" );
                push( @read_abs_site, $abs_site );
                if ( exists( $read_exist{$read_name} ) ) {
                    $read_exist{$read_name}++;
                }
                else {
                    $read_exist{$read_name} = 1;
                }
            }
        }
        else {
            warn "Warning: Cannot extract gene ID from: $map_name\n";
            next;
        }
    }
    else {
        warn "Warning: Cannot extract transcript ID from: $map_name\n";
        next;
    }
}

my $OUT_FH;
if ( defined $stdout ) {
    $OUT_FH = *STDOUT;
}
else {
    open( $OUT_FH, ">", $out_fq );
}
foreach ( 0 .. $#read_out ) {
    my $count = 1 / $read_exist{ $read_out[$_] };
    print $OUT_FH (
"$read_out[$_]\t$read_seq_info[$_]\t$count\t$read_abs_site[$_]\t$gene_id{$read_abs_site[$_]}\n"
    );
}
close($OUT_FH) unless defined $stdout;

__END__
