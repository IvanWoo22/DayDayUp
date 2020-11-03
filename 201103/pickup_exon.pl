use strict;
use warnings;
use autodie;

use Getopt::Long;
use AlignDB::IntSpan;

sub GET_ID {
    my ( $FEATURE, $INFO ) = @_;
    my $ID;
    if ( $INFO =~ m/$FEATURE(\w+[0-9]*)/ ) {
        $ID = $1;
    }
    else {
        warn("There is a problem in $INFO;\n");
    }
    return $ID;
}

my ( %trans, %gene, %no, %chr_dir, %exon, %exon_site );

open( my $IN_FH1, "<", $ARGV[0] );
while (<$IN_FH1>) {
    chomp;
    my ( $trans, $gene, $chr, $dir, undef, undef, $no ) = split "\t";
    $trans{$trans}   = 1;
    $gene{$trans}    = $gene;
    $chr_dir{$trans} = "$chr\t$dir";
    $no{$trans}      = $no;
}
close($IN_FH1);

while (<STDIN>) {
    chomp;
    my ( undef, undef, $type, $start, $end, undef, undef, undef, $info ) =
      split /\t/;
    my $trans_id = GET_ID( "transcript_id=", $info );
    if ( exists( $trans{$trans_id} ) ) {
        if ( $type eq "exon" ) {
            if ( exists( $exon{$trans_id} ) ) {
                $exon{$trans_id}->AlignDB::IntSpan::add_range( $start, $end );
                $exon_site{$trans_id}->AlignDB::IntSpan::add( $start, $end );
            }
            else {
                $exon{$trans_id} = AlignDB::IntSpan->new;
                $exon{$trans_id}->AlignDB::IntSpan::add_range( $start, $end );
                $exon_site{$trans_id} = AlignDB::IntSpan->new;
                $exon_site{$trans_id}->AlignDB::IntSpan::add( $start, $end );
            }
        }
    }
}

foreach my $trans ( keys(%trans) ) {
    print
      "$trans\t$gene{$trans}\t$no{$trans}\t$exon{$trans}\t$exon_site{$trans}\n";
}
