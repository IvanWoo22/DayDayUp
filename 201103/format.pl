use strict;
use warnings;
use autodie;

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

while (<STDIN>) {
    chomp;
    my ( $chr, $start, $end, $strand, $info ) = split /\t/;
    $chr =~ s/chr//;
    my $trans = GET_ID( "transcript_id=", $info );
    my $gene  = GET_ID( "gene_id=",       $info );
    print "$trans\t$gene\t$chr\t$strand\t$start\t$end\n";
}
