use strict;
use warnings;
use autodie;

my ( %trans_id, %info, %trans_count );
open( my $IN1, "<", $ARGV[0] );
while (<$IN1>) {
    chomp;
    my ( $trans, $gene, $chr, $strand, $start, $end ) = split "\t";
    $trans_id{$trans} = 1;
    $info{"$gene\t$chr\t$strand\t$start\t$end"} = $trans;
}
close($IN1);

open( my $IN2, "<", $ARGV[1] );
while (<$IN2>) {
    chomp;
    my ( $trans, $gene, $chr, $strand, $start, $end ) = split "\t";
    if ( exists( $trans_id{$trans} ) ) {
        print "$_\t1\n";
    }
    else {
        if ( not exists( $info{"$gene\t$chr\t$strand\t$start\t$end"} ) ) {
            if ( exists( $trans_count{$gene} ) ) {
                $trans_count{$gene}++;
            }
            else {
                $trans_count{$gene} = 2;
            }
            print "$_\t$trans_count{$gene}\n";
            $info{"$gene\t$chr\t$strand\t$start\t$end"} = $trans;
        }
    }
}
close($IN2);

__END__
