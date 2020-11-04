use strict;
use warnings;
use autodie;
use POSIX;
use AlignDB::IntSpan;

sub DISTANCE {
    my ( $SITE, $SET_STRING ) = @_;
    my $SET = AlignDB::IntSpan->new;
    $SET->add($SET_STRING);
    my ( $ISLAND, $DIS );
    if ( $SET->find_islands($SITE) eq "-" ) {
        $ISLAND = $SET->nearest_island($SITE)->as_string;
        if ( $ISLAND =~ /,/ ) {
            ( $ISLAND, undef ) = split( ",", $ISLAND );
        }
        if ( $ISLAND =~ /\-/ ) {
            my ( $START, $END ) = split( "-", $ISLAND );
            if ( $START > $SITE ) {
                $DIS = $SITE - $START;
            }
            else {
                $DIS = $SITE - $END;
            }
        }
        else {
            $DIS = $SITE - $ISLAND;
        }
    }
    else {
        $DIS = 0;
    }
    return $DIS;
}

my ( %gene, %exon, %exon_site );
open( my $IN1, "<", $ARGV[0] );
while (<$IN1>) {
    chomp;
    my ( $trans, $gene, $no, $exon, $exon_site ) = split "\t";
    $gene{$gene}[ $no - 1 ] = $trans;
    $exon{$trans} = AlignDB::IntSpan->new;
    $exon{$trans}->AlignDB::IntSpan::add($exon);
    $exon_site{$trans} = AlignDB::IntSpan->new;
    $exon_site{$trans}->AlignDB::IntSpan::add($exon_site);
}

open( my $IN2, "<", $ARGV[1] );
while (<$IN2>) {
    chomp;
    my ( undef, $site, $strand, undef, undef, undef, $gene, undef ) =
      split "\t";
    foreach my $no ( 0 .. $#{ $gene{$gene} } ) {
        if ( $exon{ $gene{$gene}[$no] }->AlignDB::IntSpan::contains_all($site) )
        {
            my $distance =
              DISTANCE( $site,
                $exon_site{ $gene{$gene}[$no] }->AlignDB::IntSpan::as_string );
            my $size =
              $exon{ $gene{$gene}[$no] }->AlignDB::IntSpan::find_islands($site)
              ->AlignDB::IntSpan::size;
            if ( $strand eq "-" ) {
                $distance = -$distance;
            }
            print
"$site\t$gene\t$gene{$gene}[$no]\t$exon{$gene{$gene}[$no]}\t$strand\t$exon_site{$gene{$gene}[$no]}\t$distance\t$size";
            foreach ( 1 .. 100 ) {
                my $randis = int( rand( $size + 1 ) );
                if ( ( $size - $randis ) < $randis ) {
                    $randis = $randis - $size;
                }
                if ( $strand eq "-" ) {
                    $randis = -$randis;
                }
                print "\t$randis";
            }
            print "\n";
            last;
        }
    }
}

__END__
