use strict;
use warnings;
use autodie;
use AlignDB::IntSpan;

sub DISTANCE {
    my ($SITE, $SET_STRING) = @_;
    my $SET = AlignDB::IntSpan->new;
    $SET->add($SET_STRING);
    my ($ISLAND, $DIS, $DIR);
    if ($SET->find_islands($SITE) eq "-") {
        $ISLAND = $SET->nearest_island($SITE);
        if ($ISLAND =~ /,/) {
            ($ISLAND, undef) = split(",", $ISLAND);
        }
        if ($ISLAND =~ /-/) {
            my ($START, $END) = split("-", $ISLAND);
            if ($START > $SITE) {
                $DIS = "$SITE-$START";
                $DIR = -1;
            }
            else {
                $DIS = "$END-$SITE";
                $DIR = 1;
            }
        }
        else {
            if ($SITE > $ISLAND->as_string) {
                $DIS = "$ISLAND-$SITE";
                $DIR = 1;
            }
            else {
                $DIS = "$SITE-$ISLAND";
                $DIR = -1;
            }
        }
    }
    else {
        $DIS = 0;
        $DIR = 0;
    }
    return ($DIS, $DIR);
}

my (%gene, %exon, %intron, %start_codon);
open(my $IN1, "<", $ARGV[0]);
while (<$IN1>) {
    chomp;
    my ($trans, $gene, $no, undef, $exon, $intron, $start_codon) = split "\t";
    $gene{$gene}[ $no - 1 ] = $trans;
    $exon{$trans} = AlignDB::IntSpan->new;
    $exon{$trans}->AlignDB::IntSpan::add($exon);
    $intron{$trans} = AlignDB::IntSpan->new;
    $intron{$trans}->AlignDB::IntSpan::add($intron);
    $start_codon{$trans} = AlignDB::IntSpan->new;
    $start_codon{$trans}->AlignDB::IntSpan::add_runlist($start_codon);
}
close($IN1);

open(my $IN2, "<", $ARGV[1]);
while (<$IN2>) {
    chomp;
    my (undef, $site, $strand, undef, undef, undef, $gene, undef) =
        split "\t";
    foreach my $no (0 .. $#{$gene{$gene}}) {
        if ($exon{ $gene{$gene}[$no] }->AlignDB::IntSpan::contains_all($site)) {
            my $distance;
            my ($distance_neo, $direction) =
                DISTANCE($site, $start_codon{ $gene{$gene}[$no] });
            if ($distance_neo ne "0") {
                my $set1 = AlignDB::IntSpan->new;
                $set1->add($distance_neo);
                my $set2 = $set1->diff($intron{ $gene{$gene}[$no] });
                if ($strand eq "+") {
                    $distance = ($set2->cardinality - 1) * $direction * 1;
                }
                else {
                    $distance = ($set2->cardinality - 1) * $direction * (-1);
                }
            }
            else {
                $distance = 0;
            }
            print
                "$site\t$gene\t$gene{$gene}[$no]\t$exon{$gene{$gene}[$no]}\t$start_codon{$gene{$gene}[$no]}\t$distance\n";
            last;
        }
    }
}
close($IN2);

__END__
