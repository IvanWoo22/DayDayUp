use strict;
use warnings;
use autodie;

my (
    %gid_chr, %gid_strand, %gid_name, %eid_gid, %eid_se,
    %aid_psi, %eid_psi, %eid_chara, %eid_eind, %gid_eid
);

open(my $gid_info_in, "<", $ARGV[1]);
readline($gid_info_in);
while (<$gid_info_in>) {
    chomp;
    my ($gid, $name, $chr, $strand) = split "\t";
    $gid_chr{$gid} = $chr;
    $gid_name{$gid} = $name;
    $gid_strand{$gid} = $strand;
}
close($gid_info_in);

open(my $gid_eid_in, "<", $ARGV[2]);
readline($gid_eid_in);
while (<$gid_eid_in>) {
    chomp;
    my ($gid, $eid, $eid_s, $eid_e, $eind) = split "\t";
    $eid_gid{$eid} = $gid;
    if ($eid_s < $eid_e) {
        $eid_se{$eid} = $eid_s . "\t" . $eid_e;
    }
    else {
        $eid_se{$eid} = $eid_e . "\t" . $eid_s;
    }
    $eid_eind{$eid} = $eind;
    $gid_eid{ $gid . "\t" . $eind } = $eid;
}
close($gid_eid_in);

open(my $aid_psi_in, "<", $ARGV[3]);
readline($aid_psi_in);
while (<$aid_psi_in>) {
    chomp;
    my ($aid, $sample, $psi) = split "\t";
    if ($sample == $ARGV[0]) {
        $aid_psi{$aid} = $psi;
    }
}
close($aid_psi_in);

sub UPDATE_PSI {
    my ($PSI, $TYPE, $NOVEL, $GENE, $EIND) = @_;
    foreach my $IND (@{$EIND}) {
        if (exists($eid_psi{ $gid_eid{ $GENE . "\t" . $IND } })) {
            if (
                (
                    ($eid_psi{ $gid_eid{ $GENE . "\t" . $IND } } < $PSI)
                        and ($PSI != 1)
                )
                    or (($PSI != 1)
                    and ($eid_psi{ $gid_eid{ $GENE . "\t" . $IND } } == 1)
                    and ($PSI != 0))
            ) {
                $eid_psi{ $gid_eid{ $GENE . "\t" . $IND } } = $PSI;
                $eid_chara{ $gid_eid{ $GENE . "\t" . $IND } } =
                    $TYPE . "\t" . $NOVEL;
            }
        }
        else {
            $eid_psi{ $gid_eid{ $GENE . "\t" . $IND } } = $PSI;
            $eid_chara{ $gid_eid{ $GENE . "\t" . $IND } } =
                $TYPE . "\t" . $NOVEL;
        }
    }
}

open(my $aid_info_in, "<", $ARGV[4]);
readline($aid_info_in);
while (<$aid_info_in>) {
    chomp;
    my ($aid, $gid, $type, $novel, $exon) = split "\t";
    if (exists($aid_psi{$aid})) {
        if ($type eq "ME") {
            my ($exon1, $exon2) = split /\|/, $exon;
            my @eind1 = split ":", $exon1;
            my @eind2 = split ":", $exon2;
            my $another_psi = 1 - $aid_psi{$aid};
            if ($eind1[0] > $eind2[0]) {
                UPDATE_PSI($aid_psi{$aid}, $type, $novel, $gid, \@eind2);
                UPDATE_PSI($another_psi, $type, $novel, $gid, \@eind1);
            }
            else {
                UPDATE_PSI($aid_psi{$aid}, $type, $novel, $gid, \@eind1);
                UPDATE_PSI($another_psi, $type, $novel, $gid, \@eind2);
            }
        }
        else {
            my @eind = split ":", $exon;
            UPDATE_PSI($aid_psi{$aid}, $type, $novel, $gid, \@eind);
        }
    }
}
close($aid_info_in);

foreach my $exon_id (keys(%eid_psi)) {
    print("chr"
        . $gid_chr{ $eid_gid{$exon_id} } . "\t"
        . $eid_se{$exon_id} . "\t"
        . $exon_id . "\t"
        . $eid_gid{$exon_id} . "\t"
        . $gid_strand{ $eid_gid{$exon_id} } . "\t"
        . $gid_name{ $eid_gid{$exon_id} } . "\t"
        . $eid_eind{$exon_id} . "\t"
        . $eid_psi{$exon_id} . "\t"
        . $eid_chara{$exon_id}
        . "\n");
}

__END__
