use strict;
use warnings;
use autodie;

my (
    %gid_chr, %gid_strand, %gid_name, %eid_gid, %eid_se,
    %me, %es, %ri, %ad, %aa,
    %at, %ap, %eid_eind, %gid_eid
);

open(my $GID_INFO_IN, "<", $ARGV[0]);
readline($GID_INFO_IN);
while (<$GID_INFO_IN>) {
    chomp;
    my ($gid, $name, $chr, $strand) = split "\t";
    $gid_chr{$gid} = $chr;
    $gid_name{$gid} = $name;
    $gid_strand{$gid} = $strand;
}
close($GID_INFO_IN);

open(my $GID_EID_IN, "<", $ARGV[1]);
readline($GID_EID_IN);
while (<$GID_EID_IN>) {
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
close($GID_EID_IN);

open(my $AID_INFO_IN, "<", $ARGV[2]);
readline($AID_INFO_IN);
while (<$AID_INFO_IN>) {
    chomp;
    my (undef, $gid, $type, undef, $exon) = split "\t";
    if ($type eq "ME") {
        my ($exon1, $exon2) = split /\|/, $exon;
        my @eind1 = split ":", $exon1;
        my @eind2 = split ":", $exon2;
        foreach my $ind (@eind1) {
            $me{ $gid_eid{ $gid . "\t" . $ind } } = 1
                unless exists($me{ $gid_eid{ $gid . "\t" . $ind } });
        }
        foreach my $ind (@eind2) {
            $me{ $gid_eid{ $gid . "\t" . $ind } } = 1
                unless exists($me{ $gid_eid{ $gid . "\t" . $ind } });
        }
    }
    elsif ($type eq "ES") {
        my @eind = split ":", $exon;
        foreach my $ind (@eind) {
            $es{ $gid_eid{ $gid . "\t" . $ind } } = 1
                unless exists($es{ $gid_eid{ $gid . "\t" . $ind } });
        }
    }
    elsif ($type eq "RI") {
        my @eind = split ":", $exon;
        foreach my $ind (@eind) {
            $ri{ $gid_eid{ $gid . "\t" . $ind } } = 1
                unless exists($ri{ $gid_eid{ $gid . "\t" . $ind } });
        }
    }
    elsif ($type eq "AD") {
        my @eind = split ":", $exon;
        foreach my $ind (@eind) {
            $ad{ $gid_eid{ $gid . "\t" . $ind } } = 1
                unless exists($ad{ $gid_eid{ $gid . "\t" . $ind } });
        }
    }
    elsif ($type eq "AA") {
        my @eind = split ":", $exon;
        foreach my $ind (@eind) {
            $aa{ $gid_eid{ $gid . "\t" . $ind } } = 1
                unless exists($aa{ $gid_eid{ $gid . "\t" . $ind } });
        }
    }
    elsif ($type eq "AT") {
        my @eind = split ":", $exon;
        foreach my $ind (@eind) {
            $at{ $gid_eid{ $gid . "\t" . $ind } } = 1
                unless exists($at{ $gid_eid{ $gid . "\t" . $ind } });
        }
    }
    elsif ($type eq "AP") {
        my @eind = split ":", $exon;
        foreach my $ind (@eind) {
            $ap{ $gid_eid{ $gid . "\t" . $ind } } = 1
                unless exists($ap{ $gid_eid{ $gid . "\t" . $ind } });
        }
    }
}
close($AID_INFO_IN);

foreach my $exon_id (sort {$a <=> $b} keys(%eid_gid)) {
    print("chr"
        . $gid_chr{ $eid_gid{$exon_id} } . "\t"
        . $eid_se{$exon_id} . "\t"
        . $exon_id . "\t"
        . $eid_gid{$exon_id} . "\t"
        . $gid_strand{ $eid_gid{$exon_id} } . "\t"
        . $gid_name{ $eid_gid{$exon_id} } . "\t"
        . $eid_eind{$exon_id});
    if (exists($es{$exon_id})) {
        print("\tES");
    }
    else {
        print("\tNULL");
    }
    if (exists($ad{$exon_id})) {
        print("\tAD");
    }
    else {
        print("\tNULL");
    }
    if (exists($aa{$exon_id})) {
        print("\tAA");
    }
    else {
        print("\tNULL");
    }
    if (exists($ri{$exon_id})) {
        print("\tRI");
    }
    else {
        print("\tNULL");
    }
    if (exists($me{$exon_id})) {
        print("\tME");
    }
    else {
        print("\tNULL");
    }
    if (exists($at{$exon_id})) {
        print("\tAT");
    }
    else {
        print("\tNULL");
    }
    if (exists($ap{$exon_id})) {
        print("\tAP");
    }
    else {
        print("\tNULL");
    }
    print("\n");
}

__END__
