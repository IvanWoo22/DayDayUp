use strict;
use warnings;

open IN, "<", $ARGV[0];

while (<IN>) {
    $_ =~ s/[\r\n]//g;
    my @tmp = split /\t/, $_;
    my $jug;
    my $length;

    my @lsc_l_s;
    my @lsc_l_e;
    my @lsc_r_s;
    my @lsc_r_e;

    my $lsc_5_s;
    my $lsc_5_e;

    my @lsc;

    print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7]\t$tmp[8]\t$tmp[9]\t$tmp[10]\t$tmp[11]\t";
    if ($tmp[6] =~ /,/) {
        my @intspan = split /,/, $tmp[6];
        my @first = split /-/, $intspan[0];
        if ($intspan[1] =~ /-/) {
            my @second = split /-/, $intspan[1];
            $length = $first[1] - $first[0] + $second[1] - $second[0] + 2;
            my $t = $length % 10;
            my $lth = ($length - $t) / 10;
            $jug = $second[1];
            my @dat = ($second[0] .. $second[1], $first[0] .. $first[1]);
            foreach (0 .. 3) {
                $lsc_l_s[$_] = $dat[$_ * $lth];
                $lsc_l_e[$_] = $dat[($_ + 1) * $lth - 1];
                $lsc_r_e[$_] = $dat[$length - 1 - $_ * $lth];
                $lsc_r_s[$_] = $dat[$length - (1 + $_) * $lth];
                if ($lsc_l_s[$_] > $lsc_l_e[$_]) {
                    $lsc[$_] = $lsc_l_s[$_] . "-" . $jug . ",1-" . $lsc_l_e[$_] . "," . $lsc_r_s[$_] . "-" . $lsc_r_e[$_];
                }
                else {
                    $lsc[$_] = $lsc_l_s[$_] . "-" . $lsc_l_e[$_] . "," . $lsc_r_s[$_] . "-" . $lsc_r_e[$_];
                }
            }
            $lsc_5_s = $dat[4 * $lth];
            $lsc_5_e = $dat[$length - 4 * $lth - 1];
            $lsc[4] = $lsc_5_s . "-" . $lsc_5_e;

        }
        else {
            $length = $first[1] - $first[0] + 2;
            my $t = $length % 10;
            my $lth = ($length - $t) / 10;

            my @dat = ($intspan[1], $first[0] .. $first[1]);
            foreach (0 .. 3) {
                $lsc_l_s[$_] = $dat[$_ * $lth];
                $lsc_l_e[$_] = $dat[($_ + 1) * $lth - 1];
                $lsc_r_e[$_] = $dat[$length - 1 - $_ * $lth];
                $lsc_r_s[$_] = $dat[$length - (1 + $_) * $lth];
                if ($_ == 0) {
                    $lsc[$_] = $lsc_l_s[$_] . ",1-" . $lsc_l_e[$_] . "," . $lsc_r_s[$_] . "-" . $lsc_r_e[$_];
                }
                else {
                    $lsc[$_] = $lsc_l_s[$_] . "-" . $lsc_l_e[$_] . "," . $lsc_r_s[$_] . "-" . $lsc_r_e[$_];
                }
            }
            $lsc_5_s = $dat[4 * $lth];
            $lsc_5_e = $dat[$length - 4 * $lth - 1];
            $lsc[4] = $lsc_5_s . "-" . $lsc_5_e;

        }
    }
    else {
        my @first = split /-/, $tmp[6];
        $length = $first[1] - $first[0];
        my $t = $length % 10;
        my $lth = ($length - $t) / 10;
        my @dat = ($first[0] .. $first[1]);
        foreach (0 .. 3) {
            $lsc_l_s[$_] = $dat[$_ * $lth];
            $lsc_l_e[$_] = $dat[($_ + 1) * $lth - 1];
            $lsc_r_e[$_] = $dat[$length - 1 - $_ * $lth];
            $lsc_r_s[$_] = $dat[$length - (1 + $_) * $lth];
            $lsc[$_] = $lsc_l_s[$_] . "-" . $lsc_l_e[$_] . "," . $lsc_r_s[$_] . "-" . $lsc_r_e[$_];
        }
        $lsc_5_s = $dat[4 * $lth];
        $lsc_5_e = $dat[$length - 4 * $lth - 1];
        $lsc[4] = $lsc_5_s . "-" . $lsc_5_e;
    }

    foreach (0 .. 3) {
        print "$lsc[$_]\t";
    }
    print "$lsc[4]\n";
}

close IN;
