use strict;
use warnings;

open IN, "<", $ARGV[0];


while(<IN>){
	$_ =~ s/[\r\n]//g;
	my @tmp = split /\t/, $_;
	my $start;
	my $end;
	my $length;
	print "$tmp[0]\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t$tmp[7]\t$tmp[8]\t$tmp[9]\t$tmp[10]\t$tmp[11]\t";
	if ($tmp[6] =~ /,/){
		my @intspan = split/,/,$tmp[6];
		my @first = split/-/,$intspan[0];
		if ($intspan[1] =~ /-/){
			my @second = split/-/,$intspan[1];
			$length = $first[1] - $first[0] + $second[1] - $second[0] + 2;
			my $t = $length % 4;
			my $lth = ($length - $t) / 4;
			my $lsc_1_l_e = $lth - ($second[1] - $second[0] + 1);
			my $lsc_1_r_s = $first[1] - $lth + 1;
			my $lsc_1 = $intspan[1].",1-".$lsc_1_l_e.",".$lsc_1_r_s."-".$first[1];
			my $lsc_2_s = $lsc_1_l_e + 1;
			my $lsc_2_e = $first[1] - $lth;
			my $lsc_2 = $lsc_2_s."-".$lsc_2_e;			
			print "$lsc_1\t$lsc_2\n";
			
		}else{
                        $length = $first[1] - $first[0] + 2;
			my $t = $length % 4;
			my $lth = ($length - $t) / 4;
			my $lsc_1_l_e = $lth - 1;
			my $lsc_1_r_s = $first[1] - $lth + 1;
			my $lsc_1 = $intspan[1].",1-".$lsc_1_l_e.",".$lsc_1_r_s."-".$first[1];
			my $lsc_2_s = $lsc_1_l_e + 1;
			my $lsc_2_e = $first[1] - $lth;
			my $lsc_2 = $lsc_2_s."-".$lsc_2_e;			
			print "$lsc_1\t$lsc_2\n";
		}
	}else{
		my @first = split/-/,$tmp[6];
		$start = $first[0];
		$end = $first[1];
		$length = $first[1] - $first[0];
		my $t = $length % 4;
		my $lth = ($length - $t) / 4;
		my $lsc_1_l_e = $start + $lth - 1;
		my $lsc_1_r_s = $end - $lth + 1;
		my $lsc_1 = "1-".$lsc_1_l_e.",".$lsc_1_r_s."-".$end;
		my $lsc_2_s = $start + $lth;
		my $lsc_2_e = $end - $lth;
		my $lsc_2 = $lsc_2_s."-".$lsc_2_e;
		print "$lsc_1\t$lsc_2\n";
	}
}

close IN;
