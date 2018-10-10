#!/usr/bin/perl -w

print "How much money do you take?\n";
my $m = <STDIN>;

print "How many caps can exchange to a bottle?\n";
my $n = <STDIN>;

print "How much is a bottle of Coca?\n";
my $l = <STDIN>;



my $s = int($m / $l);
$m = $s;

while ($m >= $n){
	my $i = $m % $n;
	my $j = int($m / $n);
	$s += $j;
	$m = $i + $j;
}	
	


print "You can drink $s bottles at most.\n";
print "You're so poor.\n" if ($s == 0);
