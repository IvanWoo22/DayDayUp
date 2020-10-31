use strict;
use warnings;
use autodie;

my ( @chr, @site, @strand );
open( my $IN1, "<", $ARGV[0] );
while (<$IN1>) {
    chomp;
    my @tmp = split "\t";
    push( @chr,    $tmp[0] );
    push( @site,   $tmp[1] );
    push( @strand, $tmp[2] );
}
close($IN1);

open( my $IN2, "<", $ARGV[1] );
while (<$IN2>) {
    chomp;
    s/\r//g;
    my @tmp = split "\t";
    $tmp[0] =~ s/chr//;
    foreach my $i ( 0 .. $#chr ) {
        if (    ( $chr[$i] eq $tmp[0] )
            and ( $strand[$i] eq $tmp[5] )
            and ( $site[$i] > $tmp[1] )
            and ( $site[$i] < $tmp[2] ) )
        {
            print( join( "\t", @tmp ) . "\t$site[$i]\n" );
        }
    }
}
close($IN2);

__END__
