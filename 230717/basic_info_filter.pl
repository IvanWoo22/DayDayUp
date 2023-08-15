#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use utf8;

my %data;

while (<STDIN>) {
    chomp;
    my @fields = split /\t/;

    my $id = $fields[0];
    for ( my $i = 1 ; $i < scalar(@fields) ; $i++ ) {
        my $value = $fields[$i];
        $data{$id}[$i]{$value}++;
    }
}

foreach my $id ( sort keys %data ) {
    print "$id";
    for ( my $i = 1 ; $i < scalar( @{ $data{$id} } ) ; $i++ ) {
        my $max_value  = "";
        my $max_length = 0;
        my $max_count  = 0;

        foreach my $value ( keys %{ $data{$id}[$i] } ) {
            my $len   = length($value);
            my $count = $data{$id}[$i]{$value};
            if (    ( $len > $max_length )
                and ( $value ne "其他" )
                and ( $value ne "未提供" )
                and ( $value ne "" )
                and ( $value ne "无" )
                and ( $value ne "-" ) )
            {
                $max_count  = $count;
                $max_value  = $value;
                $max_length = $len;
            }
            elsif ( ( $count > $max_count )
                and ( $len == $max_length )
                and ( $value ne "其他" )
                and ( $value ne "未提供" )
                and ( $value ne "" )
                and ( $value ne "无" )
                and ( $value ne "-" ) )
            {
                $max_count  = $count;
                $max_value  = $value;
                $max_length = $len;
            }
            elsif ( $value eq "其他" ) {
                $data{$id}[$i]{$value} = 2;
            }
            elsif (( $value eq "未提供" )
                or ( $value eq "" )
                or ( $value eq "无" )
                or ( $value eq "-" ) )
            {
                $data{$id}[$i]{$value} = 1;
            }
        }
        if ( $max_count == 0 ) {
            foreach my $value ( keys %{ $data{$id}[$i] } ) {
                my $count = $data{$id}[$i]{$value};
                if ( $count > $max_count ) {
                    $max_count = $count;
                    $max_value = $value;
                }
            }
        }
        print "\t$max_value";
    }
    print "\n";
}

__END__
