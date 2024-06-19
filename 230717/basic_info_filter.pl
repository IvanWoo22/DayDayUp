#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use utf8;
use open ':std', ':encoding(UTF-8)';

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
        my $max_value  = "未提供";
        my $max_length = 0;
        my $max_count  = 0;

        foreach my $value ( keys %{ $data{$id}[$i] } ) {
            my $len   = length($value);
            my $count = $data{$id}[$i]{$value};
            if (   ( $value eq "未提供" )
                or ( $value eq "" )
                or ( $value eq "无" )
                or ( $value eq "-" )
                or ( $value eq "其他" ) )
            {
                $data{$id}[$i]{$value} = 1;
            }
            else {
                $len = $len + 9;
                if ( ( $len > $max_length ) or ( $count > $max_count ) ) {
                    $max_count  = $count;
                    $max_value  = $value;
                    $max_length = $len;
                }

            }
        }
        print "\t$max_value";
        if ( $max_count == 0 ) {
            foreach my $value ( keys %{ $data{$id}[$i] } ) {
                my $count = $data{$id}[$i]{$value};
                if ( $count > $max_count ) {
                    $max_count = $count;
                    $max_value = $value;
                }
            }
        }

    }
    print "\n";
}

__END__
