#!/usr/bin/perl
use strict;
use warnings;
use utf8;
use open ':std', ':encoding(UTF-8)';

my $input_file  = $ARGV[0];
my $output_file = $ARGV[1];

open( my $in,  '<:encoding(UTF-8)', $input_file )  or die "无法打开输入文件: $!";
open( my $out, '>:encoding(UTF-8)', $output_file ) or die "无法打开输出文件: $!";

while (<$in>) {
    chomp;
    s/\s*$//;
    if (
/^(\d+),([^,]+),.*,(20\d{2}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d+),(20\d{2}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d+|NULL)$/
      )
    {
        my $first_col  = $1;
        my $second_col = $2;
        my $last_col1  = $3;
        my $last_col2  = $4;

        my $middle_content = $_;
        $middle_content =~ s/^\Q$first_col\E,\Q$second_col\E,//;
        $middle_content =~ s/,\Q$last_col1\E,\Q$last_col2\E\s*$//;
        $middle_content =~ s/\s*$//;

        my $comma_count = () = $middle_content =~ /,/g;
        if ( $comma_count == 1 ) {
            my ( $third_col, $fourth_col ) = split /,/, $middle_content, 2;
            print $out
"$first_col,$second_col,$third_col,$fourth_col,$last_col1,$last_col2\n";
        }
        elsif ( $comma_count > 1 ) {
            if ( $middle_content =~ /(.+)\. ,(.+)/ ) {
                my $third_col  = $1;
                my $fourth_col = $2;
                $third_col  =~ s/\s*,/，/g;
                $fourth_col =~ s/\s*,/，/g;
                print $out
"$first_col,$second_col,$third_col,$fourth_col,$last_col1,$last_col2\n";
            }
            else {
                if ( $middle_content =~ /(.*),(.*)/ ) {
                    my $third_col  = $1;
                    my $fourth_col = $2;
                    $third_col  =~ s/\s*,/，/g;
                    $fourth_col =~ s/\s*,/，/g;
                    print $out
"$first_col,$second_col,$third_col,$fourth_col,$last_col1,$last_col2\n";
                }
                else {
                    print STDERR "无法解析的行: $_\n";
                }
            }
        }
    }
    else {
        print "未匹配行: $_\n";
    }
}

close($in);
close($out);
print "文件处理完毕，输出已写入 $output_file。\n";

__END__
