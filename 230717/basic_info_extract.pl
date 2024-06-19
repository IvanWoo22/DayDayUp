#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use utf8;

binmode STDOUT, ":encoding(utf8)";
binmode STDIN,  ":encoding(utf8)";
binmode STDERR, ":encoding(utf8)";

while (<STDIN>) {
    my @col = split( "\t", $_ );
    if ( $col[4] =~ /<(?:基本信息|其他基本信息)>\s*(.*)\s*<\/(?:基本信息|其他基本信息)>/ ) {
        my $basic_info = $1;
        if ( $basic_info =~
/姓名:(.*?)\s*职业:(.*?)\s*性别:(.*?)\s*工作单位:(.*?)\s*年龄:(.*?)\s*住址:(.*?)\s*婚姻:(.*?)\s*供史者\(与患者关系\):(.*?)\s*出生地:(.*?)\s*入院日期:(.*?)\s*民族:(.*?)\s*记录日期:(.*)/
          )
        {
            print
"$col[0]\t$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$11\t$10\t$col[3]\t$col[2]\t$col[1]\n";
        }
        elsif ( $basic_info =~
/姓名:(.*?)\s*职业:(.*?)\s*性别:(.*?)\s*工作单位:(.*?)\s*年龄:(.*?)\s*住址:(.*?)\s*婚姻:(.*?)\s*供史者\(与患者关系\):(.*?)\s*出生地:(.*?)\s*入院时间:(.*?)\s*民族:(.*?)\s*记录日期:(.*)/
          )
        {
            print
"$col[0]\t$1\t$2\t$3\t$4\t$5\t$6\t$7\t$8\t$9\t$11\t$10\t$col[3]\t$col[2]\t$col[1]\n";
        }
        elsif ( $basic_info =~
            /姓名:(.*?)\s*职业:(.*?)\s*性别:(.*?)\s*工作单位:(.*?)\s*年龄:(.*?)\s*发病节气:(.*)/
          )
        {
            print
"$col[0]\t$1\t$2\t$3\t$4\t$5\t\t\t\t\t\t\t$col[3]\t$col[2]\t$col[1]\n";
        }
        elsif ( $basic_info =~
/姓名:(.*?)\s*性别:(.*?)\s*年龄:(.*?)\s*出生地:(.*?)\s*职业:(.*?)\s*婚姻:(.*?)\s*民族:(.*?)\s*工作单位:(.*?)\s*住址:(.*)/
          )
        {
            print
"$col[0]\t$1\t$5\t$2\t$8\t$3\t$9\t$6\t\t$4\t$7\t\t$col[3]\t$col[2]\t$col[1]\n";
        }
        elsif ( $basic_info =~
            /姓名:(.*?)\s*性别:(.*?)\s*年龄:(.*?)\s*出生地:(.*?)\s*职业:(.*?)\s*民族:(.*)/ )
        {
            print
"$col[0]\t$1\t$5\t$2\t\t$3\t\t\t\t$4\t$6\t\t$col[3]\t$col[2]\t$col[1]\n";
        }
        elsif ( $basic_info =~
            /姓名:(.*?)\s*性别:(.*?)\s*年龄:(.*?)\s*婚姻:(.*?)\s*职业:(.*)/ )
        {
            print
"$col[0]\t$1\t$5\t$2\t\t$3\t\t$4\t\t\t\t\t$col[3]\t$col[2]\t$col[1]\n";
        }
        elsif ( $basic_info =~ /姓名:(.*?)\s*性别:(.*?)\s*年龄:(.*?)\s*婚姻:(.*)/ ) {
            print
"$col[0]\t$1\t\t$2\t\t$3\t\t$4\t\t\t\t\t$col[3]\t$col[2]\t$col[1]\n";
        }
        elsif ( $basic_info =~ /姓名:(.*?)\s*性别:(.*?)\s*年龄:(.*?)/ ) {
            print
"$col[0]\t$1\t\t$2\t\t$3\t\t\t\t\t\t\t$col[3]\t$col[2]\t$col[1]\n";
        }
        else {
            warn "No sufficient info in: $col[0]: $basic_info\n";
        }
    }
}

__END__
