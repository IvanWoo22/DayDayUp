#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV;

# CSV文件路径，请根据实际情况修改
my $csv_file    = $ARGV[0];
my $output_file = $ARGV[1];

# 初始化CSV阅读器和写入器
my $csv     = Text::CSV->new( { binary => 1, auto_diag => 1, eol => "\n" } );
my $out_csv = Text::CSV->new( { binary => 1, auto_diag => 1, eol => "\n" } );

# 打开输入和输出文件句柄
open my $fh,  "<:encoding(utf8)", $csv_file    or die "$csv_file: $!";
open my $ofh, ">:encoding(utf8)", $output_file or die "$output_file: $!";

# 写入头部，假设原始CSV有标题行
my $header = $csv->getline($fh);
$header->[ scalar(@$header) ] = 'Extracted_Content';    # 添加新列标题
$out_csv->print( $ofh, $header );

while ( my $row = $csv->getline($fh) ) {
    my $col2 = $row->[1];
    if ( $col2 =~ /\(([^)]+)(?=\)$)/ ) {
        my $extracted_str = $1;    # 提取末尾括号内的内容
        push @$row, $extracted_str;
        $row->[1] =~ s/\([^()]*\)$//;
    }
    else {
        push @$row, '';            # 如果没有匹配到末尾括号内容，则新列为空
    }

    # 写入处理后的行
    $out_csv->print( $ofh, $row );
}

close $fh;
close $ofh;

print "处理完成。\n";
