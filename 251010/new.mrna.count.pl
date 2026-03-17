#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# 定义参数变量
my $dedup = 0;

# 获取命令行参数
GetOptions( 'dedup' => \$dedup )
  or die "Error in command line arguments\n";

# %data 结构:
# $data{KEY} -> {
#    sum_weight   => 数值,
#    longest_gene => "GeneID",
#    seq_data     => { "SEQUENCE" => { count => N, len => 第二列数值 } },
#    longest_seq  => "SEQUENCE",
#    max_seq_len  => 0
# }
my %data;

while ( my $line = <> ) {
    chomp $line;
    next if $line =~ /^\s*$/;

    my @cols = split /\t/, $line;
    next if scalar(@cols) < 8;

    my $raw_len = $cols[1];    # 第二列长度
    my $seq     = $cols[2];    # 第三列序列
    my $score   = $cols[3];    # 第四列分值
    my $chr     = $cols[4];    # 第五列
    my $pos     = $cols[5];    # 第六列
    my $str     = $cols[6];    # 第七列
    my $gene    = $cols[7];    # 第八列

    # --dedup 过滤
    next if ( $dedup && $score != 1 );

    my $key = join( "\t", $chr, $pos, $str );

    unless ( exists $data{$key} ) {
        $data{$key} = {
            sum_weight   => 0,
            longest_gene => "",
            seq_data     => {},
            longest_seq  => "",
            max_seq_len  => -1
        };
    }

    # 1. 权重加和
    $data{$key}->{sum_weight} += $score;

    # 2. 更新最长 Gene (此处仍使用 length() 判断 Gene ID 字符串长度)
    if ( length($gene) > length( $data{$key}->{longest_gene} ) ) {
        $data{$key}->{longest_gene} = $gene;
    }

    # 3. 统计 Seq 频数及对应的原始长度
    $data{$key}->{seq_data}->{$seq}->{count}++;
    $data{$key}->{seq_data}->{$seq}->{len} = $raw_len;

    # 4. 寻找最长 Seq (使用输入的 raw_len 进行比较)
    if ( $raw_len > $data{$key}->{max_seq_len} ) {
        $data{$key}->{max_seq_len} = $raw_len;
        $data{$key}->{longest_seq} = $seq;
    }
}

# 遍历结果并输出
foreach my $k ( keys %data ) {
    my $info     = $data{$k};
    my $seq_hash = $info->{seq_data};

    # 计算该位置下不同序列的数量 (Uniq Count)
    my $uniq_seq_count = scalar( keys %$seq_hash );

    # 计算最高频 Seq
    my $most_freq_seq = "";
    my $max_freq      = 0;
    foreach my $s ( keys %$seq_hash ) {
        if ( $seq_hash->{$s}->{count} > $max_freq ) {
            $max_freq      = $seq_hash->{$s}->{count};
            $most_freq_seq = $s;
        }
    }
    my $most_freq_len = $seq_hash->{$most_freq_seq}->{len};

    # 处理最长 Seq 及其最后三位
    my $longest_seq_str = $info->{longest_seq};
    my $longest_len     = $info->{max_seq_len};

    # 截取碱基 (注意：这里必须基于字符串实际内容，否则无法提取)
    my $actual_str_len = length($longest_seq_str);
    my $b3 =
      ( $actual_str_len >= 3 ) ? substr( $longest_seq_str, -3, 1 ) : "NA";
    my $b2 =
      ( $actual_str_len >= 2 ) ? substr( $longest_seq_str, -2, 1 ) : "NA";
    my $b1 =
      ( $actual_str_len >= 1 ) ? substr( $longest_seq_str, -1, 1 ) : "NA";

    # 输出列对照：
    # 1-3: Chr, Pos, Strand
    # 4: Sum_Weight
    # 5: Uniq_Seq_Count (新增)
    # 6: Max_Gene
    # 7-9: Mode_Seq, Mode_Len, Mode_Freq
    # 10-11: Longest_Seq, Longest_Len
    # 12-14: Last 3 Bases
    print join( "\t",
        $k,               $info->{sum_weight}, $uniq_seq_count,
        $most_freq_seq,   $most_freq_len,      $max_freq,
        $longest_seq_str, $longest_len,        $info->{longest_gene},
        $b3,              $b2,                 $b1 )
      . "\n";
}
