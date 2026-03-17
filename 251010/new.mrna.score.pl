#!/usr/bin/env perl
use strict;
use warnings;
use autodie;
use POSIX;
use File::Basename;

sub SCORE {
    my ( $TR_END_COUNT, $NC_END_COUNT, $TR_TOTAL, $NC_TOTAL ) = @_;
    my %SCORE;
    my %END_COR;

    for my $CURRENT ( keys %$TR_END_COUNT ) {
        my ( $CHR, $DIR, $POS ) = split( /\t/, $CURRENT );
        my $FORMAL_POS = ( $DIR eq "+" ) ? $POS - 1 : $POS + 1;

        my $T_END = $TR_END_COUNT->{ $CHR . "\t" . $DIR . "\t" . $FORMAL_POS }
          // 1;
        my $N_END = $NC_END_COUNT->{ $CHR . "\t" . $DIR . "\t" . $FORMAL_POS }
          // 1;

        my $T_END_P1 = $TR_END_COUNT->{$CURRENT};
        my $END_COR =
          POSIX::ceil( $TR_END_COUNT->{$CURRENT} * $NC_TOTAL / $TR_TOTAL );
        my $N_END_P1 = $NC_END_COUNT->{$CURRENT} // 1;

        my $SCORE = $T_END_P1 / $T_END - $N_END_P1 / $N_END;
        $SCORE{$CURRENT}   = $SCORE;
        $END_COR{$CURRENT} = $END_COR;
    }
    return ( \%SCORE, \%END_COR );
}

# --- 处理样本名称 ---
my @sample_names;
for my $i ( 0 .. $#ARGV ) {
    my $dir    = dirname( $ARGV[$i] );
    my $prefix = basename($dir);
    push @sample_names, $prefix;
}

# 存储结构
my ( @end_count, @score,       @end_count_cor, @total );
my ( @site_meta, %all_site_id, %score_all,     %base_info );

# 3. 循环读取所有文件 (0 为 NC, 1..n 为 TR)
for my $i ( 0 .. $#ARGV ) {
    $total[$i] = 0;
    open( my $IN, "<", $ARGV[$i] );
    while ( my $line = <$IN> ) {
        chomp $line;
        next if $line =~ /^\s*$/;
        my @tmp = split /\t/, $line;

        # ID: Chr \t Strand \t Pos
        my $id = "$tmp[0]\t$tmp[2]\t$tmp[1]";

        # 记录该位点在哪些样本中出现过
        $all_site_id{$id}++ if $i > 0;

        # 存储该样本下该位点的特有序列统计信息
        # UniqSeqCount, ModeSeq, ModeLen, ModeFreq, LongSeq, LongLen
        $site_meta[$i]{$id} = join( "\t", @tmp[ 4 .. 9 ] );

        # 存储基础信息 (GeneID, N-1, N, N+1)，仅存一次即可
        if ( !exists $base_info{$id} ) {
            $base_info{$id} =
              join( "\t", $tmp[11], $tmp[12], $tmp[13], $tmp[10] );
        }

        # 存储计数
        $end_count[$i]{$id} = $tmp[3];
        $total[$i] += $tmp[3];
    }
    close($IN);
}

# 对每个处理组计算 Score
for my $sample ( 1 .. $#ARGV ) {
    my ( $score_ref, $end_count_cor_ref ) = SCORE(
        \%{ $end_count[$sample] },
        \%{ $end_count[0] },
        $total[$sample], $total[0]
    );
    $score[$sample]         = $score_ref;
    $end_count_cor[$sample] = $end_count_cor_ref;
}

# --- 输出表头 ---
print "#Chr\tPos\tStrand\tN-1\tN\tN+1\tGeneID";

# 为每个样本动态生成其序列信息的表头
my @meta_cols = qw(UniqCount ModeSeq ModeLen ModeFreq LongSeq LongLen Count);
foreach my $name (@sample_names) {
    my $label = ( $name eq $sample_names[0] ) ? "${name}(NC)" : $name;
    foreach my $col (@meta_cols) {
        print "\t${label}_$col";
    }
    if ( $name ne $sample_names[0] ) {
        print "\t${label}_count_corr\t${label}_score";
    }
}
print "\tTotalFoldChange\tAveFC\tTotalScore\tAveScore\n";

# --- 输出内容 ---
for my $id ( keys %all_site_id ) {

    # 仅输出在所有 Treatment 样本中都出现的位点
    if ( $all_site_id{$id} == $#ARGV ) {
        my ( $chr, $dir, $pos ) = split( /\t/, $id );
        my $nc_val = exists $end_count[0]{$id} ? $end_count[0]{$id} : 0;

        my ( $SoaS, $SoaC ) = ( 0, 0 );
        for my $sample ( 1 .. $#ARGV ) {
            $SoaC += $end_count_cor[$sample]{$id} // 0;
            $SoaS += $score[$sample]{$id}         // 0;
        }

        # 阈值过滤 (保持原脚本逻辑)
        if ( $SoaS >= ( 30 * $#ARGV ) && $SoaC >= ( 3 * $#ARGV * $nc_val ) ) {

            # 1. 基础位点信息
            my $output = "$chr\t$pos\t$dir\t$base_info{$id}";

            # 2. NC 样本的信息 (元数据 + 原始Count)
            my $nc_meta = $site_meta[0]{$id} // "0\tNA\t0\t0\tNA\t0";
            $output .= "\t$nc_meta\t$nc_val";

            # 3. 各个 Treatment 样本的信息
            for my $s ( 1 .. $#ARGV ) {
                my $tr_meta  = $site_meta[$s]{$id}     // "0\tNA\t0\t0\tNA\t0";
                my $tr_count = $end_count[$s]{$id}     // 0;
                my $tr_corr  = $end_count_cor[$s]{$id} // 0;
                my $tr_score = $score[$s]{$id}         // 0;
                $output .= "\t$tr_meta\t$tr_count\t$tr_corr\t$tr_score";
            }

            # 4. 计算总计 FC
            my $nc_denom = ( $nc_val == 0 ) ? 1 : $nc_val;
            my $fc       = $SoaC / $nc_denom;
            my $fca      = $fc / $#ARGV;
            $output .= "\t$fc\t$fca";

            $score_all{$output} = $SoaS;
        }
    }
}

# 按 Total Score 降序排序输出
foreach my $row ( sort { $score_all{$b} <=> $score_all{$a} } keys %score_all ) {
    my $ave_score = $score_all{$row} / $#ARGV;
    print "$row\t$score_all{$row}\t$ave_score\n";
}
