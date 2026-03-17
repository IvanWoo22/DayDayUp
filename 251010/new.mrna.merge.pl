#!/usr/bin/env perl
use strict;
use warnings;

if ( scalar @ARGV < 2 ) {
    die "Usage: perl $0 <File_A> <File_B> > merged.tsv\n";
}

my ( $file_a, $file_b ) = @ARGV;

my %data;      # 存储所有数据: $data{id} = "full_line_content"
my %source;    # 存储来源标记: 1 (A), 2 (B), 3 (Both)
my @order;     # 记录出现的先后顺序（先 A 后 B）
my $header;

# 1. 处理文件 A
open( my $fh_a, "<", $file_a ) or die "Cannot open $file_a: $!";
while (<$fh_a>) {
    chomp;
    next if /^\s*$/;
    if (/^#/) {
        $header = $_ unless defined $header;    # 记录表头
        next;
    }
    my @cols = split /\t/;
    my $id   = join( "\t", @cols[ 0 .. 2 ] );    # Chr, Pos, Strand

    $data{$id}   = $_;
    $source{$id} = 1;                            # 标记为来自 A
    push @order, $id;
}
close($fh_a);

# 2. 处理文件 B
open( my $fh_b, "<", $file_b ) or die "Cannot open $file_b: $!";
while (<$fh_b>) {
    chomp;
    next if /^\s*$/;
    if (/^#/) {
        $header = $_ unless defined $header;
        next;
    }
    my @cols = split /\t/;
    my $id   = join( "\t", @cols[ 0 .. 2 ] );

    if ( exists $source{$id} ) {

        # 如果 A 已经有了，标记为 Both
        $source{$id} = 3;
    }
    else {
        # 如果 A 没有，记录新数据
        $data{$id}   = $_;
        $source{$id} = 2;    # 标记为来自 B
        push @order, $id;
    }
}
close($fh_b);

# 3. 输出结果
# 打印新的表头
print "$header\tMerge_Status\n";

# 按照出现的顺序打印
foreach my $id (@order) {
    my $status;
    if    ( $source{$id} == 1 ) { $status = "dedup"; }
    elsif ( $source{$id} == 2 ) { $status = "dup"; }
    elsif ( $source{$id} == 3 ) { $status = "both"; }
    print "$data{$id}\t$status\n";
}

__END__
