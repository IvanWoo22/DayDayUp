#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import pandas as pd
from pathlib import Path


def parse_bowtie2_log(log_path):
    total_reads = None
    exactly_one = None
    overall_rate = None

    with open(log_path, 'r') as f:
        content = f.read()

    m_total = re.search(r'^(\d+) reads; of these:$', content, re.MULTILINE)
    if m_total:
        total_reads = int(m_total.group(1))

    m_exact = re.search(r'^ +(\d+) \(\d+\.\d+%\) aligned concordantly exactly 1 time$', content, re.MULTILINE)
    if m_exact:
        exactly_one = int(m_exact.group(1))

    m_rate = re.search(r'^(\d+\.\d+)% overall alignment rate$', content, re.MULTILINE)
    if m_rate:
        overall_rate = float(m_rate.group(1))

    if total_reads is None or exactly_one is None or overall_rate is None:
        raise ValueError("未能完整解析所需字段")

    return {
        'Total_reads': total_reads,
        'Concordant_exact_1': exactly_one,
        'Overall_alignment_rate_%': overall_rate
    }


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("用法: python summarize_bowtie2_logs.py [-o output.tsv] <log1> <log2> ...")
        print("示例:")
        print("  python summarize_bowtie2_logs.py NJU*/rrna.origin.bowtie2.log")
        print("  python summarize_bowtie2_logs.py -o my_summary.tsv NJU*/rrna.origin.bowtie2.log")
        sys.exit(1)

    output_file = "bowtie2_rrna_alignment_summary.tsv"
    log_files = []

    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg in ['-o', '--output']:
            if i + 1 >= len(sys.argv):
                sys.exit("错误：-o/--output 后必须指定输出文件名")
            output_file = sys.argv[i + 1]
            i += 2
        else:
            log_files.append(arg)
            i += 1

    if not log_files:
        sys.exit("错误：没有指定任何 log 文件！")

    data = {}

    for log_path_str in log_files:
        path = Path(log_path_str)
        if not path.exists():
            print("警告：文件不存在，跳过: {}".format(log_path_str))
            continue

        sample_name = path.parent.name
        if not sample_name or sample_name == '.':
            sample_name = path.stem.replace('.bowtie2', '').replace('.origin', '')

        try:
            result = parse_bowtie2_log(path)
            data[sample_name] = result
            print("成功解析: {} <- {}".format(sample_name, path.name))
        except Exception as e:
            print("解析失败 {}: {}".format(log_path_str, e))

    if not data:
        sys.exit("错误：没有成功解析任何文件！")

    df = pd.DataFrame.from_dict(data, orient='index')
    df = df.sort_index()

    if 'Overall_alignment_rate_%' in df.columns:
        df['Overall_alignment_rate_%'] = df['Overall_alignment_rate_%'].round(2)

    df.to_csv(output_file, sep='\t')

    print("\n" + "=" * 60)
    print("Bowtie2 对齐汇总完成！")
    print("共处理 {} 个样本".format(len(df)))
    print("结果保存至: {}".format(output_file))
    print("\n预览：")
    print(df.to_string())