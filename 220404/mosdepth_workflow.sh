for i in NJU92{21,26,33,38}; do
  mosdepth --threads 4 --by 1000 ${i}_dep ${i}_bwa_deduplicated.bam
  pigz -dc ${i}_dep.regions.bed.gz |
    perl -ne '@a=split;if($a[0]=~m/chr[0-9]*$/){print join("\t",@a);print"\n"}' \
      >${i}_dep.regions.bed
  perl filterN.pl hg38_N.bed ${i}_dep.regions.bed |
    awk '$4<1' \
      >${i}_filter.bed
  perl rangemerge.pl \
    <${i}_filter.bed |
    awk '$4>4000&&$5<0.2' \
      >${i}_filter_N.tsv
  perl filterbed.pl \
    centromeres.tsv ${i}_filter_N.tsv \
    >${i}_gap.tsv
  perl -ne '@a=split;$a[0]=~s/chr//;print"$a[0]\t$a[1]\t$a[2]\t$a[4]\n"' \
    ${i}_gap.tsv \
    >${i}_gap_input.tsv
done
