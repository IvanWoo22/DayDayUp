for PREFIX in NBG{0013..0028}; do
  mkdir ${PREFIX}
  rm ../log/${PREFIX}_cutadapt.log
  bsub -n 24 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
-O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1.fq.gz -p ${PREFIX}/R2.fq.gz \
NBG/${PREFIX}_*R1*.gz NBG/${PREFIX}_*R2*.gz -j 20
"
done

for PREFIX in NBG0016; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  bsub -n 3 -o log/${PREFIX}_peco.log -J "${PREFIX}" "
perl NJU_seq/quality_control/pe_consistency.pl \
data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
temp/${PREFIX}.fq.gz
"
done

for PREFIX in NBG{{0013..0015},{0017..0028}}; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  bsub -n 3 -o log/${PREFIX}_peco.log -J "${PREFIX}" "
perl NJU_seq/quality_control/pe_consistency.pl \
data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
temp/${PREFIX}.fq.gz
"
done

time perl NJU_seq/quality_control/fastq_qc.pl \
  temp/NBG{0013..0028}.fq.gz \
  output \
  NBG_group_4to7

for PREFIX in NBG{0013..0028}; do
  bsub -o log/${PREFIX}_rrna_alignment.log -n 24 -J "${PREFIX}" "bash NJU_seq/log/Gmax_scripts/step1.sh ${PREFIX}"
done

for PREFIX in NBG{0013..0028}; do
  rm log/${PREFIX}_rrna_filter.log
  bsub -n 2 -o log/${PREFIX}_rrna_filter.log -J "${PREFIX}" "
samtools view -h -f 97 -F 144 output/${PREFIX}/rrna.raw.sam.gz > output/${PREFIX}/rrna.temp.sam
samtools view -f 145 -F 96 output/${PREFIX}/rrna.raw.sam.gz >> output/${PREFIX}/rrna.temp.sam
samtools sort -n output/${PREFIX}/rrna.temp.sam | samtools view >output/${PREFIX}/rrna.filter.sam
rm output/${PREFIX}/rrna.temp.sam
"
done

for PREFIX in NBG{0013..0028}; do
  bsub -n 4 -o log/${PREFIX}_rrna_judge.log -J "${PREFIX}" "
time awk '\$6!=\"*\"&&\$7==\"=\"{print \$1 \"\t\" \$3 \"\t\" \$4 \"\t\" \$6 \"\t\" \$10}' output/${PREFIX}/rrna.filter.sam |
perl NJU_seq/rrna_analysis/matchquality_judge.pl |
perl NJU_seq/rrna_analysis/multimatch_judge.pl \
>temp/${PREFIX}/rrna.out.tmp
"
done

bsub -n 24 -J "rrna2mrna" "
parallel --xapply --keep-order -j 10 '
time bash NJU_seq/tool/extract_fastq.sh \
temp/{}/rrna.out.tmp \
data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz \
data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
' ::: NBG{0013..0028}
"

for PREFIX in NBG{0013..0028}; do
  bsub -q fat_768 -n 80 -J "${PREFIX}" "bash NJU_seq/log/Gmax_scripts/step2.sh ${PREFIX}"
done

for PREFIX in NBG{0013..0028}; do
  rm log/${PREFIX}_mrna_dedup.log
  bsub -n 24 -o log/${PREFIX}_mrna_dedup.log -J "${PREFIX}" "bash NJU_seq/log/Gmax_scripts/step3.sh ${PREFIX} 16"
done

bsub -n 24 -J "almostunique" '
parallel --keep-order --xapply -j 8 '\''
time pigz -dc temp/{}/mrna_basic.dedup.tmp.gz > temp/{}/mrna_basic.dedup.tmp
time bash NJU_seq/mrna_analysis/almostunique.sh \
temp/{}/mrna_basic.dedup.tmp \
data/{}/R1.mrna.fq.gz \
temp/{} \
temp/{}/mrna_basic.almostunique.tmp
rm temp/{}/mrna_basic.dedup.tmp
'\'' ::: NBG{0013..0028}
'

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 12 '\''
time perl NJU_seq/mrna_analysis/count.pl \
temp/{}/mrna_basic.almostunique.tmp \
>temp/{}/mrna.count.tmp
'\'' ::: NBG{0013..0028}
'

for PREFIX in NBG{0013..0028}; do
  time gzip -dcf data/gmax.gff3.gz |
    awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
    perl NJU_seq/log/Gmax_scripts/merge.pl \
      --refstr "Parent=" \
      --geneid "Glyma." \
      --transid "Glyma." \
      -i temp/"${PREFIX}"/mrna.count.tmp \
      -o output/"${PREFIX}"/mrna.tsv
done

parallel --xapply --keep-order -j 4 '
echo -e "{1}\t{2}\t{3}"
pigz -dc /Users/ivanwoo/Downloads/ENCFF804NTQ.bed.gz |
awk -va={1} -vb={2} -vc={3} '\''$1==a&&$2>b&&$3<c'\''
pigz -dc /Users/ivanwoo/Downloads/ENCFF696OLO.bed.gz |
awk -va={1} -vb={2} -vc={3} '\''$1==a&&$2>b&&$3<c'\''
echo
' ::: chr7 chr20 chr5 chr10 chr11 chr13 chr13 chr1 chr1 chr6 chr18 chr8 chr12 chr5 chr4 chr11 chr12 chr10 chr3 chr12 ::: 55038057 41138550 79111386 124641976 128716170 113034808 40629951 203307167 234230213 132538204 22417745 29061004 65955523 79112056 7697584 128480440 65876922 1517976 165196908 53189917 ::: 55038658 41139151 79111987 124642577 128716771 113035409 40630552 203307768 234230814 132538805 22418346 29061605 65956124 79112657 7698185 128481041 65877523 1518577 165197509 53190518 >>point.list.tsv
