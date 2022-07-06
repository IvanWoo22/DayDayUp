bsub -n 24 -J "sample" '
parallel --keep-order --xapply -j 4 '\''
ALL_STOP={2}
for sampling in 10000 20000 30000 50000 100000 200000 500000 1000000 2000000 5000000; do
COV=$(echo "scale=6; ${sampling}/${ALL_STOP}" | bc)
picard DownsampleSam \
I=output/{1}/mrna.raw.sam.gz \
O=output/{1}/mrna.${sampling}downsampled.bam \
P=${COV}
done
'\'' ::: NJU61{96..99} ::: 14683470 39470116 20594479 21995831
'

#530034465
#201834710
#148016720
#189594931

for PREFIX in 10000 20000 30000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000; do
  bsub -n 6 -J "${PREFIX}" "bash sample1.sh ${PREFIX}"
done

###
MYVAR=$1
SAMPLE=$2
export MYVAR
export SAMPLE
parallel --env MYVAR --env SAMPLE --keep-order -j 12 '
POINT={}
COV=$(echo "scale=6; ${POINT}/${MYVAR}" | bc)
picard DownsampleSam \
I=output/${SAMPLE}/mrna.raw.sam.gz \
O=output/${SAMPLE}/mrna.{}downsampled.bam \
A=1.0E-6 \
P=${COV}
' ::: 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000
###

#360639129
#
#754324830
#194917423
#183102817
#207883507
#409529701
#275644215
#286050617
#286926093

parallel --keep-order --xapply -j 4 '
pigz -dc output/{}/mrna.raw.sam.gz | grep -v "^@" | wc -l
' ::: NJU61{92..95}

bsub -n 24 -J "Sample" "
parallel --keep-order --xapply -j 12 '
bash sample3.sh {2} {1}
' ::: NJU61{84..95} ::: 754324830 194917423 183102817 207883507 409529701 275644215 286050617 286926093 360639129 139754563 149984624 109624503
"

parallel --keep-order -j 4 '
bsub -n 24 -J "{1}" "bash temp_step3.sh {1} {2} 16"
' ::: NJU61{84..95} ::: 100000000

parallel --keep-order -j 8 '
bsub -n 24 -J "{1}" "bash temp_step3.sh {1} {2} 16"
' ::: NJU61{84..91} ::: 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000

bsub -n 24 -J "almostunique" '
parallel --keep-order -j 12 '\''
time pigz -dc temp/{1}/mrna.{2}dedup.tmp.gz > temp/{1}/mrna.{2}dedup.tmp
time bash NJU_seq/mrna_analysis/almostunique_sample.sh \
  temp/{1}/mrna.{2}dedup.tmp \
  data/{1}/R1.mrna.fq.gz \
  temp/{1} \
  temp/{1}/mrna.{2}almostunique.tmp \
  {2}
rm temp/{1}/mrna.{2}dedup.tmp
'\'' ::: NJU61{84..91} ::: 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000
'

bsub -n 12 -J "almostunique" '
parallel --keep-order -j 6 '\''
time pigz -dc temp/{1}/mrna.{2}dedup.tmp.gz > temp/{1}/mrna.{2}dedup.tmp
time bash NJU_seq/mrna_analysis/almostunique_sample.sh \
  temp/{1}/mrna.{2}dedup.tmp \
  data/{1}/R1.mrna.fq.gz \
  temp/{1} \
  temp/{1}/mrna.{2}almostunique.tmp \
  {2}
rm temp/{1}/mrna.{2}dedup.tmp
'\'' ::: NJU61{84..95} ::: 100000000
'

bsub -n 24 -J "almostunique" '
parallel --keep-order -j 12 '\''
time pigz -dc temp/{1}/mrna.{2}dedup.tmp.gz > temp/{1}/mrna.{2}dedup.tmp
time bash NJU_seq/mrna_analysis/almostunique_sample.sh \
temp/{1}/mrna.{2}dedup.tmp \
data/{1}/R1.mrna.fq.gz \
temp/{1} \
temp/{1}/mrna.{2}almostunique.tmp \
{2}
rm temp/{1}/mrna.{2}dedup.tmp
'\'' ::: NJU61{92..99} ::: 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000
'

bsub -n 24 -J "count" '
parallel --keep-order -j 16 '\''
time perl NJU_seq/mrna_analysis/count.pl \
temp/{1}/mrna.{2}almostunique.tmp \
>temp/{1}/mrna.{2}count.tmp
'\'' ::: NJU61{84..99} ::: 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000
'

bsub -n 8 -J "count" '
parallel --keep-order -j 4 '\''
time perl NJU_seq/mrna_analysis/count.pl \
temp/{1}/mrna.{2}almostunique.tmp \
>temp/{1}/mrna.{2}count.tmp
'\'' ::: NJU61{84..95} ::: 100000000
'

bsub -n 24 -J "count" '
parallel --keep-order -j 16 '\''
time pigz -dcf data/ath.gff3.gz |
awk '\''\'\'''\'' $3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9} '\''\'\'''\'' |
perl NJU_seq/mrna_analysis/merge.pl --refstr "Parent=transcript:" --geneid "AT" --transid "AT" -i temp/{1}/mrna.{2}count.tmp -o output/{1}/{2}mrna.tsv
'\'' ::: NJU61{84..99} ::: 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000
'

bsub -n 4 -J "count" '
parallel --keep-order -j 4 '\''
time pigz -dcf data/ath.gff3.gz |
awk '\''\'\'''\'' $3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9} '\''\'\'''\'' |
perl NJU_seq/mrna_analysis/merge.pl --refstr "Parent=transcript:" --geneid "AT" --transid "AT" -i temp/{1}/mrna.{2}count.tmp -o output/{1}/{2}mrna.tsv
'\'' ::: NJU61{84..95} ::: 100000000
'

parallel --keep-order -j 8 '
echo -e "{1}\n{2}"
bash seq_depth.sh \
temp/{1}/mrna.{2}almostunique.tmp \
output/{1}/{2}mrna.tsv
' ::: NJU61{84..99} ::: 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000

parallel --keep-order -j 8 '
echo -e "{1}\n{2}"
bash seq_depth.sh \
temp/{1}/mrna.{2}almostunique.tmp \
output/{1}/{2}mrna.tsv
' ::: NJU61{84..95} ::: 100000000

parallel -j 3 "
perl score_neo.pl \\
output/NJU6196/mrna.tsv \\
output/{}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_scored.tsv
" ::: NJU61{97..99}

parallel --keep-order --xapply -j 8 "
perl NJU_seq/mrna_analysis/extract_point.pl \
  output/NJU{1}/mrna_basic_scored.tsv \
  output/NJU{2}/mrna_basic_scored.tsv \
  output/NJU{3}/mrna_basic_scored.tsv \
  1000 1 >output/{4}_mrna_scored_1000p.tsv
  " ::: $(seq 6185 4 6199) $(seq 6221 4 6235) ::: $(seq 6186 4 6199) $(seq 6222 4 6235) ::: $(seq 6187 4 6199) $(seq 6223 4 6235) ::: Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF

perl NJU_seq/mrna_analysis/extract_point.pl \
  output/NJU6197/mrna_scored.tsv \
  output/NJU6198/mrna_scored.tsv \
  output/NJU6199/mrna_scored.tsv \
  1000 1 >output/Ath_fruit_RF_mrna_scored_1000p.tsv

parallel --keep-order --xapply -j 4 "
perl NJU_seq/mrna_analysis/score_neo.pl \
output/NJU6188/mrna.tsv \
output/{}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \
>output/{}/mrna_scored_test.tsv
" ::: NJU61{85..87}

parallel --keep-order --xapply -j 4 "
perl NJU_seq/mrna_analysis/score_neo.pl \
output/NJU6188/mrna.tsv \
output/{}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \
>output/{}/mrna_scored_test.tsv
" ::: NJU61{97..99}

perl extract_point.pl \
  output/NJU6185/mrna_scored_test.tsv \
  output/NJU6186/mrna_scored_test.tsv \
  output/NJU6187/mrna_scored_test.tsv \
  1000 1 >output/Ath_sixl_RF_mrna_scored_500p_test.tsv

perl extract_point.pl \
  output/NJU6197/mrna_scored_test.tsv \
  output/NJU6198/mrna_scored_test.tsv \
  output/NJU6199/mrna_scored_test.tsv \
  1000 1 >output/Ath_fruit_RF_mrna_scored_500p_test.tsv

parallel --keep-order -j 8 '
wc -l temp/{1}/mrna.{2}almostunique.tmp | awk '\''{print $1}'\''
' ::: NJU61{84..95} ::: 100000000

bowtie2-build data/SMV.fa index/SMV
