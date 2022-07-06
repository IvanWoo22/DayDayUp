for PREFIX in NJU{6276..6355}; do
  time parallel -j 3 "
perl NJU_seq/rrna_analysis/readend_count.pl \\
NJU_seq/data/ath_rrna/{}.fa temp/${PREFIX}/rrna.out.tmp {} \\
>output/${PREFIX}/rrna_{}.tsv
" ::: 25s 18s 5-8s
done

time parallel -j 8 "
perl NJU_seq/rrna_analysis/score.pl \\
output/{1}/rrna_18s.tsv \\
output/{2}/rrna_18s.tsv \\
output/{3}/rrna_18s.tsv \\
output/{4}/rrna_18s.tsv \\
>output/{1}/rrna_18s_scored.tsv
" ::: NJU{6276..6355..4} ::: NJU{6277..6355..4} ::: NJU{6278..6355..4} ::: NJU{6279..6355..4}

for PREFIX in NJU{6276..6355..4}; do
  bash NJU_seq/presentation/point_venn.sh \
    Sample1 output/${PREFIX}/rrna_25s_scored.tsv 14Sample2 output/${PREFIX}/rrna_25s_scored.tsv 15Sample3 output/${PREFIX}/rrna_25s_scored.tsv 16output/${PREFIX}/rrna_25s_venn.png 50
done

bsub -n 24 -q largemem -J "almostunique" '
parallel --keep-order --xapply -j 6 '\''
time pigz -dc temp/NJU{}/mrna.dedup.tmp.gz > temp/NJU{}/mrna.dedup.tmp
time bash NJU_seq/mrna_analysis/almostunique.sh \
temp/NJU{}/mrna.dedup.tmp \
data/NJU{}/R1.mrna.fq.gz \
temp/NJU{} \
temp/NJU{}/mrna.almostunique.tmp
rm temp/NJU{}/mrna.dedup.tmp
'\'' ::: {6276..6355}
'

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 16 '\''
time perl NJU_seq/mrna_analysis/count.pl \
temp/NJU{}/mrna.almostunique.tmp \
>temp/NJU{}/mrna.count.tmp
time pigz -dcf data/ath.gff3.gz |
awk '\''\'\'''\'' $3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9} '\''\'\'''\'' |
perl NJU_seq/mrna_analysis/merge.pl \
--refstr "Parent=transcript:" \
--geneid "AT" \
--transid "AT" \
-i temp/NJU{}/mrna.count.tmp \
-o output/NJU{}/mrna.tsv
'\'' ::: {6276..6355}
'

parallel --keep-order -j 4 '
echo NJU{} >>output/NJU{}/mrna.cov
bash NJU_seq/presentation/seq_depth.sh \
temp/NJU{}/mrna.almostunique.tmp \
output/NJU{}/mrna.tsv \
>>output/NJU{}/mrna.cov
' ::: {6276..6355}

parallel --keep-order --xapply -j 8 "
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna.tsv \\
output/{2}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{2}/mrna_scored.tsv
perl NJU_seq/mrna_analysis/score_neo.pl \\
output/{1}/mrna.tsv \\
output/{3}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{3}/mrna_scored.tsv
perl NJU_seq/mrna_analysis/score_neo.pl \\
output/{1}/mrna.tsv \\
output/{4}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{4}/mrna_scored.tsv
" ::: NJU{6276..6355..4} ::: NJU{6277..6355..4} ::: NJU{6278..6355..4} ::: NJU{6279..6355..4}

parallel --keep-order --xapply -j 8 "
perl extract_point.pl \
output/{1}/mrna_scored.tsv \
output/{2}/mrna_scored.tsv \
output/{3}/mrna_scored.tsv \
500 1 >output/{4}_mrna_scored_500p.tsv
" ::: NJU{6277..6355..4} ::: NJU{6278..6355..4} ::: NJU{6279..6355..4} ::: Ath_ck_leaf_RF Ath_ck_root_RF Ath_ck_leaf Ath_ck_root Ath_Al_20mM_leaf_RF Ath_Al_20mM_leaf Ath_Al_20mM_root_RF Ath_Al_20mM_root Ath_Al_60mM_leaf_RF Ath_Al_60mM_leaf Ath_Al_60mM_root_RF Ath_Al_60mM_root Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_leaf Ath_Cd_5mM_root_RF Ath_Cd_5mM_root Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_leaf Ath_Cd_15mM_root_RF Ath_Cd_15mM_root

for PREFIX in NJU{6184..6215} NJU{6220..6251} NJU{6276..6355}; do
  if [ -e output/${PREFIX}/mrna_basic.cov ]; then
    mv output/${PREFIX}/mrna_basic.cov output/${PREFIX}/mrna.cov
  fi
done

parallel --keep-order --xapply -j 8 "
perl NJU_seq/mrna_analysis/score_neo.pl \\
output/{1}/mrna_basic.tsv \\
output/{2}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{2}/mrna_scored_neo.tsv
perl NJU_seq/mrna_analysis/score_neo.pl \\
output/{1}/mrna_basic.tsv \\
output/{3}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{3}/mrna_scored_neo.tsv
perl NJU_seq/mrna_analysis/score_neo.pl \\
output/{1}/mrna_basic.tsv \\
output/{4}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{4}/mrna_scored_neo.tsv
" ::: NJU{6092..6099..4} ::: NJU{6093..6099..4} ::: NJU{6094..6099..4} ::: NJU{6095..6099..4}

parallel --keep-order --xapply -j 2 "
perl NJU_seq/mrna_analysis/extract_point.pl \
output/{1}/mrna_scored_neo.tsv \
output/{2}/mrna_scored_neo.tsv \
output/{3}/mrna_scored_neo.tsv \
1000 1 1>/dev/null
" ::: NJU{6093..6099..4} ::: NJU{6094..6099..4} ::: NJU{6095..6099..4}

parallel --keep-order --xapply -j 8 "
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna_basic.tsv \\
output/{2}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{2}/mrna_scored.tsv
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna_basic.tsv \\
output/{3}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{3}/mrna_scored.tsv
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna_basic.tsv \\
output/{4}/mrna_basic.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{4}/mrna_scored.tsv
" ::: NJU{6092..6099..4} ::: NJU{6093..6099..4} ::: NJU{6094..6099..4} ::: NJU{6095..6099..4}

parallel --keep-order --xapply -j 2 "
perl NJU_seq/mrna_analysis/extract_point.pl \
output/{1}/mrna_scored.tsv \
output/{2}/mrna_scored.tsv \
output/{3}/mrna_scored.tsv \
1000 1 1>/dev/null
" ::: NJU{6093..6099..4} ::: NJU{6094..6099..4} ::: NJU{6095..6099..4}

bsub -n 24 -q largemem -J "almostunique" '
parallel --keep-order --xapply -j 6 '\''
time pigz -dc temp/NJU{}/mrna.dedup.tmp.gz > temp/NJU{}/mrna.dedup.tmp
time bash NJU_seq/mrna_analysis/almostunique.sh \
temp/NJU{}/mrna.dedup.tmp \
data/NJU{}/R1.mrna.fq.gz \
temp/NJU{} \
temp/NJU{}/mrna.almostunique.tmp
rm temp/NJU{}/mrna.dedup.tmp
'\'' ::: {6200..6215} {6236..6251}
'

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 16 '\''
time perl NJU_seq/mrna_analysis/count.pl \
temp/NJU{}/mrna.almostunique.tmp \
>temp/NJU{}/mrna.count.tmp
time pigz -dcf data/ath.gff3.gz |
awk '\''\'\'''\'' $3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9} '\''\'\'''\'' |
perl NJU_seq/mrna_analysis/merge.pl \
--refstr "Parent=transcript:" \
--geneid "AT" \
--transid "AT" \
-i temp/NJU{}/mrna.count.tmp \
-o output/NJU{}/mrna.tsv
'\'' ::: {6200..6215} {6236..6251}
'

parallel --keep-order -j 4 '
echo NJU{} >>output/NJU{}/mrna.cov
bash NJU_seq/presentation/seq_depth.sh \
temp/NJU{}/mrna.almostunique.tmp \
output/NJU{}/mrna.tsv \
>>output/NJU{}/mrna.cov
' ::: {6200..6215} {6236..6251}

parallel --keep-order --xapply -j 4 "
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna.tsv \\
output/{2}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{2}/mrna_scored.tsv
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna.tsv \\
output/{3}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{3}/mrna_scored.tsv
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna.tsv \\
output/{4}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{4}/mrna_scored.tsv
" ::: NJU{6184..6215..4} ::: NJU{6185..6215..4} ::: NJU{6186..6215..4} ::: NJU{6187..6215..4}

parallel --keep-order --xapply -j 4 "
perl NJU_seq/mrna_analysis/score_neo.pl \\
output/{1}/mrna.tsv \\
output/{2}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{2}/mrna_scored_neo.tsv
perl NJU_seq/mrna_analysis/score_neo.pl \\
output/{1}/mrna.tsv \\
output/{3}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{3}/mrna_scored_neo.tsv
perl NJU_seq/mrna_analysis/score_neo.pl \\
output/{1}/mrna.tsv \\
output/{4}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{4}/mrna_scored_neo.tsv
" ::: NJU{6184..6215..4} ::: NJU{6185..6215..4} ::: NJU{6186..6215..4} ::: NJU{6187..6215..4}

parallel --keep-order --xapply -j 4 "
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna.tsv \\
output/{2}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{2}/mrna_scored.tsv
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna.tsv \\
output/{3}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{3}/mrna_scored.tsv
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna.tsv \\
output/{4}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{4}/mrna_scored.tsv
" ::: NJU{6220..6251..4} ::: NJU{6221..6251..4} ::: NJU{6222..6251..4} ::: NJU{6223..6251..4}

parallel --keep-order --xapply -j 4 "
perl NJU_seq/mrna_analysis/score_neo.pl \\
output/{1}/mrna.tsv \\
output/{2}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{2}/mrna_scored_neo.tsv
perl NJU_seq/mrna_analysis/score_neo.pl \\
output/{1}/mrna.tsv \\
output/{3}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{3}/mrna_scored_neo.tsv
perl NJU_seq/mrna_analysis/score_neo.pl \\
output/{1}/mrna.tsv \\
output/{4}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{4}/mrna_scored_neo.tsv
" ::: NJU{6220..6251..4} ::: NJU{6221..6251..4} ::: NJU{6222..6251..4} ::: NJU{6223..6251..4}

parallel --keep-order --xapply -j 1 "
cat output/{1}/mrna.cov
cat output/{2}/mrna.cov
cat output/{3}/mrna.cov
" ::: NJU{6185..6199..4} ::: NJU{6186..6199..4} ::: NJU{6187..6199..4}

parallel --keep-order --xapply -j 1 "
cat output/{1}/mrna.cov
cat output/{2}/mrna.cov
cat output/{3}/mrna.cov
" ::: NJU{6221..6251..4} ::: NJU{6222..6251..4} ::: NJU{6223..6251..4}

parallel --keep-order --xapply -j 1 '
cat output/{1}/mrna_scored.tsv > output/temp1_mrna_scored.tsv
cat output/{2}/mrna_scored.tsv > output/temp2_mrna_scored.tsv
cat output/{3}/mrna_scored.tsv > output/temp3_mrna_scored.tsv
for TOP in 50 100 300 500; do
awk -v a=`head -${TOP} output/temp1_mrna_scored.tsv | tail -1 | awk '\''{print $12}'\''` \
'\''$12>=a {print $1 $3 $2}'\'' output/temp1_mrna_scored.tsv \
>temp/sample1.txt
awk -v a=`head -${TOP} output/temp2_mrna_scored.tsv | tail -1 | awk '\''{print $12}'\''` \
'\''$12>=a {print $1 $3 $2}'\'' output/temp2_mrna_scored.tsv \
>temp/sample2.txt
awk -v a=`head -${TOP} output/temp3_mrna_scored.tsv | tail -1 | awk '\''{print $12}'\''` \
'\''$12>=a {print $1 $3 $2}'\'' output/temp3_mrna_scored.tsv \
>temp/sample3.txt
Rscript NJU_seq/presentation/point_venn.R \
Sample1 temp/sample1.txt \
Sample2 temp/sample2.txt \
Sample3 temp/sample3.txt \
output/{4}_mrna_top${TOP}_venn.png
rm temp/sample1.txt temp/sample2.txt temp/sample3.txt output/{4}_mrna_top${TOP}_venn.png*.log
done
' ::: NJU{6221..6235..4} ::: NJU{6222..6235..4} ::: NJU{6223..6235..4} ::: Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF

parallel --keep-order --xapply -j 1 "
perl extract_point.pl \
output/{1}/mrna_scored.tsv \
output/{2}/mrna_scored.tsv \
output/{3}/mrna_scored.tsv \
1000 1 >output/{4}_mrna_scored_500p.tsv
" ::: NJU{6185..6215..4} ::: NJU{6186..6215..4} ::: NJU{6187..6215..4} ::: Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_sixl Ath_bolt Ath_bloom Ath_fruit

parallel --keep-order --xapply -j 1 "
perl extract_point.pl \
output/{1}/mrna_scored.tsv \
output/{2}/mrna_scored.tsv \
output/{3}/mrna_scored.tsv \
500 1 >output/{4}_mrna_scored_500p.tsv
" ::: NJU{6221..6251..4} ::: NJU{6222..6251..4} ::: NJU{6223..6251..4} ::: Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_root Ath_stem Ath_flower Ath_leaf

parallel --keep-order --xapply -j 1 "
perl extract_point.pl \
output/{1}/mrna_scored_neo.tsv \
output/{2}/mrna_scored_neo.tsv \
output/{3}/mrna_scored_neo.tsv \
1000 1 >output/{4}_mrna_scored_500p.tmp
" ::: NJU{6221..6251..4} ::: NJU{6222..6251..4} ::: NJU{6223..6251..4} ::: Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_root Ath_stem Ath_flower Ath_leaf

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF; do
  perl NJU_seq/presentation/signature_count.pl \
    output/${PREFIX}_mrna_scored_500p.tsv \
    output/${PREFIX}_mrna_signature.pdf
done

parallel --keep-order --xapply -j 4 "
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna.tsv \\
output/{2}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{2}/mrna_scored.tsv
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna.tsv \\
output/{3}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{3}/mrna_scored.tsv
perl NJU_seq/mrna_analysis/score.pl \\
output/{1}/mrna.tsv \\
output/{4}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{4}/mrna_scored.tsv
" ::: NJU{6184..6215..4} ::: NJU{6185..6215..4} ::: NJU{6186..6215..4} ::: NJU{6187..6215..4}

for PREFIX in Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  perl NJU_seq/presentation/signature_count.pl \
    output/${PREFIX}_mrna_scored_500p.tsv \
    output/${PREFIX}_mrna_signature.pdf
done

parallel --keep-order --xapply -j 1 '
cat output/{1}/mrna_scored.tsv > output/temp1_mrna_scored.tsv
cat output/{2}/mrna_scored.tsv > output/temp2_mrna_scored.tsv
cat output/{3}/mrna_scored.tsv > output/temp3_mrna_scored.tsv
for TOP in 50 100 300 500; do
awk -v a=`head -${TOP} output/temp1_mrna_scored.tsv | tail -1 | awk '\''{print $12}'\''` \
'\''$12>=a {print $1 $3 $2}'\'' output/temp1_mrna_scored.tsv \
>temp/sample1.txt
awk -v a=`head -${TOP} output/temp2_mrna_scored.tsv | tail -1 | awk '\''{print $12}'\''` \
'\''$12>=a {print $1 $3 $2}'\'' output/temp2_mrna_scored.tsv \
>temp/sample2.txt
awk -v a=`head -${TOP} output/temp3_mrna_scored.tsv | tail -1 | awk '\''{print $12}'\''` \
'\''$12>=a {print $1 $3 $2}'\'' output/temp3_mrna_scored.tsv \
>temp/sample3.txt
Rscript NJU_seq/presentation/point_venn.R \
Sample1 temp/sample1.txt \
Sample2 temp/sample2.txt \
Sample3 temp/sample3.txt \
output/{4}_mrna_top${TOP}_venn.png
rm temp/sample1.txt temp/sample2.txt temp/sample3.txt output/{4}_mrna_top${TOP}_venn.png*.log
done
' ::: NJU{6221..6235..4} ::: NJU{6222..6235..4} ::: NJU{6223..6235..4} ::: Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF

parallel --keep-order --xapply -j 1 '
cat output/{1}/mrna_scored.tsv > output/temp1_mrna_scored.tsv
cat output/{2}/mrna_scored.tsv > output/temp2_mrna_scored.tsv
cat output/{3}/mrna_scored.tsv > output/temp3_mrna_scored.tsv
for TOP in 50 100 300 500; do
awk -v a=`head -${TOP} output/temp1_mrna_scored.tsv | tail -1 | awk '\''{print $12}'\''` \
'\''$12>=a {print $1 $3 $2}'\'' output/temp1_mrna_scored.tsv \
>temp/sample1.txt
awk -v a=`head -${TOP} output/temp2_mrna_scored.tsv | tail -1 | awk '\''{print $12}'\''` \
'\''$12>=a {print $1 $3 $2}'\'' output/temp2_mrna_scored.tsv \
>temp/sample2.txt
awk -v a=`head -${TOP} output/temp3_mrna_scored.tsv | tail -1 | awk '\''{print $12}'\''` \
'\''$12>=a {print $1 $3 $2}'\'' output/temp3_mrna_scored.tsv \
>temp/sample3.txt
Rscript NJU_seq/presentation/point_venn.R \
Sample1 temp/sample1.txt \
Sample2 temp/sample2.txt \
Sample3 temp/sample3.txt \
output/{4}_mrna_top${TOP}_venn.png
rm temp/sample1.txt temp/sample2.txt temp/sample3.txt output/{4}_mrna_top${TOP}_venn.png*.log
done
' ::: NJU{6277..6355..4} ::: NJU{6278..6355..4} ::: NJU{6279..6355..4} ::: Ath_ck_leaf_RF Ath_ck_root_RF Ath_ck_leaf Ath_ck_root Ath_Al_20mM_leaf_RF Ath_Al_20mM_leaf Ath_Al_20mM_root_RF Ath_Al_20mM_root Ath_Al_60mM_leaf_RF Ath_Al_60mM_leaf Ath_Al_60mM_root_RF Ath_Al_60mM_root Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_leaf Ath_Cd_5mM_root_RF Ath_Cd_5mM_root Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_leaf Ath_Cd_15mM_root_RF Ath_Cd_15mM_root

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  pigz -dc data/ath.gff3.gz |
    awk '$3=="gene"' |
    perl NJU_seq/tool/add_gene_name.pl \
      --id "gene_id=" \
      --name "Name=" \
      --col "8" \
      --file "output/${PREFIX}_mrna_scored_500p.tsv" \
      >output/${PREFIX}_mrna_scored_500p_name.tsv
done

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  pigz -dc data/ath.gff3.gz |
    awk '$3=="gene"' |
    perl NJU_seq/tool/add_gene_name.pl \
      --id "gene_id=" \
      --name "description=" \
      --col "9" \
      --file "output/${PREFIX}_mrna_scored_500p_name.tsv" \
      >output/${PREFIX}_mrna_scored_500p_name_dis.tsv
done

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/ath.fa.gz \
    output/${PREFIX}_mrna_scored_500p_name_dis.tsv 10 10 \
    >output/${PREFIX}_mrna_scored_500p_name_dis_motif1.tsv
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/ath.fa.gz \
    output/${PREFIX}_mrna_scored_500p_name_dis_motif1.tsv 20 20 \
    >output/${PREFIX}_mrna_scored_500p_name_dis_motif2.tsv
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/ath.fa.gz \
    output/${PREFIX}_mrna_scored_500p_name_dis_motif2.tsv 50 50 \
    >output/${PREFIX}_mrna_scored_500p_name_dis_motif3.tsv
done

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  mv output/${PREFIX}_mrna_scored_500p_name_dis_motif3.tsv output/${PREFIX}_mrna_500p.tsv
done

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  rm output/${PREFIX}_mrna_scored_500p_name_dis.tsv output/${PREFIX}_mrna_scored_500p_name.tsv
done

pigz -dc data/ath.gff3.gz |
  awk '(($3=="mRNA")&&($9~/biotype=protein_coding/))' |
  perl main_transcript_1.pl \
    --geneid "Parent=gene:" \
    --transid "transcript_id=" \
    >data/ath_main_transcript.txt

pigz -dc data/ath.gff3.gz |
  awk '(($3=="five_prime_UTR") || ($3=="three_prime_UTR") || ($3=="CDS"))' |
  perl NJU_seq/mrna_analysis/main_transcript_2.pl \
    --transid "Parent=transcript:" \
    --rep_trans "data/ath_main_transcript.txt" \
    >data/ath_transcript_region.tsv

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  perl NJU_seq/mrna_analysis/main_transcript_3.pl \
    data/ath_transcript_region.tsv \
    output/${PREFIX}_mrna_500p.tsv \
    output/${PREFIX}_mrna_500p.bak \
    >output/${PREFIX}_mrna_500p_dist.tsv
done

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  perl NJU_seq/mrna_analysis/main_transcript_4.pl \
    <output/${PREFIX}_mrna_500p_dist.tsv \
    >output/${PREFIX}_mrna_500p_dist_norm.tsv
  Rscript NJU_seq/presentation/point_distribution.R \
    output/${PREFIX}_mrna_500p_dist_norm.tsv \
    output/${PREFIX}_mrna_500p_distribution.pdf
done

pigz -dc data/ath.gff3.gz |
  grep "protein_coding" |
  awk '$3=="exon"' |
  perl NJU_seq/mrna_analysis/exon_distance_1.pl \
    data/ath_main_transcript.txt \
    >data/ath_transcript_exon.tsv

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  mv output/${PREFIX}_mrna_500p.bak output/${PREFIX}_mrna_500p.tsv
done

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  perl NJU_seq/mrna_analysis/exon_distance_2.pl \
    data/ath_transcript_exon.tsv \
    output/${PREFIX}_mrna_500p.tsv \
    output/${PREFIX}_mrna_500p_exon_site_bar.tsv \
    output/${PREFIX}_mrna_500p_exon_site_porta.tsv \
    >output/${PREFIX}_mrna_500p_exon_site.tsv
  Rscript NJU_seq/presentation/exon_distance.R \
    output/${PREFIX}_mrna_500p_exon_site_bar.tsv \
    output/${PREFIX}_mrna_500p_exon_site_porta.tsv \
    output/${PREFIX}_mrna_500p_exon_site.pdf
  rm output/${PREFIX}_mrna_500p_exon_site_bar.tsv
  rm output/${PREFIX}_mrna_500p_exon_site_porta.tsv
done

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  cat output/${PREFIX}_mrna_500p_exon_site.tsv |
    perl NJU_seq/mrna_analysis/exon_distance_patch.pl \
      >output/${PREFIX}_mrna_500p_exon_site_filtered.tsv
done

pigz -dc data/ath.gff3.gz |
  awk '$3=="CDS"{print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
  perl codon_distance_1.pl data/ath_main_transcript.txt \
    >data/ath_main_transcript_start_codon.tsv

pigz -dc data/ath.gff3.gz |
  awk '$3=="CDS"{print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
  perl codon_distance_1.pl data/ath_main_transcript.txt \
    >data/ath_main_transcript_stop_codon.tsv

pigz -dc data/ath.gff3.gz |
  awk '$3=="mRNA"||$3=="exon"' |
  perl codon_distance_2.pl \
    data/ath_main_transcript_start_codon.tsv \
    >data/ath_main_transcript_start_codon.yml

pigz -dc data/ath.gff3.gz |
  awk '$3=="mRNA"||$3=="exon"' |
  perl codon_distance_2.pl \
    data/ath_main_transcript_stop_codon.tsv \
    >data/ath_main_transcript_stop_codon.yml

for PREFIX in Ath_root_RF Ath_stem_RF Ath_flower_RF Ath_leaf_RF Ath_sixl_RF Ath_bolt_RF Ath_bloom_RF Ath_fruit_RF Ath_ck_leaf_RF Ath_ck_root_RF Ath_Al_20mM_leaf_RF Ath_Al_20mM_root_RF Ath_Al_60mM_leaf_RF Ath_Al_60mM_root_RF Ath_Cd_5mM_leaf_RF Ath_Cd_5mM_root_RF Ath_Cd_15mM_leaf_RF Ath_Cd_15mM_root_RF; do
  perl NJU_seq/mrna_analysis/codon_distance_3.pl \
    data/ath_main_transcript_start_codon.yml \
    output/${PREFIX}_mrna_500p.tsv \
    output/${PREFIX}_mrna_500p_start_codon_distance_bar.tsv \
    >output/${PREFIX}_mrna_500p_start_codon_distance.tsv
  perl NJU_seq/mrna_analysis/codon_distance_3.pl \
    data/ath_main_transcript_stop_codon.yml \
    output/${PREFIX}_mrna_500p.tsv \
    output/${PREFIX}_mrna_500p_stop_codon_distance_bar.tsv \
    >output/${PREFIX}_mrna_500p_stop_codon_distance.tsv
done

Rscript codon_distance_4group.R \
  output/Ath_root_RF_mrna_500p_start_codon_distance_bar.tsv \
  output/Ath_stem_RF_mrna_500p_start_codon_distance_bar.tsv \
  output/Ath_flower_RF_mrna_500p_start_codon_distance_bar.tsv \
  output/Ath_leaf_RF_mrna_500p_start_codon_distance_bar.tsv \
  output/Ath_root_RF_mrna_500p_stop_codon_distance_bar.tsv \
  output/Ath_stem_RF_mrna_500p_stop_codon_distance_bar.tsv \
  output/Ath_flower_RF_mrna_500p_stop_codon_distance_bar.tsv \
  output/Ath_leaf_RF_mrna_500p_stop_codon_distance_bar.tsv \
  output/Ath_tissue_mrna_500p_codon_distance.pdf

Rscript codon_distance_4group.R \
  output/Ath_sixl_RF_mrna_500p_start_codon_distance_bar.tsv \
  output/Ath_bolt_RF_mrna_500p_start_codon_distance_bar.tsv \
  output/Ath_bloom_RF_mrna_500p_start_codon_distance_bar.tsv \
  output/Ath_fruit_RF_mrna_500p_start_codon_distance_bar.tsv \
  output/Ath_sixl_RF_mrna_500p_stop_codon_distance_bar.tsv \
  output/Ath_bolt_RF_mrna_500p_stop_codon_distance_bar.tsv \
  output/Ath_bloom_RF_mrna_500p_stop_codon_distance_bar.tsv \
  output/Ath_fruit_RF_mrna_500p_stop_codon_distance_bar.tsv \
  output/Ath_develop_mrna_500p_codon_distance.pdf
