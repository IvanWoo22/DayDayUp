parallel --xapply -j 15 "
pigz -dc CpG_sites_merge.tsv.gz |
awk -va={1} 'NR==1{for(i=1;i<=2;i++){printf \"\\t\" \$i};
for(i=a;i<=(a+13);i++){printf \"\\t\" \$i}; print \"\"};
NR>1{for(i=1;i<=3;i++) printf \$i \"\\t\";
for(i=a+1;i<=a+13;i++) printf \$i \"\\t\"; print \$(a+14)}' > raw_P{2}_merge.tsv
" ::: {3..199..14} ::: {1..15}

parallel --xapply -j 7 "
pigz -dc CpG_sites_merge.tsv.gz |
awk -va={1} 'NR==1{printf \"\\t\" \$1 \"\\t\" \$2;
for(i=a;i<=NF;i+=14) printf \"\\t\" \$i \"\\t\" \$(i+1); print \"\"};
NR>1{printf \$1 \"\\t\" \$2 \"\\t\" \$3;
for(i=a+1;i<=NF;i+=14) printf \"\\t\" \$i \"\\t\" \$(i+1); print \"\"}' > raw_L{2}_merge.tsv
" ::: {3..15..2} ::: {1..7}

parallel --keep-order --xapply -j 15 "
awk '{for(i=5;i<=17;i=i+2){if(\$i!=\"NA\"&&\$i>=15){sum[(i-3)/2]++}}};
END{print sum[1] \"\\t\" sum[2] \"\\t\" sum[3] \"\t\" sum[4] \"\\t\" sum[5] \"\\t\" sum[6] \"\\t\" sum[7]}' raw_P{}_merge.tsv
" ::: {1..15} >raw_site_each_sample.tsv

pigz -dc CpG_sites_merge.tsv.gz |
  awk '{sum=0;
    {for(i=5;i<=NF;i=i+2){if($i!="NA"&&$i>=15){sum++}}};
    {a[sum]++}};END{for(i=1;i<=105;i++){print a[i]}}' >each_sample.tsv

for P in {1..15}; do
  parallel --xapply --keep-order -j 7 "
    awk -vk={} '{sum=0;{for(i=5;i<=17;i=i+2){if(\$i!=\"NA\"&&\$i>=15){sum++}}};
    {if(sum==k){print}}}' raw_P${P}_merge.tsv | wc -l
    " ::: {7..1}
done >each_loc_each_sample.tsv

for L in {1..7}; do
  parallel --xapply --keep-order -j 15 "
    awk -vk={} '{sum=0;{for(i=5;i<=NF;i=i+2){if(\$i!=\"NA\"&&\$i>=15){sum++}}};
    {if(sum==k){print}}}' raw_L${L}_merge.tsv | wc -l
    " ::: {15..1}
done >each_sample_each_loc.tsv

pigz -dc CpG_sites_merge.tsv.gz |
  awk 'NR==1{for(i=1;i<=16;i++){printf "\t" $i};
    for(i=59;i<=NF;i++){printf "\t" $i}; printf "\n"};
    NR>1{for(i=1;i<=17;i++){printf $i "\t"};
    for(i=60;i<NF;i++){printf $i "\t"}; printf $NF "\n"}' |
  pigz >CpG_sites_merge_3loss.tsv.gz

pigz -dc CpG_sites_merge_3loss.tsv.gz |
  awk '{sum=0;
    {for(i=5;i<=NF;i=i+2){if($i!="NA"&&$i>=15){sum++}}};
    {a[sum]++}};END{for(i=1;i<=84;i++){print a[i]}}' >each_sample_3loss.tsv

pigz -dc CpG_sites_merge_3loss.tsv.gz |
  awk -vk=68 '{sum=0;{for(i=5;i<=NF;i=i+2){if($i!="NA"&&$i>=15){sum++}}};{if(sum>=k){print}}}' |
  pigz >CpG_sites_merge_3loss_68sample.tsv.gz

parallel --xapply -j 12 "
  pigz -dc CpG_sites_merge_3loss_68sample.tsv.gz |
    awk -va={1} 'NR==1{for(i=1;i<=2;i++){printf \"\\t\" \$i};
  for(i=a;i<=(a+13);i++){printf \"\\t\" \$i}; printf \"\\n\"};
  NR>1{for(i=1;i<=3;i++) printf \$i \"\\t\";
  for(i=a+1;i<=a+13;i++) printf \$i \"\\t\"; print \$(a+14)}' >neo_P{2}_merge.tsv
" ::: {3..157..14} ::: {1..12}

parallel --keep-order --xapply -j 12 "
awk '{for(i=5;i<=17;i=i+2){if(\$i!=\"NA\"&&\$i>=15){sum[(i-3)/2]++}}};
END{print sum[1] \"\\t\" sum[2] \"\\t\" sum[3] \"\\t\" sum[4] \"\\t\" sum[5] \"\\t\" sum[6] \"\\t\" sum[7]}' neo_P{}_merge.tsv
" ::: {1..12} >filter_site_each_sample.tsv

pigz -dc CpG_sites_merge_3loss_68sample.tsv.gz |
  awk '$2!="chrM"{print $2 ":" $3}' |
  spanr cover stdin >RRBS_sites_merge_filterM.yml

parallel --keep-order --xapply -j 18 "
  V=\$(({}*2+3)); pigz -dc CpG_sites_merge_3loss_68sample.tsv.gz |
  awk -va=\${V} 'NR>1&&\$a>=15{print \$2 \"\\t\" \$3}' > neo_S{}_merge.tsv
" ::: {1..84}

parallel --keep-order --xapply -j 18 "
awk '\$1!=\"chrM\"{print \$1 \":\" \$2}' neo_S{}_merge.tsv | spanr cover stdin > neo_S{}_merge.yml
" ::: {1..84}

parallel --xapply --keep-order -j 12 "
 spanr compare --op intersect {}.yml eu_open1.yml |
 spanr stat chr.sizes stdin --all |
 awk -F ',' 'NR==2{printf \$2 \"\\t\"}'
 spanr compare --op intersect {}.yml eu_open2.yml |
 spanr stat chr.sizes stdin --all |
 awk -F ',' 'NR==2{printf \$2 \"\\t\"}'
 spanr compare --op intersect {}.yml eu_open3.yml |
 spanr stat chr.sizes stdin --all |
 awk -F ',' 'NR==2{printf \$2 \"\\t\"}'
 spanr compare --op intersect {}.yml het_open1.yml |
 spanr stat chr.sizes stdin --all |
 awk -F ',' 'NR==2{printf \$2 \"\\t\"}'
 spanr compare --op intersect {}.yml het_open2.yml |
 spanr stat chr.sizes stdin --all |
 awk -F ',' 'NR==2{printf \$2 \"\\t\"}'
 spanr compare --op intersect {}.yml het_open3.yml |
 spanr stat chr.sizes stdin --all |
 awk -F ',' 'NR==2{printf \$2 \"\\t\"}'
 echo
 " ::: neo_S{1..84}_merge

parallel --xapply --keep-order -j 12 "
  for j in cpgi_unmasked exon enhancer promoter ctcf; do
    for i in open{1..3}; do
      awk 'BEGIN{sum=0};{sum+=\$4};END{printf sum \"\\t\"}' \${i}_\${j}/\${i}_\${j}.withS{}_merge.tsv
    done
  done
  echo
  " ::: {1..84}

parallel --xapply -j 15 "awk '{sum=0;{for(i=5;i<=17;i=i+2){if(\$i!=\"NA\"&&\$i>=15){sum++}}};
{if(sum==7){print}}}' raw_P{}_merge.tsv >filter15_P{}_merge.tsv" ::: {1..15}

perl venn_site.pl filter15_P{{5..10},{12..15}}_merge.tsv \
  >filter_site_venn_merge.tsv

awk 'sum=0;{for(i=4;i<=NF;i++)sum+=$i;if(sum>=5){printf $1;for(j=4;j<=NF;j++)printf "\t" $j;printf "\n"}}' filter_site_venn_merge.tsv |
  sort -nk 1,1 >site_venn_merge.tsv

for i in NJU61{48..51}; do
  bsub -n 24 -J "${i}" "
    parallel --pipe --block 1000M --no-run-if-empty --linebuffer --keep-order -j 16 '
        perl NJU_seq/mrna_analysis/dedup.pl --refstr \"Parent=transcript:\" --transid \"ENST\" --info data/hsa_exon.info
      ' <output/${i}/mrnas.out.tmp |
      perl NJU_seq/mrna_analysis/dedup.pl \
        --refstr \"Parent=transcript:\" \
        --transid \"ENST\" \
        --info data/hsa_exon.info \
        >output/${i}/mrna.dedup.tmp
    "
done

bsub -n 24 -J "almostunique" '
parallel --xapply -j 4 "
time bash NJU_seq/mrna_analysis/almostunique.sh output/{}/mrnas.dedup.tmp data/{}/R1.mrnas.fq.gz output/{} output/{}/mrnas.almostunique.tmp
" ::: NJU61{48..51}
'

bsub -n 24 -J "almostunique" '
parallel --xapply -j 4 "
time perl NJU_seq/mrna_analysis/count.pl output/{}/mrnas.almostunique.tmp >output/{}/mrnas.count.tmp
" ::: NJU61{48..51}
'

for PREFIX in NJU61{48..51}; do
  time pigz -dc data/hsa.gff3.gz |
    awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
    perl NJU_seq/mrna_analysis/merge.pl \
      --refstr "Parent=" \
      --geneid "ENSG" \
      --transid "ENST" \
      -i output/"${PREFIX}"/mrnas.count.tmp \
      -o output/"${PREFIX}"/mrnas.tsv
done

perl NJU_seq/mrna_analysis/score.pl \
  output/NJU6148/mrnas.tsv \
  output/NJU6149/mrnas.tsv \
  output/NJU6150/mrnas.tsv \
  output/NJU6151/mrnas.tsv \
  >output/HeLa_RF_mrnas_Nm_score.tsv
