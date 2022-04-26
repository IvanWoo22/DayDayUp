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

parallel --xapply -j 15 "awk '{sum=0;{for(i=5;i<=17;i=i+2){if(\$i!=\"NA\"&&\$i>=15){sum++}}};
{if(sum==7){print}}}' raw_P{}_merge.tsv >filter15_P{}_merge.tsv" ::: {1..15}
