for i in NJU81{{00..07},16,17}; do
  mkdir ${i}
done

for i in NJU81{{00..07},16,17}; do
  ln -sf /home/ivan/cater_data/${i}/${i}_R1.fq.gz /home/ivan/fat/APOBEC_EDIT/${i}/R1.fq.gz
  ln -sf /home/ivan/cater_data/${i}/${i}_R2.fq.gz /home/ivan/fat/APOBEC_EDIT/${i}/R2.fq.gz
done

parallel --keep-order --xapply -j 4 \
  'pear -j 4 -f {}/R1.fq.gz -r {}/R2.fq.gz -o {}/merge' \
  ::: NJU81{{00..07},16,17}

parallel --keep-order --xapply -j 4 \
  'perl ~/perbool/fastq2count.pl <{}/merge.assembled.fastq | sort -nrk 2 >{}/count.tsv' \
  ::: NJU81{{00..07},16,17}

for i in NJU81{{00..07},16,17}; do
  parallel --keep-order --xapply -j 8 \
    "awk -va={1} -vb={2} '\$1~a&&\$1~b' ${i}/count.tsv >${i}/{3}.count.tsv" \
    ::: ^GTCAGTCA ^CCTTCCAT ^AGGAACAC ^CTTACAGC ^TACCTGCA ^GTCAGTCA ^GTCAGTCA ^GTCAGTCA ^GTCAGTCA ^CCTTCCAT ^CCTTCCAT ^CCTTCCAT ^AGGAACAC \
    ::: CGGCATTA$ CACGCAAT$ GGAATGTC$ TGGTGAAG$ GGACATCA$ CACGCAAT$ GGAATGTC$ TGGTGAAG$ GGACATCA$ CGGCATTA$ GGAATGTC$ TGGTGAAG$ CACGCAAT$ \
    ::: AA BB CC DD EE AB AC AD AE BA BC BD CB
done

perl C2CT.pl GGGGCGGACCGCGTGCGCTCGGCGGCTG >search.list.tsv
perl C2CT.pl GCTGGCCACGGCCGCGGCCCGGGGTC >>search.list.tsv

for i in NJU81{{00..07},16,17}; do
  for j in AA BB CC DD EE AB AC AD AE BA BC BD CB; do
    parallel --keep-order --xapply --colsep '\t' -j 8 \
      "awk -va={1} -vb={2} 'BEGIN{sum1=0;sum2=0};\$1~a{sum1=sum1+\$2};\$1~b{sum2=sum2+\$2};END{print sum1 \"\\t\" sum2 }' ${i}/${j}.count.tsv" \
      <search.list.tsv >${i}/${j}.cov
    echo -e "${i}\t${j}\t$(awk 'BEGIN{sum=0};{sum=sum+$2};END{print sum}' ${i}/${j}.count.tsv)" >>temp1.tsv
  done
done

for i in NJU81{{00..07},16,17}; do
  for j in AA BB CC DD EE AB AC AD AE BA BC BD CB; do
    awk -va=${i} -vb=${j} 'BEGIN{printf a "\t" b};{printf "\t" $1 "\t" $2 "\t" ($2/($1+$2+1)) "\t"};END{printf "\n"}' ${i}/${j}.cov
  done
done >temp2.tsv
tsv-join temp1.tsv -f temp2.tsv -k 1,2 -a 3-77 >out0315.tsv
