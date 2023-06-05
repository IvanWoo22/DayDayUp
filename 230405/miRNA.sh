for i in NJU1013 NJU1014; do
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
    -O 6 -m 10 --discard-untrimmed -e 0.1 -j 8 \
    ${i}/*_R1_001.fastq.gz ${i}/*_R2_001.fastq.gz \
    -o ${i}/R1.fq.gz -p ${i}/R2.fq.gz \
    >${i}/cutadapt.log 2>&1
  pear -j 4 -f ${i}/R1.fq.gz -r ${i}/R2.fq.gz -o ${i} -n 10
  perl ../../perbool/fastq2count.pl <${i}.assembled.fastq >${i}/reads.count
  perl mirna_count.pl ${i}/reads.count hsa_mirna.fa >${i}/mir.count.tsv
done

for i in NJU10{21..46}; do
  mkdir ${i}
  ln -sf /home/ivan/cater_data/${i}/${i}*_R1_*.fastq.gz /home/ivan/fat/miR_Nm/${i}/R1.fq.gz
  ln -sf /home/ivan/cater_data/${i}/${i}*_R2_*.fastq.gz /home/ivan/fat/miR_Nm/${i}/R2.fq.gz
  cutadapt -a AGATCGGAAGAGC -A GATCGTCGGACTGTAG \
    -O 6 -m 10 -e 0.1 -j 8 \
    ${i}/R1.fq.gz ${i}/R2.fq.gz \
    -o ${i}/clean_R1.fq.gz \
    -p ${i}/clean_R2.fq.gz \
    >${i}/cutadapt.log 2>&1
  pear -j 4 -f ${i}/R1.fq.gz -r ${i}/R2.fq.gz -o ${i} -n 10
  perl ../../perbool/fastq2count.pl <${i}.assembled.fastq >${i}/reads.count
  perl mirna_count.pl ${i}/reads.count hsa_mirna.fa >${i}/mir.count.tsv
  rm ${i}.assembled.fastq ${i}.forword.fastq ${i}.reverse.fastq
done
