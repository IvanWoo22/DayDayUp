#!/usr/bin/env bash

for PREFIX in NJU{6420..6427}; do
  mkdir ${PREFIX}
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
    -O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1.fq.gz -p ${PREFIX}/R2.fq.gz \
    ${PREFIX}_*R1*.gz ${PREFIX}_*R2*.gz -j 12 | tee ../log/${PREFIX}_cutadapt.log
done

for PREFIX in NJU{6420..6427}; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  perl ~/NJU_seq/quality_control/pe_consistency.pl \
    data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
    temp/${PREFIX}.fq.gz
done

time perl ~/NJU_seq/quality_control/fastq_qc.pl \
  temp/NJU6420.fq.gz \
  temp/NJU6421.fq.gz \
  temp/NJU6422.fq.gz \
  temp/NJU6423.fq.gz \
  temp/NJU6424.fq.gz \
  temp/NJU6425.fq.gz \
  temp/NJU6426.fq.gz \
  temp/NJU6427.fq.gz \
  output \
  HeLa_FBL

for PREFIX in NJU{6420..6427}; do
  mkdir output/${PREFIX}
  time bowtie2 -p 20 -a -t \
    --end-to-end -D 20 -R 3 \
    -N 0 -L 10 -i S,1,0.50 --np 0 \
    --xeq -x index/hsa_rrna \
    -1 data/${PREFIX}/R1.fq.gz \
    -2 data/${PREFIX}/R2.fq.gz \
    -S output/${PREFIX}/rrna.raw.sam 2>output/${PREFIX}.rrna_align.log
done

for PREFIX in NJU{6420..6427}; do
  samtools view -@ 20 -h -f 97 -F 144 \
    output/${PREFIX}/rrna.raw.sam >output/${PREFIX}/rrna.temp.sam \
    2>output/${PREFIX}.rrna_filter.log
  samtools view -@ 20 -f 145 -F 96 \
    output/${PREFIX}/rrna.raw.sam >>output/${PREFIX}/rrna.temp.sam \
    2>>output/${PREFIX}.rrna_filter.log
  samtools sort -@ 12 -n output/${PREFIX}/rrna.temp.sam |
    samtools view -@ 8 >output/${PREFIX}/rrna.filter.sam \
      2>>output/${PREFIX}.rrna_filter.log
  rm output/${PREFIX}/rrna.temp.sam
  pigz output/${PREFIX}/rrna.raw.sam
done

for PREFIX in NJU{6420..6427}; do
  time awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' output/${PREFIX}/rrna.filter.sam |
    perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl \
      >temp/${PREFIX}/rrna.out.tmp
done

for PREFIX in NJU{6420..6427}; do
  for CHR in 28s 18s 5-8s; do
    perl ~/NJU_seq/rrna_analysis/readend_count.pl \
      ~/NJU_seq/data/hsa_rrna/${CHR}.fa temp/${PREFIX}/rrna.out.tmp ${CHR} \
      >output/${PREFIX}/rrna_${CHR}.tsv
  done
done

for CHR in 28s 18s 5-8s; do
  perl ~/NJU_seq/rrna_analysis/score.pl \
    output/NJU6420/rrna_${CHR}.tsv \
    output/NJU6421/rrna_${CHR}.tsv \
    output/NJU6422/rrna_${CHR}.tsv \
    output/NJU6423/rrna_${CHR}.tsv \
    >output/HeLa_HC_rrna_${CHR}_scored.tsv
done

for CHR in 28s 18s 5-8s; do
  perl ~/NJU_seq/rrna_analysis/score.pl \
    output/NJU6424/rrna_${CHR}.tsv \
    output/NJU6425/rrna_${CHR}.tsv \
    output/NJU6426/rrna_${CHR}.tsv \
    output/NJU6427/rrna_${CHR}.tsv \
    >output/HeLa_sR_rrna_${CHR}_scored.tsv
done

for CHR in 28s 18s; do
  perl NmMark.pl ${CHR}_Nm.txt output/HeLa_HC_rrna_${CHR}_scored.tsv > output/HeLa_HC_rrna_${CHR}.tsv
  perl NmMark.pl ${CHR}_Nm.txt output/HeLa_sR_rrna_${CHR}_scored.tsv > output/HeLa_sR_rrna_${CHR}.tsv
done


for PREFIX in NJU{6420..6425}; do
  kraken2 --threads 16 --use-names \
    --gzip-compressed --paired \
    --report output/${PREFIX}/class.report \
    --classified-out output/${PREFIX}/class#.fq \
    --db index/human_bacteria_viral_fungi \
    data/${PREFIX}/R1.fq.gz \
    data/${PREFIX}/R2.fq.gz \
    --output output/${PREFIX}/class.tsv
done
