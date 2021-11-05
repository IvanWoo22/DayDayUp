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

for PREFIX in NJU{6412..6419}; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  perl NJU_seq/quality_control/pe_consistency.pl \
    data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
    temp/${PREFIX}.fq.gz
done

for PREFIX in NJU{6412..6419}; do
  mkdir output/${PREFIX}
  bsub -n 24 -o log/${PREFIX}_rrna_alignment.log -J "${PREFIX}" "
time bowtie2 -p 22 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/hsa_rrna -1 data/${PREFIX}/R1.fq.gz -2 data/${PREFIX}/R2.fq.gz -S output/${PREFIX}/rrna.raw.sam
"
done

for PREFIX in NJU{6412..6419}; do
  mkdir temp/${PREFIX}
  bsub -n 4 -o log/${PREFIX}_rrna_filter.log -J "${PREFIX}" "
    samtools view -h -f 97 -F 144 output/${PREFIX}/rrna.raw.sam > output/${PREFIX}/rrna.temp.sam
    samtools view -f 145 -F 96 output/${PREFIX}/rrna.raw.sam >> output/${PREFIX}/rrna.temp.sam
    samtools sort -n output/${PREFIX}/rrna.temp.sam | samtools view >output/${PREFIX}/rrna.filter.sam
    rm output/${PREFIX}/rrna.temp.sam
    pigz output/${PREFIX}/rrna.raw.sam
"
done

for PREFIX in NJU{6412..6419}; do
  bsub -n 4 -o log/${PREFIX}_rrna_judge.log -J "${PREFIX}" "
time awk '\$6!=\"*\"&&\$7==\"=\"{print \$1 \"\t\" \$3 \"\t\" \$4 \"\t\" \$6 \"\t\" \$10}' output/${PREFIX}/rrna.filter.sam |
perl NJU_seq/rrna_analysis/matchquality_judge.pl |
perl NJU_seq/rrna_analysis/multimatch_judge.pl \
>temp/${PREFIX}/rrna.out.tmp
"
done

for PREFIX in NJU{6412..6419}; do
  for CHR in 28s 18s 5-8s; do
    perl NJU_seq/rrna_analysis/readend_count.pl \
      NJU_seq/data/hsa_rrna/${CHR}.fa temp/${PREFIX}/rrna.out.tmp ${CHR} \
      >output/${PREFIX}/rrna_${CHR}.tsv
  done
done

for CHR in 28s 18s 5-8s; do
  perl NJU_seq/rrna_analysis/score.pl \
    output/NJU6415/rrna_${CHR}.tsv \
    output/NJU6412/rrna_${CHR}.tsv \
    output/NJU6413/rrna_${CHR}.tsv \
    output/NJU6414/rrna_${CHR}.tsv \
    >output/HeLa_G1_rrna_${CHR}_scored.tsv
done

for CHR in 28s 18s 5-8s; do
  perl NJU_seq/rrna_analysis/score.pl \
    output/NJU6419/rrna_${CHR}.tsv \
    output/NJU6416/rrna_${CHR}.tsv \
    output/NJU6417/rrna_${CHR}.tsv \
    output/NJU6418/rrna_${CHR}.tsv \
    >output/HeLa_G2_rrna_${CHR}_scored.tsv
done
