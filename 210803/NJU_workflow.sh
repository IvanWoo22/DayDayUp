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
  perl NmMark.pl ${CHR}_Nm.txt output/HeLa_HC_rrna_${CHR}_scored.tsv >output/HeLa_HC_rrna_${CHR}.tsv
  perl NmMark.pl ${CHR}_Nm.txt output/HeLa_sR_rrna_${CHR}_scored.tsv >output/HeLa_sR_rrna_${CHR}.tsv
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

bsub -n 24 -J "NJU6176" 'time bowtie2 -p 20 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/mmu_rrna -1 data/NJU6176/R1.fq.gz -2 data/NJU6176/R2.fq.gz -S output/NJU6176/rrna.raw.sam 2>&1 | tee output/NJU6176/rrna.bowtie2_1114.log'
bsub -n 24 -J "NJU6177" 'time bowtie2 -p 20 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/mmu_rrna -1 data/NJU6177/R1.fq.gz -2 data/NJU6177/R2.fq.gz -S output/NJU6177/rrna.raw.sam 2>&1 | tee output/NJU6177/rrna.bowtie2_1114.log'
bsub -n 24 -J "NJU6178" 'time bowtie2 -p 20 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/mmu_rrna -1 data/NJU6178/R1.fq.gz -2 data/NJU6178/R2.fq.gz -S output/NJU6178/rrna.raw.sam 2>&1 | tee output/NJU6178/rrna.bowtie2_1114.log'
bsub -n 24 -J "NJU6179" 'time bowtie2 -p 20 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/mmu_rrna -1 data/NJU6179/R1.fq.gz -2 data/NJU6179/R2.fq.gz -S output/NJU6179/rrna.raw.sam 2>&1 | tee output/NJU6179/rrna.bowtie2_1114.log'
bsub -n 24 -J "NJU6180" 'time bowtie2 -p 20 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/mmu_rrna -1 data/NJU6180/R1.fq.gz -2 data/NJU6180/R2.fq.gz -S output/NJU6180/rrna.raw.sam 2>&1 | tee output/NJU6180/rrna.bowtie2_1114.log'
bsub -n 24 -J "NJU6181" 'time bowtie2 -p 20 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/mmu_rrna -1 data/NJU6181/R1.fq.gz -2 data/NJU6181/R2.fq.gz -S output/NJU6181/rrna.raw.sam 2>&1 | tee output/NJU6181/rrna.bowtie2_1114.log'
bsub -n 24 -J "NJU6182" 'time bowtie2 -p 20 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/mmu_rrna -1 data/NJU6182/R1.fq.gz -2 data/NJU6182/R2.fq.gz -S output/NJU6182/rrna.raw.sam 2>&1 | tee output/NJU6182/rrna.bowtie2_1114.log'
bsub -n 24 -J "NJU6183" 'time bowtie2 -p 20 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/mmu_rrna -1 data/NJU6183/R1.fq.gz -2 data/NJU6183/R2.fq.gz -S output/NJU6183/rrna.raw.sam 2>&1 | tee output/NJU6183/rrna.bowtie2_1114.log'

bsub -n 80 -q fat_768 -J "NJU6176" 'time bowtie2 -p 70 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 --score-min C,0,0 --xeq -x index/n2a_basic_protein_coding -1 data/NJU6176/R1.mrna.fq.gz -2 data/NJU6176/R2.mrna.fq.gz -S output/NJU6176/mrna.raw.sam 2>&1 | tee output/NJU6176/mrna.bowtie2_1114.log'
bsub -n 80 -q fat_768 -J "NJU6177" 'time bowtie2 -p 70 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 --score-min C,0,0 --xeq -x index/n2a_basic_protein_coding -1 data/NJU6177/R1.mrna.fq.gz -2 data/NJU6177/R2.mrna.fq.gz -S output/NJU6177/mrna.raw.sam 2>&1 | tee output/NJU6177/mrna.bowtie2_1114.log'
bsub -n 80 -q fat_768 -J "NJU6178" 'time bowtie2 -p 70 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 --score-min C,0,0 --xeq -x index/n2a_basic_protein_coding -1 data/NJU6178/R1.mrna.fq.gz -2 data/NJU6178/R2.mrna.fq.gz -S output/NJU6178/mrna.raw.sam 2>&1 | tee output/NJU6178/mrna.bowtie2_1114.log'
bsub -n 80 -q fat_768 -J "NJU6179" 'time bowtie2 -p 70 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 --score-min C,0,0 --xeq -x index/n2a_basic_protein_coding -1 data/NJU6179/R1.mrna.fq.gz -2 data/NJU6179/R2.mrna.fq.gz -S output/NJU6179/mrna.raw.sam 2>&1 | tee output/NJU6179/mrna.bowtie2_1114.log'
bsub -n 80 -q fat_768 -J "NJU6180" 'time bowtie2 -p 70 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 --score-min C,0,0 --xeq -x index/n2a_basic_protein_coding -1 data/NJU6180/R1.mrna.fq.gz -2 data/NJU6180/R2.mrna.fq.gz -S output/NJU6180/mrna.raw.sam 2>&1 | tee output/NJU6180/mrna.bowtie2_1114.log'
bsub -n 80 -q fat_768 -J "NJU6181" 'time bowtie2 -p 70 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 --score-min C,0,0 --xeq -x index/n2a_basic_protein_coding -1 data/NJU6181/R1.mrna.fq.gz -2 data/NJU6181/R2.mrna.fq.gz -S output/NJU6181/mrna.raw.sam 2>&1 | tee output/NJU6181/mrna.bowtie2_1114.log'
bsub -n 80 -q fat_768 -J "NJU6182" 'time bowtie2 -p 70 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 --score-min C,0,0 --xeq -x index/n2a_basic_protein_coding -1 data/NJU6182/R1.mrna.fq.gz -2 data/NJU6182/R2.mrna.fq.gz -S output/NJU6182/mrna.raw.sam 2>&1 | tee output/NJU6182/mrna.bowtie2_1114.log'
bsub -n 80 -q fat_768 -J "NJU6183" 'time bowtie2 -p 70 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 --score-min C,0,0 --xeq -x index/n2a_basic_protein_coding -1 data/NJU6183/R1.mrna.fq.gz -2 data/NJU6183/R2.mrna.fq.gz -S output/NJU6183/mrna.raw.sam 2>&1 | tee output/NJU6183/mrna.bowtie2_1114.log'

bsub -n 24 -J "samtools1" "
for PREFIX in NJU{6176..6177}; do
  samtools view -@ 20 -h -f 97 -F 144 \
    output/${PREFIX}/mrna.raw.sam >output/${PREFIX}/mrna.temp.sam \
    2>output/${PREFIX}/mrna_filter.log
  samtools view -@ 20 -f 145 -F 96 \
    output/${PREFIX}/mrna.raw.sam >>output/${PREFIX}/mrna.temp.sam \
    2>>output/${PREFIX}/mrna_filter.log
  samtools sort -@ 12 -n output/${PREFIX}/mrna.temp.sam |
    samtools view -@ 8 >output/${PREFIX}/mrna.filter.sam \
      2>>output/${PREFIX}/mrna_filter.log
  rm output/${PREFIX}/mrna.temp.sam
  pigz output/${PREFIX}/mrna.raw.sam
done"

bsub -n 24 -J "samtools1" "
for PREFIX in NJU{6178..6179}; do
  samtools view -@ 20 -h -f 97 -F 144 \
    output/${PREFIX}/mrna.raw.sam >output/${PREFIX}/mrna.temp.sam \
    2>output/${PREFIX}/mrna_filter.log
  samtools view -@ 20 -f 145 -F 96 \
    output/${PREFIX}/mrna.raw.sam >>output/${PREFIX}/mrna.temp.sam \
    2>>output/${PREFIX}/mrna_filter.log
  samtools sort -@ 12 -n output/${PREFIX}/mrna.temp.sam |
    samtools view -@ 8 >output/${PREFIX}/mrna.filter.sam \
      2>>output/${PREFIX}/mrna_filter.log
  rm output/${PREFIX}/mrna.temp.sam
  pigz output/${PREFIX}/mrna.raw.sam
done"

time gzip -dcf output/"${PREFIX}"/mrna.raw.sam.gz |
  parallel --pipe --block 1G --no-run-if-empty --linebuffer --keep-order -j "${THREAD}" '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}
    '\'' |
    perl NJU_seq/mrna_analysis/multimatch_judge.pl
  ' | perl NJU_seq/mrna_analysis/multimatch_judge.pl \
    >temp/"${PREFIX}"/mrna.out.tmp

cat temp/"${PREFIX}"/mrna.out.tmp |
  parallel --pipe --block 1G --no-run-if-empty --linebuffer --keep-order -j "${THREAD}" '
    perl NJU_seq/mrna_analysis/dedup.pl \
      --refstr "Parent=" \
      --transid "ENSMUST" \
      --info data/mmu_exon.info
  ' |
  perl NJU_seq/mrna_analysis/dedup.pl \
    --refstr "Parent=" \
    --transid "ENSMUST" \
    --info data/mmu_exon.info \
    >temp/"${PREFIX}"/mrna.dedup.tmp

time bash NJU_seq/mrna_analysis/almostunique.sh \
  temp/"${PREFIX}"/mrna.dedup.tmp \
  data/"${PREFIX}"/R1.mrna.fq.gz \
  temp/"${PREFIX}" \
  temp/"${PREFIX}"/mrna.almostunique.tmp

