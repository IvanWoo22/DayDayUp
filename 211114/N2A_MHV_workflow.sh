#!/usr/bin/env bash

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

bsub -n 24 -J "sam_filter" "
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

time perl NJU_seq/mrna_analysis/count.pl \
  temp/"${PREFIX}"/mrna.almostunique.tmp \
  >temp/"${PREFIX}"/mrna.count.tmp

time gzip -dcf data/mmu.gff3.gz |
  awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
  perl NJU_seq/mrna_analysis/merge.pl \
    --refstr "Parent=" \
    --geneid "ENSMUSG" \
    --transid "ENSMUST" \
    -i temp/"${PREFIX}"/mrna.count.tmp \
    -o output/"${PREFIX}"/mrna.tsv

pigz -dc data/mmu.gff3.gz |
  awk '$3=="gene"' |
  perl NJU_seq/tool/add_gene_name.pl \
    --id "gene_id=" --name "gene_name=" \
    --col "8" --file "output/N2A_MHV_mrna_all_scored.tsv" \
    >output/N2A_MHV_mrna_all_scored_name.tsv

pigz -dc data/mmu_ens.gff3.gz |
  awk '$3=="gene"' |
  perl NJU_seq/tool/add_gene_name.pl \
    --id "gene_id=" --name "description=" \
    --col "9" --file "output/N2A_MHV_mrna_all_scored_name.tsv" \
    >output/N2A_MHV_mrna_all_scored_name_dis.tsv

pigz -dc data/mmu.gff3.gz |
  awk '$3=="gene"' |
  perl NJU_seq/tool/add_gene_name.pl \
    --id "gene_id=" --name "gene_name=" \
    --col "8" --file "output/N2A_HC_mrna_all_scored.tsv" \
    >output/N2A_HC_mrna_all_scored_name.tsv

pigz -dc data/mmu_ens.gff3.gz |
  awk '$3=="gene"' |
  perl NJU_seq/tool/add_gene_name.pl \
    --id "gene_id=" --name "description=" \
    --col "9" --file "output/N2A_HC_mrna_all_scored_name.tsv" \
    >output/N2A_HC_mrna_all_scored_name_dis.tsv