#!/usr/bin/env bash

for PREFIX in HeLa A549 293T; do
  pigz -dc data/hsa.gff3.gz |
    awk '$3=="gene"' | perl NJU_seq/tool/add_gene_name.pl \
      --id "gene_id=" --name "gene_name=" \
      --col "8" --file "output/Hsa_${PREFIX}_mrna_all_1.tsv" \
      >output/Hsa_${PREFIX}_mrna_all_2.tsv
  pigz -dc data/hsa_ens.gff3.gz |
    awk '$3=="gene"' | perl NJU_seq/tool/add_gene_name.pl \
      --id "gene_id=" --name "description=" \
      --col "9" --file "output/Hsa_${PREFIX}_mrna_all_2.tsv" \
      >output/Hsa_${PREFIX}_mrna_all_3.tsv
done

for PREFIX in HeLa A549 293T; do
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/GRCm38.primary_assembly.genome.fa output/Hsa_${PREFIX}_mrna_all_3.tsv \
    10 10 >output/Hsa_${PREFIX}_mrna_all_4.tsv
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/GRCm38.primary_assembly.genome.fa output/Hsa_${PREFIX}_mrna_all_4.tsv \
    20 20 >output/Hsa_${PREFIX}_mrna_all_5.tsv
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/GRCm38.primary_assembly.genome.fa output/Hsa_${PREFIX}_mrna_all_5.tsv \
    50 50 >output/Hsa_${PREFIX}_mrna_all_6.tsv
done

for PREFIX in HeLa A549 293T; do
  perl NJU_seq/mrna_analysis/main_transcript_3.pl \
    data/hsa_basic_transcript_region.tsv \
    output/Hsa_${PREFIX}_mrna_all_6.tsv \
    output/Hsa_${PREFIX}_mrna_all_7.tsv \
    >output/Hsa_${PREFIX}_mrna_all_6.bak
done

for PREFIX in HeLa A549 293T; do
  perl NJU_seq/mrna_analysis/main_transcript_4.pl \
    <output/Hsa_${PREFIX}_mrna_all_6.bak \
    >output/Hsa_${PREFIX}_mrna_all_6_norm.bak
  Rscript NJU_seq/presentation/point_distribution.R \
    output/Hsa_${PREFIX}_mrna_all_6_norm.bak \
    output/Hsa_${PREFIX}_mrna_all_distribution.pdf
done

for PREFIX in HeLa A549 293T; do
  perl NJU_seq/mrna_analysis/exon_distance_2.pl \
    data/hsa_basic_transcript_exon.tsv \
    output/Hsa_${PREFIX}_mrna_all_7.tsv \
    output/Hsa_${PREFIX}_mrna_all_exon_site_bar.tsv \
    output/Hsa_${PREFIX}_mrna_all_exon_site_porta.tsv \
    >output/Hsa_${PREFIX}_mrna_all_exon_site.tsv
  Rscript NJU_seq/presentation/exon_distance.R \
    output/Hsa_${PREFIX}_mrna_all_exon_site_bar.tsv \
    output/Hsa_${PREFIX}_mrna_all_exon_site_porta.tsv \
    output/Hsa_${PREFIX}_mrna_all_exon_site.pdf
  rm output/Hsa_${PREFIX}_mrna_all_exon_site_bar.tsv
  rm output/Hsa_${PREFIX}_mrna_all_exon_site_porta.tsv
done

for PREFIX in HeLa A549 293T; do
  perl NJU_seq/mrna_analysis/exon_distance_patch.pl \
    <output/Hsa_${PREFIX}_mrna_all_exon_site.tsv \
    >output/Hsa_${PREFIX}_mrna_all_exon_site_filtered.tsv
done
