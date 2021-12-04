#!/usr/bin/env bash

for PREFIX in HeLa A549 293T; do
  pigz -dc data/hsa.gff3.gz | awk '$3=="gene"' | perl NJU_seq/tool/add_gene_name.pl \
    --id "gene_id=" --name "gene_name=" \
    --col "8" --file "output/${PREFIX}_mrna_Nm_tmp1.tsv" \
    >output/${PREFIX}_mrna_Nm_tmp2.tsv
  pigz -dc data/hsa_ens.gff3.gz |
    awk '$3=="gene"' | perl NJU_seq/tool/add_gene_name.pl \
    --id "gene_id=" --name "description=" \
    --col "9" --file "output/${PREFIX}_mrna_Nm_tmp2.tsv" \
    >output/${PREFIX}_mrna_Nm_tmp3.tsv
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/GRCh38.primary_assembly.genome.fa output/${PREFIX}_mrna_Nm_tmp3.tsv \
    10 10 >output/${PREFIX}_mrna_Nm_tmp4.tsv
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/GRCh38.primary_assembly.genome.fa output/${PREFIX}_mrna_Nm_tmp4.tsv \
    20 20 >output/${PREFIX}_mrna_Nm_tmp5.tsv
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/GRCh38.primary_assembly.genome.fa output/${PREFIX}_mrna_Nm_tmp5.tsv \
    50 50 >output/${PREFIX}_mrna_Nm_tmp6.tsv

  perl NJU_seq/mrna_analysis/main_transcript_3.pl \
    data/hsa_basic_transcript_region.tsv \
    output/${PREFIX}_mrna_Nm_tmp6.tsv \
    output/${PREFIX}_mrna_Nm_tmp7.tsv |
    perl NJU_seq/mrna_analysis/main_transcript_4.pl \
      >output/${PREFIX}_mrna_Nm_distribution_norm.bak
  Rscript NJU_seq/presentation/point_distribution.R \
    output/${PREFIX}_mrna_Nm_distribution_norm.bak \
    output/${PREFIX}_mrna_Nm_distribution.pdf
  perl NJU_seq/mrna_analysis/exon_distance_2.pl \
    data/hsa_basic_transcript_exon.tsv \
    output/${PREFIX}_mrna_Nm_tmp7.tsv \
    output/${PREFIX}_mrna_Nm_exon_site_bar.tsv \
    output/${PREFIX}_mrna_Nm_exon_site_porta.tsv \
    >output/${PREFIX}_mrna_Nm_exon_site.tsv
  Rscript NJU_seq/presentation/exon_distance.R \
    output/${PREFIX}_mrna_Nm_exon_site_bar.tsv \
    output/${PREFIX}_mrna_Nm_exon_site_porta.tsv \
    output/${PREFIX}_mrna_Nm_exon_site.pdf
done

pigz -dc data/hsa.gff3.gz |
  grep "protein_coding" |
  awk '$3=="start_codon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
  perl NJU_seq/mrna_analysis/codon_distance_1.pl data/hsa_basic_main_transcript.txt \
    >data/hsa_basic_main_transcript_start_codon.tsv

pigz -dc data/hsa.gff3.gz |
  grep "protein_coding" |
  awk '$3=="stop_codon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
  perl NJU_seq/mrna_analysis/codon_distance_1.pl data/hsa_basic_main_transcript.txt \
    >data/hsa_basic_main_transcript_stop_codon.tsv

for PREFIX in HeLa A549 293T; do
  perl NJU_seq/mrna_analysis/codon_distance_3.pl \
    data/hsa_basic_main_transcript_start_codon.yml \
    output/${PREFIX}_mrna_Nm_tmp7.tsv \
    output/${PREFIX}_mrna_Nm_start_codon_distance_bar.tsv \
    >output/${PREFIX}_mrna_Nm_start_codon_distance.tsv
  perl NJU_seq/mrna_analysis/codon_distance_3.pl \
    data/hsa_basic_main_transcript_stop_codon.yml \
    output/${PREFIX}_mrna_Nm_tmp7.tsv \
    output/${PREFIX}_mrna_Nm_stop_codon_distance_bar.tsv \
    >output/${PREFIX}_mrna_Nm_stop_codon_distance.tsv
done

for PREFIX in HeLa A549 293T; do
  Rscript NJU_seq/presentation/codon_distance.R \
    output/${PREFIX}_mrna_Nm_start_codon_distance_bar.tsv \
    output/${PREFIX}_mrna_Nm_stop_codon_distance_bar.tsv \
    output/${PREFIX}_mrna_Nm_codon_distance.pdf
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

for PREFIX in NJU615{{2..3},{5..7},9} NJU616{0,1}; do
  echo ${PREFIX}
  perl NJU_seq/mrna_analysis/score.pl output/NJU6162/mrna_basic.tsv output/${PREFIX}/mrna_basic.tsv >output/${PREFIX}/mrna_Nm_tmp1.tsv
done

for PREFIX in NJU616{{4,5},{7..9}} NJU617{1..3}; do
  echo ${PREFIX}
  perl NJU_seq/mrna_analysis/score.pl output/NJU6174/mrna_basic.tsv output/${PREFIX}/mrna_basic.tsv >output/${PREFIX}/mrna_Nm_tmp1.tsv
done

for PREFIX in N2A_MHV N2A_HC; do
  pigz -dc data/mmu.gff3.gz | awk '$3=="gene"' | perl NJU_seq/tool/add_gene_name.pl \
    --id "gene_id=" --name "gene_name=" \
    --col "8" --file "output/${PREFIX}_mrna_Nm_tmp1.tsv" \
    >output/${PREFIX}_mrna_Nm_tmp2.tsv
  pigz -dc data/mmu_ens.gff3.gz |
    awk '$3=="gene"' | perl NJU_seq/tool/add_gene_name.pl \
    --id "gene_id=" --name "description=" \
    --col "9" --file "output/${PREFIX}_mrna_Nm_tmp2.tsv" \
    >output/${PREFIX}_mrna_Nm_tmp3.tsv
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/GRCm38.primary_assembly.genome.fa output/${PREFIX}_mrna_Nm_tmp3.tsv \
    10 10 >output/${PREFIX}_mrna_Nm_tmp4.tsv
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/GRCm38.primary_assembly.genome.fa output/${PREFIX}_mrna_Nm_tmp4.tsv \
    20 20 >output/${PREFIX}_mrna_Nm_tmp5.tsv
  perl NJU_seq/mrna_analysis/motif_nm.pl \
    data/GRCm38.primary_assembly.genome.fa output/${PREFIX}_mrna_Nm_tmp5.tsv \
    50 50 >output/${PREFIX}_mrna_Nm_tmp6.tsv
  perl NJU_seq/mrna_analysis/main_transcript_3.pl \
    data/mmu_basic_transcript_region.tsv \
    output/${PREFIX}_mrna_Nm_tmp6.tsv \
    output/${PREFIX}_mrna_Nm_tmp7.tsv |
    perl NJU_seq/mrna_analysis/main_transcript_4.pl \
      >output/${PREFIX}_mrna_Nm_distribution_norm.bak
  Rscript NJU_seq/presentation/point_distribution.R \
    output/${PREFIX}_mrna_Nm_distribution_norm.bak \
    output/${PREFIX}_mrna_Nm_distribution.pdf
  perl NJU_seq/mrna_analysis/exon_distance_2.pl \
    data/mmu_basic_transcript_exon.tsv \
    output/${PREFIX}_mrna_Nm_tmp7.tsv \
    output/${PREFIX}_mrna_Nm_exon_site_bar.tsv \
    output/${PREFIX}_mrna_Nm_exon_site_porta.tsv \
    >output/${PREFIX}_mrna_Nm_exon_site.tsv
  Rscript NJU_seq/presentation/exon_distance.R \
    output/${PREFIX}_mrna_Nm_exon_site_bar.tsv \
    output/${PREFIX}_mrna_Nm_exon_site_porta.tsv \
    output/${PREFIX}_mrna_Nm_exon_site.pdf
done
