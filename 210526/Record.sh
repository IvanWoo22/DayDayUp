for PREFIX in NJU60{92..99} NJU61{48..51}; do
  bsub -n 24 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
cutadapt -a AGATCGGAAGAGC -A GATCGTCGGACTG \
-O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1_test1.fq.gz -p ${PREFIX}/R2_test1.fq.gz \
Cell/${PREFIX}_*R1*.gz Cell/${PREFIX}_*R2*.gz -j 20
"
done

for PREFIX in NJU60{92..99} NJU61{48..51}; do
  bsub -n 24 -o ../log/${PREFIX}_cutadapt1.log -J "${PREFIX}" "
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
-O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1.fq.gz -p ${PREFIX}/R2.fq.gz \
Cell/${PREFIX}_*R1*.gz Cell/${PREFIX}_*R2*.gz -j 20
"
done

for PREFIX in NJU60{92..99} NJU61{48..51}; do
  mv log/${PREFIX}_cutadapt1.log log/${PREFIX}_cutadapt.log
  rm data/${PREFIX}/R1_test1.fq.gz data/${PREFIX}/R2_test1.fq.gz
done

for PREFIX in NJU60{92..99} NJU61{48..51}; do
  bsub -n 24 -o log/${PREFIX}_rrna_alignment.log -J "${PREFIX}" "
time bowtie2 -p 22 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/hsa_rrna -1 data/${PREFIX}/R1.fq.gz -2 data/${PREFIX}/R2.fq.gz -S output/${PREFIX}/rrna.raw.sam
"
done

perl NJU_seq/tool/stat_alignment.pl \
  log/NJU{60{92..99},61{48..51}}_rrna_alignment.log |
  Rscript NJU_seq/tool/draw_table.R \
    output/Celline_rrna.bowtie2.pdf

for PREFIX in NJU60{92..99} NJU61{48..51}; do
  bsub -n 4 -o log/${PREFIX}_rrna_alignment.log -J "${PREFIX}" "
samtools view -h -f 97 -F 144 output/${PREFIX}/rrna.raw.sam > output/${PREFIX}/rrna.temp.sam
samtools view -f 145 -F 96 output/${PREFIX}/rrna.raw.sam >> output/${PREFIX}/rrna.temp.sam
samtools sort -n output/${PREFIX}/rrna.temp.sam | samtools view >output/${PREFIX}/rrna.filter.sam
rm output/${PREFIX}/rrna.temp.sam
pigz output/${PREFIX}/rrna.raw.sam
"
done

for PREFIX in NJU60{12..23}; do
  mkdir data/${PREFIX} output/${PREFIX} temp/${PREFIX}
  mv data/Cell/${PREFIX}_*R1*.gz data/${PREFIX}/R1.fq.gz
  mv data/Cell/${PREFIX}_*R2*.gz data/${PREFIX}/R2.fq.gz
done

for PREFIX in NJU60{12..23}; do
  bsub -n 24 -o log/${PREFIX}_rrna_alignment.log -J "${PREFIX}" "
time bowtie2 -p 22 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/hsa_rrna -1 data/${PREFIX}/R1.fq.gz -2 data/${PREFIX}/R2.fq.gz -S output/${PREFIX}/rrna.raw.sam
"
done

perl NJU_seq/tool/stat_alignment.pl \
  log/NJU60{12..23}_rrna_alignment.log |
  Rscript NJU_seq/tool/draw_table.R \
    output/Celline_rrna.bowtie2.pdf

for PREFIX in NJU60{12..23}; do
  bsub -n 4 -o log/${PREFIX}_rrna_filter.log -J "${PREFIX}" "
samtools view -h -f 97 -F 144 output/${PREFIX}/rrna.raw.sam > temp/${PREFIX}/rrna.temp.sam
samtools view -f 145 -F 96 output/${PREFIX}/rrna.raw.sam >> temp/${PREFIX}/rrna.temp.sam
samtools sort -n temp/${PREFIX}/rrna.temp.sam | samtools view >output/${PREFIX}/rrna.filter.sam
rm temp/${PREFIX}/rrna.temp.sam
pigz output/${PREFIX}/rrna.raw.sam
"
done

for PREFIX in NJU60{12..23}; do
  bsub -n 4 -o log/${PREFIX}_rrna_rev.log -J "${PREFIX}" "
samtools view -h -f 81 -F 160 output/${PREFIX}/rrna.raw.sam.gz > temp/${PREFIX}/rrna.temp.sam
samtools view -f 161 -F 80 output/${PREFIX}/rrna.raw.sam.gz >> temp/${PREFIX}/rrna.temp.sam
samtools sort -n temp/${PREFIX}/rrna.temp.sam | samtools view >output/${PREFIX}/rrna.rev.sam
rm temp/${PREFIX}/rrna.temp.sam
"
done

for PREFIX in NJU60{12..23}; do
  bsub -n 4 -o log/${PREFIX}_rrna_judge.log -J "${PREFIX}" "
time awk '\$6!=\"*\"&&\$7==\"=\"{print \$1 \"\t\" \$3 \"\t\" \$4 \"\t\" \$6 \"\t\" \$10}' output/${PREFIX}/rrna.filter.sam |
perl NJU_seq/rrna_analysis/matchquality_judge.pl |
perl NJU_seq/rrna_analysis/multimatch_judge.pl \
>temp/${PREFIX}/rrna.out.tmp
"
done

for PREFIX in NJU60{12..23}; do
  bsub -n 4 -o log/${PREFIX}_rrna_judge.log -J "${PREFIX}" "
time awk '\$6!=\"*\"&&\$7==\"=\"{print \$1 \"\t\" \$3 \"\t\" \$4 \"\t\" \$6 \"\t\" \$10}' output/${PREFIX}/rrna.rev.sam |
perl NJU_seq/rrna_analysis/matchquality_judge.pl |
perl NJU_seq/rrna_analysis/multimatch_judge.pl \
>temp/${PREFIX}/rrna.rev.tmp
"
done

for PREFIX in NJU60{12..23}; do
  perl NJU_seq/rrna_analysis/readend_count.pl \
    NJU_seq/data/hsa_rrna/28s.fa temp/${PREFIX}/rrna.out.tmp 28s \
    >output/${PREFIX}/rrna_28s.tsv
  perl NJU_seq/rrna_analysis/readend_count.pl \
    NJU_seq/data/hsa_rrna/18s.fa temp/${PREFIX}/rrna.out.tmp 18s \
    >output/${PREFIX}/rrna_18s.tsv
done

for PREFIX in NJU60{12..23}; do
  perl NJU_seq/rrna_analysis/readend_count.pl \
    NJU_seq/data/hsa_rrna/28s.fa temp/${PREFIX}/rrna.rev.tmp 28s \
    >output/${PREFIX}/rrna_28s_art.tsv
  perl NJU_seq/rrna_analysis/readend_count.pl \
    NJU_seq/data/hsa_rrna/18s.fa temp/${PREFIX}/rrna.rev.tmp 18s \
    >output/${PREFIX}/rrna_18s_art.tsv
done

time parallel -j 3 "
perl NJU_seq/rrna_analysis/score.pl \\
output/{1}/rrna_28s.tsv \\
output/{2}/rrna_28s.tsv \\
output/{3}/rrna_28s.tsv \\
output/{4}/rrna_28s.tsv \\
>output/{5}_rrna_28s_scored.tsv
perl NJU_seq/rrna_analysis/score.pl \\
output/{1}/rrna_18s.tsv \\
output/{2}/rrna_18s.tsv \\
output/{3}/rrna_18s.tsv \\
output/{4}/rrna_18s.tsv \\
>output/{5}_rrna_18s_scored.tsv
" ::: NJU60{12..23..4} ::: NJU60{13..23..4} ::: NJU60{14..23..4} ::: NJU60{15..23..4} ::: HeLa A549 HEK293T

time parallel --xapply --keep-order -j 3 "
perl NJU_seq/rrna_analysis/score_beta3.pl \\
output/{1}/rrna_28s.tsv \\
output/{2}/rrna_28s.tsv \\
output/{3}/rrna_28s.tsv \\
output/{4}/rrna_28s.tsv \\
>output/{5}_rrna_28s_scored_b3.tsv
perl NJU_seq/rrna_analysis/score_beta3.pl \\
output/{1}/rrna_18s.tsv \\
output/{2}/rrna_18s.tsv \\
output/{3}/rrna_18s.tsv \\
output/{4}/rrna_18s.tsv \\
>output/{5}_rrna_18s_scored_b3.tsv
" ::: NJU60{12..23..4} ::: NJU60{13..23..4} ::: NJU60{14..23..4} ::: NJU60{15..23..4} ::: HeLa A549 HEK293T

time parallel --xapply --keep-order -j 3 "
perl NJU_seq/rrna_analysis/score_beta4.pl \\
output/{1}/rrna_28s.tsv \\
output/{2}/rrna_28s.tsv \\
output/{3}/rrna_28s.tsv \\
output/{4}/rrna_28s.tsv \\
>output/{5}_rrna_28s_scored_b4.tsv
perl NJU_seq/rrna_analysis/score_beta4.pl \\
output/{1}/rrna_18s.tsv \\
output/{2}/rrna_18s.tsv \\
output/{3}/rrna_18s.tsv \\
output/{4}/rrna_18s.tsv \\
>output/{5}_rrna_18s_scored_b4.tsv
" ::: NJU60{12..23..4} ::: NJU60{13..23..4} ::: NJU60{14..23..4} ::: NJU60{15..23..4} ::: HeLa A549 HEK293T

time parallel --xapply --keep-order -j 3 "
perl NJU_seq/rrna_analysis/score_beta5.pl \\
output/{1}/rrna_28s.tsv \\
output/{2}/rrna_28s.tsv \\
output/{3}/rrna_28s.tsv \\
output/{4}/rrna_28s.tsv \\
>output/{5}_rrna_28s_scored_b5.tsv
perl NJU_seq/rrna_analysis/score_beta5.pl \\
output/{1}/rrna_18s.tsv \\
output/{2}/rrna_18s.tsv \\
output/{3}/rrna_18s.tsv \\
output/{4}/rrna_18s.tsv \\
>output/{5}_rrna_18s_scored_b5.tsv
" ::: NJU60{12..23..4} ::: NJU60{13..23..4} ::: NJU60{14..23..4} ::: NJU60{15..23..4} ::: HeLa A549 HEK293T
