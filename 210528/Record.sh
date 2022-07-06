for PREFIX in NJU70{05..18} NJU7044 NJU7045; do
  mkdir ${PREFIX}
  bsub -n 24 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
-O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1.fq.gz -p ${PREFIX}/R2.fq.gz \
RBC/${PREFIX}_*R1*.gz RBC/${PREFIX}_*R2*.gz -j 20
"
done

for PREFIX in NJU70{05..18} NJU7044 NJU7045; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  bsub -n 4 -o log/${PREFIX}_peco.log -J "${PREFIX}" "
perl NJU_seq/quality_control/pe_consistency.pl \
data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
temp/${PREFIX}.fq.gz
"
done

for PREFIX in NJU70{05..18} NJU7044 NJU7045; do
  bsub -n 24 -o log/${PREFIX}_rrna_alignment.log -J "${PREFIX}" "
time bowtie2 -p 22 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/hsa_rrna -1 data/${PREFIX}/R1.fq.gz -2 data/${PREFIX}/R2.fq.gz -S output/${PREFIX}/rrna.raw.sam
"
done

perl NJU_seq/tool/stat_alignment.pl \
  log/NJU70{05..09}_rrna_alignment.log |
  Rscript NJU_seq/tool/draw_table.R \
    output/NJU7005_rrna.bowtie2.pdf

perl NJU_seq/tool/stat_alignment.pl \
  log/NJU70{10..18}_rrna_alignment.log |
  Rscript NJU_seq/tool/draw_table.R \
    output/NJU7010_rrna.bowtie2.pdf

perl NJU_seq/tool/stat_alignment.pl \
  log/NJU70{44,45}_rrna_alignment.log |
  Rscript NJU_seq/tool/draw_table.R \
    output/NJU7044_rrna.bowtie2.pdf

for PREFIX in NJU70{05..18} NJU7044 NJU7045; do
  bsub -n 4 -o log/${PREFIX}_rrna_filter.log -J "${PREFIX}" "
samtools view -h -f 97 -F 144 output/${PREFIX}/rrna.raw.sam > output/${PREFIX}/rrna.temp.sam
samtools view -f 145 -F 96 output/${PREFIX}/rrna.raw.sam >> output/${PREFIX}/rrna.temp.sam
samtools sort -n output/${PREFIX}/rrna.temp.sam | samtools view >output/${PREFIX}/rrna.filter.sam
rm output/${PREFIX}/rrna.temp.sam
pigz output/${PREFIX}/rrna.raw.sam
"
done

for PREFIX in NJU70{05..18} NJU7044 NJU7045; do
  bsub -n 4 -o log/${PREFIX}_rrna_judge.log -J "${PREFIX}" "
time awk '\$6!=\"*\"&&\$7==\"=\"{print \$1 \"\t\" \$3 \"\t\" \$4 \"\t\" \$6 \"\t\" \$10}' output/${PREFIX}/rrna.filter.sam |
perl NJU_seq/rrna_analysis/matchquality_judge.pl |
perl NJU_seq/rrna_analysis/multimatch_judge.pl \
>temp/${PREFIX}/rrna.out.tmp
"
done

for PREFIX in NJU70{05..18} NJU7044 NJU7045; do
  bsub -n 2 -J "${PREFIX}" "
time bash NJU_seq/tool/extract_fastq.sh \
temp/${PREFIX}/rrna.out.tmp \
data/${PREFIX}/R1.fq.gz data/${PREFIX}/R1.mrna.fq.gz \
data/${PREFIX}/R2.fq.gz data/${PREFIX}/R2.mrna.fq.gz
"
done

for PREFIX in NJU70{05..18} NJU7044 NJU7045; do
  rm log/${PREFIX}_mrna_alignment.log
  bsub -n 24 -q largemem -o log/${PREFIX}_mrna_alignment.log -J "${PREFIX}" "
time bowtie2 -p 22 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 --score-min C,0,0 --xeq -x index/hsa_basic_protein_coding -1 data/${PREFIX}/R1.mrna.fq.gz -2 data/${PREFIX}/R2.mrna.fq.gz | pigz > output/${PREFIX}/mrna.raw.sam.gz
"
done

for PREFIX in NJU70{05..18} NJU7044 NJU7045; do
  rm log/${PREFIX}_mrna_filter.log
  bsub -n 4 -o log/${PREFIX}_mrna_filter.log -J "${PREFIX}" "
samtools view -@ 4 -h -f 97 -F 144 output/${PREFIX}/mrna.raw.sam.gz > output/${PREFIX}/mrna.temp.sam
samtools view -@ 4 -f 145 -F 96 output/${PREFIX}/mrna.raw.sam.gz >> output/${PREFIX}/mrna.temp.sam
samtools sort -m 2G -@ 4 -n output/${PREFIX}/mrna.temp.sam -o output/${PREFIX}/mrna.filter.sam
rm output/${PREFIX}/mrna.temp.sam
pigz output/${PREFIX}/mrna.filter.sam
"
done

for PREFIX in NJU70{05..18} NJU7044 NJU7045; do
  bsub -n 24 -o log/${PREFIX}_mrna_dedup.log -J "${PREFIX}" "
time pigz -dcf output/${PREFIX}/mrna.filter.sam.gz |
parallel --pipe --block 1G --no-run-if-empty --linebuffer --keep-order -j 12 '
awk '\''\$6!=\"*\"&&\$7==\"=\"{print \$1 \"\t\" \$3 \"\t\" \$4 \"\t\" \$6 \"\t\" \$10}
'\'' | perl NJU_seq/mrna_analysis/multimatch_judge.pl
' | perl NJU_seq/mrna_analysis/multimatch_judge.pl |
parallel --pipe --block 1G --no-run-if-empty --linebuffer --keep-order -j 12 '
perl NJU_seq/mrna_analysis/dedup.pl --refstr \"Parent=transcript:\" --transid \"ENST\" --info data/hsa_exon.info
' | perl NJU_seq/mrna_analysis/dedup.pl --refstr \"Parent=transcript:\" --transid \"ENST\" --info data/hsa_exon.info >temp/${PREFIX}/mrna.dedup.tmp
"
done

bsub -n 24 -J "almostunique" "
parallel --xapply --keep-order -j 4 '
time bash NJU_seq/mrna_analysis/almostunique.sh \
temp/{}/mrna.dedup.tmp \
data/{}/R1.mrna.fq.gz \
temp/{} \
temp/{}/mrna.almostunique.tmp
' ::: NJU70{05..18} NJU7044 NJU7045
"

for PREFIX in NJU70{05..18} NJU7044 NJU7045; do
  time perl NJU_seq/mrna_analysis/count.pl \
    temp/"${PREFIX}"/mrna.almostunique.tmp \
    >temp/"${PREFIX}"/mrna.count.tmp
done

for PREFIX in NJU70{05..18} NJU7044 NJU7045; do
  time gzip -dcf data/hsa.gff3.gz |
    awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
    perl NJU_seq/mrna_analysis/merge.pl \
      --refstr "Parent=" \
      --geneid "ENSG" \
      --transid "ENST" \
      -i temp/"${PREFIX}"/mrna.count.tmp \
      -o output/"${PREFIX}"/mrna.tsv
done
