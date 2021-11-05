umi_tools extract \
  --bc-pattern=NNNNNN --bc-pattern2=NNNNNN \
  -I cater_data/NJU9233/NJU9233_R1.fq.gz \
  -S fat/NJU9233/processed.1.fastq.gz \
  --read2-in=cater_data/NJU9233/NJU9233_R2.fq.gz \
  --read2-out=fat/NJU9233/processed.2.fastq.gz

time bwa-mem2 mem -t 12 \
  fat/index/hg38_genome.fa.gz \
  fat/NJU9233/processed.1.fastq.gz \
  fat/NJU9233/processed.2.fastq.gz \
  1>fat/NJU9233/bwa_raw.sam \
  2>fat/NJU9233/bwa.log

samtools sort -@ 8 \
  fat/NJU9233/bwa_raw.sam \
  -o fat/NJU9233/bwa_sorted.bam

samtools index fat/NJU9233/bwa_sorted.bam

umi_tools dedup \
  --paired --output-stats=deduplicated \
  -I fat/NJU9233/bwa_sorted.bam \
  -S fat/NJU9233/bwa_deduplicated.bam
