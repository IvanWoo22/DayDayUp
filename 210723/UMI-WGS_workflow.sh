umi_tools extract \
  --bc-pattern=NNNNNN --bc-pattern2=NNNNNN \
  -I cater_data/NJU9221/NJU9221_R1.fq.gz \
  -S processed.1.fastq.gz \
  --read2-in=cater_data/NJU9221/NJU9221_R2.fq.gz \
  --read2-out=processed.2.fastq.gz

umi_tools dedup -I bwa_sorted.bam --paired -S bwa_deduplicated.bam --output-stats=deduplicated