gatk PreprocessIntervals \
  -R ../index/hg38_chrom.fa \
  --padding 0 \
  -imr OVERLAPPING_ONLY \
  -O hg38.preprocessed.interval_list

gatk CollectReadCounts \
  -L hg38.preprocessed.interval_list \
  -R ../index/hg38_chrom.fa \
  -imr OVERLAPPING_ONLY \
  --disable-tool-default-read-filters true \
  -I NJU9221_bwa_deduplicated.bam \
  -O NJU9221_bwa_deduplicated.hdf5

gatk AnnotateIntervals \
  -L hg38.preprocessed.interval_list \
  -R ../index/hg38_chrom.fa \
  -imr OVERLAPPING_ONLY \
  -O hg38.annotated.tsv

gatk FilterIntervals \
  -L hg38.preprocessed.interval_list \
  --annotated-intervals hg38.annotated.tsv \
  --low-count-filter-count-threshold 3 \
  -I Sample01_bwa_deduplicated.hdf5 -I Sample02_bwa_deduplicated.hdf5 \
  -I Sample03_bwa_deduplicated.hdf5 -I Sample04_bwa_deduplicated.hdf5 \
  -I Sample05_bwa_deduplicated.hdf5 -I Sample06_bwa_deduplicated.hdf5 \
  -I Sample07_bwa_deduplicated.hdf5 -I Sample08_bwa_deduplicated.hdf5 \
  -I Sample09_bwa_deduplicated.hdf5 -I Sample10_bwa_deduplicated.hdf5 \
  -I Sample11_bwa_deduplicated.hdf5 -I Sample12_bwa_deduplicated.hdf5 \
  -I Sample13_bwa_deduplicated.hdf5 -I Sample14_bwa_deduplicated.hdf5 \
  -I Sample15_bwa_deduplicated.hdf5 -I Sample16_bwa_deduplicated.hdf5 \
  -I Sample17_bwa_deduplicated.hdf5 -I NJU9221_bwa_deduplicated.hdf5 \
  -I NJU9226_bwa_deduplicated.hdf5 -I NJU9233_bwa_deduplicated.hdf5 \
  -I NJU9238_bwa_deduplicated.hdf5 \
  -imr OVERLAPPING_ONLY \
  -O hg38.cohort.gc.filtered.interval_list

gatk DetermineGermlineContigPloidy \
  -L hg38.cohort.gc.filtered.interval_list \
  --interval-merging-rule OVERLAPPING_ONLY \
  -I Sample01_bwa_deduplicated.hdf5 -I Sample02_bwa_deduplicated.hdf5 \
  -I Sample03_bwa_deduplicated.hdf5 -I Sample04_bwa_deduplicated.hdf5 \
  -I Sample05_bwa_deduplicated.hdf5 -I Sample06_bwa_deduplicated.hdf5 \
  -I Sample07_bwa_deduplicated.hdf5 -I Sample08_bwa_deduplicated.hdf5 \
  -I Sample09_bwa_deduplicated.hdf5 -I Sample10_bwa_deduplicated.hdf5 \
  -I Sample11_bwa_deduplicated.hdf5 -I Sample12_bwa_deduplicated.hdf5 \
  -I Sample13_bwa_deduplicated.hdf5 -I Sample14_bwa_deduplicated.hdf5 \
  -I Sample15_bwa_deduplicated.hdf5 -I Sample16_bwa_deduplicated.hdf5 \
  -I Sample17_bwa_deduplicated.hdf5 -I NJU9221_bwa_deduplicated.hdf5 \
  -I NJU9226_bwa_deduplicated.hdf5 -I NJU9233_bwa_deduplicated.hdf5 \
  -I NJU9238_bwa_deduplicated.hdf5 \
  --contig-ploidy-priors hg38_contig_ploidy_priors.tsv \
  --output . \
  --output-prefix ploidy \
  --verbosity DEBUG

~/miniconda3/envs/gatk/bin/gatk PostprocessGermlineCNVCalls \
  --model-shard-path cohort/01-model \
  --model-shard-path cohort/02-model \
  --model-shard-path cohort/03-model \
  --model-shard-path cohort/04-model \
  --model-shard-path cohort/05-model \
  --model-shard-path cohort/06-model \
  --model-shard-path cohort/07-model \
  --model-shard-path cohort/08-model \
  --model-shard-path cohort/09-model \
  --model-shard-path cohort/10-model \
  --model-shard-path cohort/11-model \
  --model-shard-path cohort/12-model \
  --model-shard-path cohort/13-model \
  --model-shard-path cohort/14-model \
  --model-shard-path cohort/15-model \
  --model-shard-path cohort/16-model \
  --model-shard-path cohort/17-model \
  --model-shard-path cohort/18-model \
  --model-shard-path cohort/19-model \
  --model-shard-path cohort/20-model \
  --model-shard-path cohort/21-model \
  --model-shard-path cohort/22-model \
  --model-shard-path cohort/23-model \
  --model-shard-path cohort/24-model \
  --model-shard-path cohort/25-model \
  --model-shard-path cohort/26-model \
  --model-shard-path cohort/27-model \
  --model-shard-path cohort/28-model \
  --model-shard-path cohort/29-model \
  --calls-shard-path cohort/01-calls \
  --calls-shard-path cohort/02-calls \
  --calls-shard-path cohort/03-calls \
  --calls-shard-path cohort/04-calls \
  --calls-shard-path cohort/05-calls \
  --calls-shard-path cohort/06-calls \
  --calls-shard-path cohort/07-calls \
  --calls-shard-path cohort/08-calls \
  --calls-shard-path cohort/09-calls \
  --calls-shard-path cohort/10-calls \
  --calls-shard-path cohort/11-calls \
  --calls-shard-path cohort/12-calls \
  --calls-shard-path cohort/13-calls \
  --calls-shard-path cohort/14-calls \
  --calls-shard-path cohort/15-calls \
  --calls-shard-path cohort/16-calls \
  --calls-shard-path cohort/17-calls \
  --calls-shard-path cohort/18-calls \
  --calls-shard-path cohort/19-calls \
  --calls-shard-path cohort/20-calls \
  --calls-shard-path cohort/21-calls \
  --calls-shard-path cohort/22-calls \
  --calls-shard-path cohort/23-calls \
  --calls-shard-path cohort/24-calls \
  --calls-shard-path cohort/25-calls \
  --calls-shard-path cohort/26-calls \
  --calls-shard-path cohort/27-calls \
  --calls-shard-path cohort/28-calls \
  --calls-shard-path cohort/29-calls \
  --allosomal-contig chrX --allosomal-contig chrY \
  --contig-ploidy-calls ploidy-calls \
  --sample-index 0 \
  --output-genotyped-intervals genotyped-intervals-cohort24-twelve-NA19017.vcf.gz \
  --output-genotyped-segments genotyped-segments-cohort24-twelve-NA19017.vcf.gz \
  --output-denoised-copy-ratios genotyped-ratios-cohort24-twelve-NA19017.tsv \
  -R ../index/hg38_chrom.fa
