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

