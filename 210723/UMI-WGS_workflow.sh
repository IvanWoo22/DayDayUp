PREFIX=$1
UMI=$2
KDB=$3
BDB=$4

if [ "${UMI}" = "NEB" ]; then
  ## NEB six-code umi.
  umi_tools extract \
    --bc-pattern=NNNNNN --bc-pattern2=NNNNNN \
    -I "${PREFIX}"/R1.fq.gz \
    -S "${PREFIX}"/processed.1.fq.gz \
    --read2-in="${PREFIX}"/R2.fq.gz \
    --read2-out="${PREFIX}"/processed.2.fq.gz
elif [ "${UMI}" = "VAZ" ]; then
  ## Vazyme five-or-six umi.
  umi_tools extract \
    --extract-method=regex \
    --bc-pattern="^(?P<umi_1>.{5})(?P<discard_1>.{2}).*" \
    --bc-pattern2="^(?P<umi_2>.{5})(?P<discard_2>.{2}).*" \
    -I "${PREFIX}"/R1.fq.gz -S "${PREFIX}"/processed.1.fq.gz \
    --read2-in="${PREFIX}"/R2.fq.gz --read2-out="${PREFIX}"/processed.2.fq.gz
fi

kraken2 --threads 16 --use-names \
  --gzip-compressed --paired \
  --report "${PREFIX}"/class.report \
  --classified-out "${PREFIX}"/class#.fq \
  --db "${KDB}" \
  "${PREFIX}"/processed.1.fq.gz \
  "${PREFIX}"/processed.2.fq.gz \
  --output "${PREFIX}"/class.tsv

time bwa-mem2 mem -t 12 \
  "${BDB}" \
  "${PREFIX}"/processed.1.fq.gz \
  "${PREFIX}"/processed.2.fq.gz \
  2>"${PREFIX}"/bwa.log | pigz >"${PREFIX}"/bwa_raw.sam.gz

samtools sort -@ 8 \
  "${PREFIX}"/bwa_raw.sam \
  -o "${PREFIX}"/bwa_sorted.bam

samtools index "${PREFIX}"/bwa_sorted.bam

umi_tools dedup \
  --paired --output-stats=deduplicated \
  -I "${PREFIX}"/bwa_sorted.bam \
  -S "${PREFIX}"/bwa_deduplicated.bam

samtools index "${PREFIX}"/bwa_deduplicated.bam
