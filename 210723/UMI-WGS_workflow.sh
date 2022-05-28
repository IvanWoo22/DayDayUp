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
  "${PREFIX}"/bwa_raw.sam.gz \
  -o "${PREFIX}"/bwa_sorted.bam

samtools index "${PREFIX}"/bwa_sorted.bam

umi_tools dedup \
  --paired --output-stats=deduplicated \
  -I "${PREFIX}"/bwa_sorted.bam \
  -S "${PREFIX}"/bwa_deduplicated.bam

samtools index "${PREFIX}"/bwa_deduplicated.bam

cnvkit.py access ../index/hg38_chrom.fa -o access.hg38.bed

cnvkit.py batch -m wgs -p 12 -r reference.cnn \
  "${PREFIX}"/bwa_deduplicated.bam \
  -d "${PREFIX}"/ --scatter --diagram

printf "chromosome\tstart\tend\tgene\tlog2\tdepth\tprobes\tweight\tci_lo\tci_hi\n" >NJU9238_bwa_deduplicated_filter.cns
perl -ne '@a=split;if($a[0]=~m/chr[0-9]*$/){print join("\t",@a);print"\n"}' <NJU9238_bwa_deduplicated.cns >NJU9238_bwa_deduplicated_filter.cns
vim NJU9238_bwa_deduplicated_filter.cns
cnvkit.py scatter -s NJU9238_bwa_deduplicated_filter.cns --y-min -2 -o NJU9238_cnv.pdf

cnvkit.py batch -m wgs NJU9226_bwa_deduplicated.bam --normal A549_bwa_deduplicated.bam -f index/hg38_chrom.fa --annotate refFlat_clean.txt -p 12 -d . --scatter --diagram
cnvkit.py scatter -s NJU9221_bwa_deduplicated_filter1.cns --y-min -2 -o NJU9221_cnv1.pdf
