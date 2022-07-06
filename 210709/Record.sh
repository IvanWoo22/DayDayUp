parallel --keep-order --xapply -j 1 '
cp output/{1}/mrna_scored.tsv output/RBC/{2}_mrna_scored.tsv
' ::: NJU7040 NJU7043 NJU7049 NJU7056 NJU7057 NJU7064 NJU7082 NJU7083 NJU7085 NJU7087 NJU7034 NJU7031 NJU7032 NJU7033 NJU7035 NJU7036 NJU7039 NJU7067 NJU7047 NJU7050 NJU7054 NJU7059 NJU7061 NJU7062 NJU7063 NJU7086 NJU7088 NJU7089 NJU7090 NJU7084 NJU7066 NJU7091 NJU7042 NJU7051 NJU7052 NJU7048 NJU7055 ::: Cancer01 Cancer02 Cancer03 Cancer04 Cancer05 Cancer06 Cancer07 Cancer08 Cancer09 Cancer10 Cancer11 Cancer12 Cancer13 Cancer14 Cancer15 Cancer16 Cancer17 Cancer18 Cancer19 Cancer20 Cancer21 Cancer22 Cancer23 Cancer24 Cancer25 Cancer26 Cancer27 Cancer28 Cancer29 Cancer30 Normal01 Normal02 Normal03 Normal04 Normal05 NA01 NA02

for PREFIX in NJU6402 NJU6403; do
  mkdir ${PREFIX}
  time cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
    -O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1.fq.gz -p ${PREFIX}/R2.fq.gz \
    ${PREFIX}*R1*.gz ${PREFIX}*R2*.gz -j 12 | tee ../log/${PREFIX}_cutadapt.log
done

for PREFIX in NJU6402 NJU6403; do
  mkdir -p data/${PREFIX} temp/${PREFIX} output/${PREFIX}
  time bowtie2 -p 16 -a -t \
    --end-to-end -D 20 -R 3 \
    -N 0 -L 10 -i S,1,0.50 --np 0 \
    --xeq -x index/ath_rrna \
    -1 data/${PREFIX}/R1.fq.gz -2 data/${PREFIX}/R2.fq.gz \
    -S output/${PREFIX}/rrna.raw.sam \
    2>&1 |
    tee log/${PREFIX}_rrna_alignment.log
done

for PREFIX in NJU6402 NJU6403; do
  time cat output/"${PREFIX}"/rrna.raw.sam |
    awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}
' |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl \
      >temp/"${PREFIX}"/rrna.out.tmp
done

for PREFIX in NJU6402 NJU6403; do
  samtools view -h -f 97 -F 144 output/${PREFIX}/rrna.raw.sam >output/${PREFIX}/rrna.temp.sam
  samtools view -f 145 -F 96 output/${PREFIX}/rrna.raw.sam >>output/${PREFIX}/rrna.temp.sam
  samtools sort -n output/${PREFIX}/rrna.temp.sam | samtools view >output/${PREFIX}/rrna.filter.sam
  rm output/${PREFIX}/rrna.temp.sam
done

for PREFIX in NJU6402 NJU6403; do
  time bash NJU_seq/tool/extract_fastq.sh \
    temp/"${PREFIX}"/rrna.out.tmp \
    data/"${PREFIX}"/R1.fq.gz data/"${PREFIX}"/R1.mrna.fq.gz \
    data/"${PREFIX}"/R2.fq.gz data/"${PREFIX}"/R2.mrna.fq.gz
done

for PREFIX in NJU6402 NJU6403; do
  time bowtie2 -p 16 -a -t \
    --end-to-end -D 20 -R 3 \
    -N 0 -L 10 --score-min C,0,0 \
    --xeq -x index/ath_protein_coding \
    -1 data/${PREFIX}/R1.mrna.fq.gz -2 data/${PREFIX}/R2.mrna.fq.gz \
    -S output/${PREFIX}/mrna.raw.sam \
    2>&1 |
    tee log/${PREFIX}_mrna_alignment.log
done

for PREFIX in NJU6402 NJU6403; do
  time awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' output/"${PREFIX}"/mrna.raw.sam |
    perl NJU_seq/mrna_analysis/multimatch_judge.pl \
      >temp/${PREFIX}/mrna.out.tmp
done

for PREFIX in NJU6402 NJU6403; do
  cat temp/${PREFIX}/mrna.out.tmp |
    perl NJU_seq/mrna_analysis/dedup.pl --refstr "Parent=transcript:" --transid "AT" --info data/ath_exon.info >temp/${PREFIX}/mrna.dedup.tmp
done

for PREFIX in NJU6402 NJU6403; do
  time bowtie2 -p 20 -a -t \
    --end-to-end -D 20 -R 3 \
    -N 0 -L 10 --score-min C,0,0 \
    --xeq -x index/bac_target \
    -1 data/${PREFIX}/R1.vrna.fq.gz -2 data/${PREFIX}/R2.vrna.fq.gz \
    -S output/${PREFIX}/bac.raw.sam \
    2>&1 |
    tee log/${PREFIX}_vrna_alignment.log
done
