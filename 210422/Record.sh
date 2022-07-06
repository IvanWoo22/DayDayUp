for PREFIX in NBG00{01..12}; do
  mkdir ${PREFIX}
  bsub -n 24 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
-O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1.fq.gz -p ${PREFIX}/R2.fq.gz \
NBG/${PREFIX}_*_R1_*.gz NBG/${PREFIX}_*_R2_*.gz -j 20
"
done

for PREFIX in NBG0{01..12}; do
  bsub -n 24 -J "${PREFIX}" "bash NJU_seq/log/Gmax_scripts/step1.sh ${PREFIX}"
done

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 10 '\''
time pigz -dcf output/{}/rrna.raw.sam.gz |
awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\''\'\'''\'' |
perl NJU_seq/rrna_analysis/matchquality_judge.pl |
perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' ::: NBG0{01..12}
'

for PREFIX in NBG0{01..12}; do
  bsub -q fat_768 -n 80 -J "${PREFIX}" "bash NJU_seq/log/Gmax_scripts/step2.sh ${PREFIX}"
done

for PREFIX in NBG0{01..12}; do
  bsub -n 24 -o log/${PREFIX}_mrnadedup.log -J "${PREFIX}" "bash NJU_seq/log/Gmax_scripts/step3.sh ${PREFIX} 16"
done

bsub -n 24 -J "almostunique" '
parallel --keep-order --xapply -j 6 '\''
time pigz -dc temp/{}/mrna_basic.dedup.tmp.gz > temp/{}/mrna_basic.dedup.tmp
time bash NJU_seq/mrna_analysis/almostunique.sh \
temp/{}/mrna_basic.dedup.tmp \
data/{}/R1.mrna.fq.gz \
temp/{} \
temp/{}/mrna_basic.almostunique.tmp
rm temp/{}/mrna_basic.dedup.tmp
'\'' ::: NBG0{01..12}
'

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 12 '\''
time perl NJU_seq/mrna_analysis/count.pl \
temp/{}/mrna_basic.almostunique.tmp \
>temp/{}/mrna_basic.count.tmp
'\'' ::: NBG0{01..12}
'

bsub -n 24 -J "extract_fastq" '
parallel --keep-order --xapply -j 10 '\''
time bash NJU_seq/tool/extract_fastq.sh temp/{}/mrna_basic.almostunique.tmp data/{}/R1.mrna.fq.gz data/{}/R1.smv.fq.gz data/{}/R2.mrna.fq.gz data/{}/R2.smv.fq.gz
'\'' ::: NBG0{01..12}
'

for PREFIX in NBG0{01..12}; do
  bsub -n 24 -o log/${PREFIX}_smv_alignment.log -J "${PREFIX}" "bash smv_step.sh ${PREFIX}"
done

perl NJU_seq/tool/stat_alignment.pl \
  log/NBG0{01..12}_smv_alignment.log |
  Rscript NJU_seq/tool/draw_table.R \
    output/smv.pdf

for PREFIX in NBG0{01..12}; do
  time pigz -dc output/"${PREFIX}"/smv.raw.sam.gz |
    awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/"${PREFIX}"/smv.out.tmp
done

for PREFIX in NBG0{01..12}; do
  perl NJU_seq/rrna_analysis/readend_count.pl \
    data/SMV.fa temp/"${PREFIX}"/smv.out.tmp "KP710867.1" \
    >output/"${PREFIX}"/smv.tsv
done

time parallel --keep-order --xapply -j 3 "
  perl NJU_seq/rrna_analysis/score4virus_beta1.pl \\
    output/{1}/smv.tsv \\
    output/{2}/smv.tsv \\
    output/{3}/smv.tsv \\
    output/{4}/smv.tsv \\
      >output/{5}_scored.tsv
  " ::: NBG0{01..12..4} ::: NBG0{02..12..4} ::: NBG0{03..12..4} ::: NBG0{04..12..4} ::: NBG_group1_SMV NBG_group2_SMV NBG_group3_SMV
