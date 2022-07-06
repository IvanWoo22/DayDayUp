parallel --keep-order --xapply -j 4 '
mkdir {2}
bsub -n 24 -o ../log/{2}_cutadapt.log -J "{2}" "
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
-O 6 -m 10 -e 0.1 --discard-untrimmed -o {2}/R1.fq.gz -p {2}/R2.fq.gz \
{1}_*R1*.gz {1}_*R2*.gz -j 20
"
' ::: NJU{45..48} ::: NJU{6042..6045}

parallel --keep-order --xapply -j 4 '
pigz -dc {1}_*R1*.gz | pigz >{2}.R1.fq.gz
pigz -dc {1}_*R2*.gz | pigz >{2}.R2.fq.gz
' ::: NJU{69..72} ::: NJU{6066..6069}

parallel --keep-order --xapply -j 4 '
mkdir {}
bsub -n 24 -o ../log/{}_cutadapt.log -J "{}" "
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
-O 6 -m 10 -e 0.1 --discard-untrimmed -o {}/R1.fq.gz -p {}/R2.fq.gz \
{}.R1.fq.gz {}.R2.fq.gz -j 20
"
' ::: NJU{6066..6069}

for PREFIX in NJU{6066..6069}; do
  mkdir output/${PREFIX} temp/${PREFIX}
  bsub -n 24 -o log/${PREFIX}_rrna_alignment.log -J "${PREFIX}" "bash NJU_seq/log/Ath_scripts/step1.sh ${PREFIX}"
done

bsub -n 24 -q largemem -J "count" '
parallel --keep-order --xapply -j 10 '\''
time pigz -dcf output/{}/rrna.raw.sam.gz |
awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\''\'\'''\'' |
perl NJU_seq/rrna_analysis/matchquality_judge.pl |
perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' ::: NJU{6042..6045} NJU{6066..6069}
'

for RNA in 25s 18s 5-8s; do
  parallel --keep-order --xapply -j 10 "
perl NJU_seq/rrna_analysis/readend_count.pl NJU_seq/data/ath_rrna/${RNA}.fa temp/{}/rrna.out.tmp ${RNA} >output/{}/rrna_${RNA}.tsv
" ::: NJU{6042..6045} NJU{6066..6069}
done

for RNA in 25s 18s 5-8s; do
  parallel --keep-order --xapply -j 2 "
  perl NJU_seq/rrna_analysis/score.pl \\
    output/NJU{1}/rrna_${RNA}.tsv \\
    output/NJU{2}/rrna_${RNA}.tsv \\
    output/NJU{3}/rrna_${RNA}.tsv \\
    output/NJU{4}/rrna_${RNA}.tsv \\
      >output/{5}_rrna_${RNA}_scored.tsv
" ::: 6042 6066 ::: 6043 6067 ::: 6044 6068 ::: 6045 6069 ::: Ath_stem_former Ath_flower_former
done

for PREFIX in NJU{6042..6045} NJU{6066..6069}; do
  rm log/${PREFIX}_mrna_alignment.log
  bsub -q largemem -n 24 -o log/${PREFIX}_mrna_alignment.log -J "${PREFIX}" "bash NJU_seq/log/Ath_scripts/step2.sh ${PREFIX}"
done

for PREFIX in NJU{6042..6045} NJU{6066..6069}; do
  bsub -q largemem -n 24 -o log/${PREFIX}_mrna_dedup.log -J "${PREFIX}" "bash NJU_seq/log/Ath_scripts/step3.sh ${PREFIX} 16"
done

bsub -n 24 -J "almostunique" '
parallel --keep-order -j 8 '\''
time pigz -dc temp/{1}/mrna.dedup.tmp.gz > temp/{1}/mrna.dedup.tmp
time bash NJU_seq/mrna_analysis/almostunique.sh \
temp/{1}/mrna.dedup.tmp \
data/{1}/R1.mrna.fq.gz \
temp/{1} \
temp/{1}/mrna.almostunique.tmp
rm temp/{1}/mrna.dedup.tmp
'\'' ::: NJU{6042..6045} NJU{6066..6069}
'

for PREFIX in NJU{6042..6045} NJU{6066..6069}; do
  bsub -q fat_384 -n 80 -o log/${PREFIX}_bac_alignment.log -J "${PREFIX}" "bash bac_target.sh ${PREFIX}"
done

bsub -n 24 -q largemem -J "count" '
parallel --keep-order --xapply -j 8 '\''
bash NJU_seq/tool/extract_fastq_1.sh temp/{}/mrna.almostunique.tmp data/{}/R1.mrna.fq.gz data/{}/R1.bac.fq.gz data/{}/R2.mrna.fq.gz data/{}/R2.bac.fq.gz
'\'' ::: NJU{6042..6045} NJU{6066..6069}
'

for PREFIX in NJU{6042..6045} NJU{6066..6069}; do
  rm log/${PREFIX}_bac_in_mrna_alignment.log
  bsub -q largemem -n 24 -o log/${PREFIX}_bac_in_mrna_alignment.log -J "${PREFIX}" "bash bac_in_mrna.sh ${PREFIX}"
done
