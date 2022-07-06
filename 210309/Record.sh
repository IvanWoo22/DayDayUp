for PREFIX in NJU6{390..401}; do
  mkdir ${PREFIX}
  bsub -n 24 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
-O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1.fq.gz -p ${PREFIX}/R2.fq.gz \
mhv/${PREFIX}_*_R1_*.gz mhv/${PREFIX}_*_R2_*.gz -j 20
"
done

for PREFIX in NJU6{390..401}; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  bsub -n 4 -o log/${PREFIX}_peco.log -J "${PREFIX}" "
perl NJU_seq/quality_control/pe_consistency.pl \
data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
temp/${PREFIX}.fq.gz
"
done

for PREFIX in NJU6{390..401}; do
  mkdir output/${PREFIX} temp/${PREFIX}
  bsub -n 24 -o log/${PREFIX}_rrna_alignment.log -J "${PREFIX}" "bash mhv_step1.sh ${PREFIX}"
done

for PREFIX in NJU6{390..401}; do
  bsub -n 24 -o log/${PREFIX}_mhv_alignment.log -J "${PREFIX}" "bash mhv_step2.sh ${PREFIX}"
done

# mhv_step3.sh
PREFIX=$1

time pigz -dc output/"${PREFIX}"/mhv.raw.sam.gz |
  parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 20 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}
    '\'' |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl
  ' \
    >temp/"${PREFIX}"/mhv.out.tmp
#

for PREFIX in NJU6{390..401}; do
  bsub -n 24 -o log/${PREFIX}_mhv_filter.log -J "${PREFIX}" "bash mhv_step3.sh ${PREFIX}"
done

for PREFIX in NJU6{390..401}; do
  perl NJU_seq/rrna_analysis/readend_count.pl \
    data/MHV.fa temp/${PREFIX}/mhv.out.tmp "NC_048217.1" \
    >output/${PREFIX}/mhv.tsv
done

time parallel --keep-order --xapply -j 3 "
  perl NJU_seq/rrna_analysis/score.pl \\
    output/{1}/mhv.tsv \\
    output/{2}/mhv.tsv \\
    output/{3}/mhv.tsv \\
    output/{4}/mhv.tsv \\
      >output/{5}_scored.tsv
  " ::: NJU6{390..401..4} ::: NJU6{391..401..4} ::: NJU6{392..401..4} ::: NJU6{393..401..4} ::: MHV_par MHV_inf MHV_uni
