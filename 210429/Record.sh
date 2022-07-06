for PREFIX in NJU{6356..6363}; do
  mkdir ${PREFIX}
  bsub -n 24 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
-O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1.fq.gz -p ${PREFIX}/R2.fq.gz \
Temp/${PREFIX}_*_R1_*.gz Temp/${PREFIX}_*_R2_*.gz -j 20
"
done

for PREFIX in NJU{6356..6363}; do
  mkdir output/${PREFIX} temp/${PREFIX}
  bsub -n 24 -o log/${PREFIX}_rrna_alignment.log -J "${PREFIX}" "bash mmu_rrna_alignment.sh ${PREFIX}"
done

bsub -n 24 -J "count" '
parallel --keep-order --xapply -j 6 '\''
time pigz -dcf output/{}/rrna.raw.sam.gz |
awk '\''\'\'''\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\''\'\'''\'' |
perl NJU_seq/rrna_analysis/matchquality_judge.pl |
perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/{}/rrna.out.tmp
time bash NJU_seq/tool/extract_fastq.sh temp/{}/rrna.out.tmp data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
'\'' ::: NJU{6356..6363}
'

for PREFIX in NJU{6356..6363}; do
  bsub -n 24 -o log/${PREFIX}_mcmv_alignment.log -J "${PREFIX}" "bash mmu_mcmv_alignment.sh ${PREFIX}"
done

for PREFIX in NJU{6356..6363}; do
  bsub -n 2 -o log/${PREFIX}_mcmv_filter.log -J "${PREFIX}" "bash mmu_mcmv_multimatch.sh ${PREFIX}"
done

for PREFIX in NJU{6356..6363}; do
  perl NJU_seq/rrna_analysis/readend_count.pl \
    data/MCMV.fa temp/${PREFIX}/MCMV.out.tmp "NC_004065.1" \
    >output/${PREFIX}/MCMV.tsv
done

time parallel --keep-order --xapply -j 2 '
perl NJU_seq/rrna_analysis/score4virus_beta1.pl output/{1}/MCMV.tsv output/{2}/MCMV.tsv output/{3}/MCMV.tsv output/{4}/MCMV.tsv >output/{5}_scored.tsv
' ::: NJU{6356..6363..4} ::: NJU{6357..6363..4} ::: NJU{6358..6363..4} ::: NJU{6359..6363..4} ::: MCMV_inf MCMV_uni

for PREFIX in NJU6{390..401}; do
  samtools view -h -f 81 -F 160 output/${PREFIX}/mhv.new_raw.sam.gz >output/${PREFIX}/mhv.new_filter.sam
  samtools view -f 161 -F 80 output/${PREFIX}/mhv.new_raw.sam.gz >>output/${PREFIX}/mhv.new_filter.sam
  samtools sort -n output/${PREFIX}/mhv.new_filter.sam | samtools view >output/${PREFIX}/mhv.new_filter_sorted.sam
  awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' output/${PREFIX}/mhv.new_filter_sorted.sam |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/${PREFIX}/mhv_rev.out.tmp
  perl NJU_seq/rrna_analysis/readend_count.pl \
    data/MHV.fa temp/${PREFIX}/mhv_rev.out.tmp "NC_048217.1" \
    >output/${PREFIX}/mhv_rev.tsv
done

cp NJU_seq/rrna_analysis/score4virus_beta1.pl score4virus_beta2.pl

time parallel --keep-order --xapply -j 3 '
perl score4virus_beta2.pl output/{1}/mhv_rev.tsv output/{2}/mhv_rev.tsv output/{3}/mhv_rev.tsv output/{4}/mhv_rev.tsv >output/{5}_scored.tsv
' ::: NJU6{390..401..4} ::: NJU6{391..401..4} ::: NJU6{392..401..4} ::: NJU6{393..401..4} ::: MHV_rev_par MHV_rev_inf MHV_rev_uni

for PREFIX in NJU{6356..6363}; do
  samtools view -h -f 81 -F 160 output/${PREFIX}/MCMV.raw.sam.gz >output/${PREFIX}/mhv.new_filter.sam
  samtools view -f 161 -F 80 output/${PREFIX}/MCMV.raw.sam.gz >>output/${PREFIX}/mhv.new_filter.sam
  samtools sort -n output/${PREFIX}/mhv.new_filter.sam | samtools view >output/${PREFIX}/mhv.new_filter_sorted.sam
  awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' output/${PREFIX}/mhv.new_filter_sorted.sam |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/${PREFIX}/mhv_rev.out.tmp
  perl NJU_seq/rrna_analysis/readend_count.pl \
    data/MHV.fa temp/${PREFIX}/mhv_rev.out.tmp "NC_004065.1" \
    >output/${PREFIX}/mhv_rev.tsv
done

time parallel --keep-order --xapply -j 2 '
perl score4virus_beta2.pl output/{1}/mhv_rev.tsv output/{2}/mhv_rev.tsv output/{3}/mhv_rev.tsv output/{4}/mhv_rev.tsv >output/{5}_scored.tsv
' ::: NJU{6356..6363..4} ::: NJU{6357..6363..4} ::: NJU{6358..6363..4} ::: NJU{6359..6363..4} ::: MCMV_rev_inf MCMV_rev_uni

Ath-stress-ck-leaf-RF
Ath-stress-ck-root-RF
Ath-stress-ck-leaf
Ath-stress-ck-root
Ath-Al-20mM-leaf-RF
Ath-Al-20mM-leaf
Ath-Al-20mM-root-RF
Ath-Al-20mM-root
Ath-Al-60mM-leaf-RF
Ath-Al-60mM-leaf
Ath-Al-60mM-root-RF
Ath-Al-60mM-root
Ath-Cd-5mM-leaf-RF
Ath-Cd-5mM-leaf
Ath-Cd-5mM-root-RF
Ath-Cd-5mM-root
Ath-Cd-15mM-leaf-RF
Ath-Cd-15mM-leaf
Ath-Cd-15mM-root-RF
Ath-Cd-15mM-root

for PREFIX in NJU{6390..6401}; do
  samtools view -h -f 97 -F 144 output/${PREFIX}/mhv.new_raw.sam.gz >output/${PREFIX}/mhv.new_filter.sam
  samtools view -f 145 -F 96 output/${PREFIX}/mhv.new_raw.sam.gz >>output/${PREFIX}/mhv.new_filter.sam
  samtools sort -n output/${PREFIX}/mhv.new_filter.sam | samtools view >output/${PREFIX}/mhv.new_filter_sorted.sam
  awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' output/${PREFIX}/mhv.new_filter_sorted.sam |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/${PREFIX}/mhv_rev.out.tmp
  perl NJU_seq/rrna_analysis/readend_count.pl \
    data/MHV.fa temp/${PREFIX}/mhv_rev.out.tmp "NC_048217.1" \
    >output/${PREFIX}/mhv_rev.tsv
done

for PREFIX in NJU{6390..6401}; do
  awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' output/${PREFIX}/mhv.new_filter_sorted.sam |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/${PREFIX}/mhv_filter.out.tmp
  perl NJU_seq/rrna_analysis/readend_count.pl \
    data/MHV.fa temp/${PREFIX}/mhv_filter.out.tmp "NC_048217.1" \
    >output/${PREFIX}/mhv_filter.tsv
done

time parallel --keep-order --xapply -j 3 '
perl score4virus_beta3.pl output/{1}/mhv_filter.tsv output/{2}/mhv_filter.tsv output/{3}/mhv_filter.tsv output/{4}/mhv_filter.tsv >output/{5}_scored_v4.tsv
' ::: NJU{6390..6401..4} ::: NJU{6391..6401..4} ::: NJU{6392..6401..4} ::: NJU{6393..6401..4} ::: MHV_particle MHV_infect MHV_uninfect

for PREFIX in NJU{6356..6363}; do
  samtools view -h -f 81 -F 160 output/{}/MCMV.raw.sam.gz >output/{}/mcmv.rev.sam
  samtools view -f 161 -F 80 output/{}/MCMV.raw.sam.gz >>output/{}/mcmv.rev.sam
  samtools sort -n output/{}/mcmv.rev.sam | samtools view >output/{}/mcmv.rev_sorted.sam
  awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' output/${PREFIX}/mhv.rev_sorted.sam |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/${PREFIX}/mhv_rev.out.tmp
  perl NJU_seq/rrna_analysis/readend_count.pl \
    data/MCMV.fa temp/${PREFIX}/mhv_rev.out.tmp "NC_004065.1" \
    >output/${PREFIX}/mcmv_rev.tsv
done

for PREFIX in NJU{6356..6363}; do
  samtools view -h -f 97 -F 144 output/${PREFIX}/MCMV.raw.sam.gz >output/${PREFIX}/mcmv.filter.sam
  samtools view -f 145 -F 96 output/${PREFIX}/MCMV.raw.sam.gz >>output/${PREFIX}/mcmv.filter.sam
  samtools sort -n output/${PREFIX}/mcmv.filter.sam | samtools view >output/${PREFIX}/mcmv.filter_sorted.sam
  awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' output/${PREFIX}/mcmv.filter_sorted.sam |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl >temp/${PREFIX}/mcmv_filter.out.tmp
  perl NJU_seq/rrna_analysis/readend_count.pl \
    data/MCMV.fa temp/${PREFIX}/mcmv_filter.out.tmp "NC_004065.1" \
    >output/${PREFIX}/mcmv_filter.tsv
done

time parallel --keep-order --xapply -j 2 '
perl NJU_seq/vRNA_analysis/score_beta3.pl output/{1}/mcmv_filter.tsv output/{2}/mcmv_filter.tsv output/{3}/mcmv_filter.tsv output/{4}/mcmv_filter.tsv >output/{5}_scored_v4.tsv
' ::: NJU{6356..6363..4} ::: NJU{6357..6363..4} ::: NJU{6358..6363..4} ::: NJU{6359..6363..4} ::: MCMV_inf MCMV_uni

time parallel --keep-order --xapply -j 2 '
perl score4virus_beta3.pl output/{1}/mcmv_rev.tsv output/{2}/mcmv_rev.tsv output/{3}/mcmv_rev.tsv output/{4}/mcmv_rev.tsv >output/{5}_rev_scored_v4.tsv
' ::: NJU{6356..6363..4} ::: NJU{6357..6363..4} ::: NJU{6358..6363..4} ::: NJU{6359..6363..4} ::: MCMV_inf MCMV_uni
