for PREFIX in NJU7{{055..057},059,{061..068}}; do
  mkdir ${PREFIX}
  rm ../log/${PREFIX}_cutadapt.log
  bsub -n 24 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGAT \
-O 6 -m 10 -e 0.1 --discard-untrimmed -o ${PREFIX}/R1.fq.gz -p ${PREFIX}/R2.fq.gz \
RBC/${PREFIX}_*R1*.gz RBC/${PREFIX}_*R2*.gz -j 20
"
done

for PREFIX in NJU7{{055..057},059,{061..068}}; do
  mkdir -p "temp/${PREFIX}" "output/${PREFIX}"
  bsub -n 3 -o log/${PREFIX}_peco.log -J "${PREFIX}" "
perl NJU_seq/quality_control/pe_consistency.pl \
data/${PREFIX}/R1.fq.gz data/${PREFIX}/R2.fq.gz \
temp/${PREFIX}.fq.gz
"
done

for PREFIX in NJU7{{055..057},059,{061..068}}; do
  bsub -n 24 -o log/${PREFIX}_rrna_alignment.log -J "${PREFIX}" "
time bowtie2 -p 22 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 --np 0 --xeq -x index/hsa_rrna -1 data/${PREFIX}/R1.fq.gz -2 data/${PREFIX}/R2.fq.gz -S output/${PREFIX}/rrna.raw.sam
"
done

for PREFIX in NJU7{{055..057},059,{061..068}}; do
  bsub -n 4 -o log/${PREFIX}_rrna_filter.log -J "${PREFIX}" "
samtools view -h -f 97 -F 144 output/${PREFIX}/rrna.raw.sam > output/${PREFIX}/rrna.temp.sam
samtools view -f 145 -F 96 output/${PREFIX}/rrna.raw.sam >> output/${PREFIX}/rrna.temp.sam
samtools sort -n output/${PREFIX}/rrna.temp.sam | samtools view >output/${PREFIX}/rrna.filter.sam
rm output/${PREFIX}/rrna.temp.sam
pigz output/${PREFIX}/rrna.raw.sam
"
done

for PREFIX in NJU7{{055..057},059,{061..068}}; do
  bsub -n 4 -o log/${PREFIX}_rrna_judge.log -J "${PREFIX}" "
time awk '\$6!=\"*\"&&\$7==\"=\"{print \$1 \"\t\" \$3 \"\t\" \$4 \"\t\" \$6 \"\t\" \$10}' output/${PREFIX}/rrna.filter.sam |
perl NJU_seq/rrna_analysis/matchquality_judge.pl |
perl NJU_seq/rrna_analysis/multimatch_judge.pl \
>temp/${PREFIX}/rrna.out.tmp
"
done

bsub -n 24 -J "rrna2mrna" "
parallel --xapply --keep-order -j 10 '
time bash NJU_seq/tool/extract_fastq.sh \
temp/{}/rrna.out.tmp \
data/{}/R1.fq.gz data/{}/R1.mrna.fq.gz \
data/{}/R2.fq.gz data/{}/R2.mrna.fq.gz
' ::: NJU7{{055..057},059,{061..068}}
"

for PREFIX in NJU7{{055..057},059,{061..068}}; do
  rm log/${PREFIX}_mrna_alignment.log
  bsub -n 80 -q fat_384 -o log/${PREFIX}_mrna_alignment.log -J "${PREFIX}" "
time bowtie2 -p 72 -a -t --end-to-end -D 20 -R 3 -N 0 -L 10 --score-min C,0,0 --xeq -x index/hsa_basic_protein_coding -1 data/${PREFIX}/R1.mrna.fq.gz -2 data/${PREFIX}/R2.mrna.fq.gz | pigz > output/${PREFIX}/mrna.raw.sam.gz
"
done

bsub -n 24 -q largemem -J "filter" "
parallel --xapply --keep-order -j 3 '
samtools view -@ 8 -h -f 97 -F 144 output/{}/mrna.raw.sam.gz > output/{}/mrna.temp.sam
samtools view -@ 8 -f 145 -F 96 output/{}/mrna.raw.sam.gz >> output/{}/mrna.temp.sam
samtools sort -m 2G -@ 8 -n output/{}/mrna.temp.sam -o output/{}/mrna.filter.sam
rm output/{}/mrna.temp.sam
pigz output/{}/mrna.filter.sam
' ::: NJU7{{055..057},059,{061..068}}
"

for PREFIX in NJU7{{055..057},059,{061..068}}; do
  rm log/${PREFIX}_mrna_dedup.log
  bsub -n 24 -o log/${PREFIX}_mrna_dedup.log -J "${PREFIX}" "
time pigz -dcf output/${PREFIX}/mrna.filter.sam.gz |
parallel --pipe --block 1G --no-run-if-empty --linebuffer --keep-order -j 12 '
awk '\''\$6!=\"*\"&&\$7==\"=\"{print \$1 \"\t\" \$3 \"\t\" \$4 \"\t\" \$6 \"\t\" \$10}
'\'' | perl NJU_seq/mrna_analysis/multimatch_judge.pl
' | perl NJU_seq/mrna_analysis/multimatch_judge.pl |
parallel --pipe --block 1G --no-run-if-empty --linebuffer --keep-order -j 12 '
perl NJU_seq/mrna_analysis/dedup.pl --refstr \"Parent=\" --transid \"ENST\" --info data/hsa_exon.info
' | perl NJU_seq/mrna_analysis/dedup.pl --refstr \"Parent=\" --transid \"ENST\" --info data/hsa_exon.info >temp/${PREFIX}/mrna.dedup.tmp
"
done

bsub -n 24 -J "almostunique" "
parallel --xapply --keep-order -j 6 '
time bash NJU_seq/mrna_analysis/almostunique.sh \
temp/{}/mrna.dedup.tmp \
data/{}/R1.mrna.fq.gz \
temp/{} \
temp/{}/mrna.almostunique.tmp
' ::: NJU7{{055..057},059,{061..068}}
"

for PREFIX in NJU7{{055..057},059,{061..068}}; do
  time perl NJU_seq/mrna_analysis/count.pl \
    temp/"${PREFIX}"/mrna.almostunique.tmp \
    >temp/"${PREFIX}"/mrna.count.tmp
done

for PREFIX in NJU7{{055..057},059,{061..068}}; do
  time gzip -dcf data/hsa.gff3.gz |
    awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
    perl NJU_seq/mrna_analysis/merge.pl \
      --refstr "Parent=" \
      --geneid "ENSG" \
      --transid "ENST" \
      -i temp/"${PREFIX}"/mrna.count.tmp \
      -o output/"${PREFIX}"/mrna.tsv
done

parallel --keep-order -j 12 '
rm output/{}/mrna.cov
echo {} >>output/{}/mrna.cov
bash NJU_seq/presentation/seq_depth.sh \
temp/{}/mrna.almostunique.tmp \
output/{}/mrna.tsv \
>>output/{}/mrna.cov
' ::: NJU7{{055..057},059,{061..068}} NJU70{05..18} NJU7044 NJU7045 NJU7{{031..036},039,040,042,043,{047..054},073,{075..102}}

parallel -j 4 "
perl score_neo.pl \\
output/NJU7068/mrna.tsv \\
output/{}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_scored.tsv
" ::: NJU7031 NJU7032 NJU7033 NJU7034 NJU7035 NJU7036 NJU7039 NJU7040 NJU7066 NJU7067

parallel -j 4 "
perl score_neo.pl \\
output/NJU7053/mrna.tsv \\
output/{}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_scored.tsv
" ::: NJU7042 NJU7043 NJU7047 NJU7048 NJU7049 NJU7050 NJU7051 NJU7052

parallel -j 4 "
perl score_neo.pl \\
output/NJU7065/mrna.tsv \\
output/{}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_scored.tsv
" ::: NJU7054 NJU7055 NJU7056 NJU7057 NJU7059 NJU7061 NJU7062 NJU7063 NJU7064

parallel -j 4 "
perl score_neo.pl \\
output/NJU7092/mrna.tsv \\
output/{}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_scored.tsv
" ::: NJU7082 NJU7083 NJU7084 NJU7085 NJU7086 NJU7087 NJU7088 NJU7089 NJU7090 NJU7091

NJU70{{31..36},39,40,42,43,{47..52},{54..57},59,{61..64},66,67,{82..91}}

parallel -j 4 "
perl score_neo.pl \\
output/NJU7013/mrna.tsv \\
output/{}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_scored.tsv
" ::: NJU7010 NJU7011 NJU7012

parallel -j 4 "
perl score_neo.pl \\
output/NJU7081/mrna.tsv \\
output/{}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_scored.tsv
" ::: NJU7014 NJU7015 NJU7016 NJU7017 NJU7018

parallel -j 4 "
perl score_neo.pl \\
output/NJU7080/mrna.tsv \\
output/{}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_scored.tsv
" ::: NJU7073 NJU7075 NJU7076 NJU7077 NJU7078 NJU7079

parallel -j 4 "
perl score_neo.pl \\
output/NJU7045/mrna.tsv \\
output/{}/mrna.tsv |
sort -t $'\t' -nrk 12,12 \\
>output/{}/mrna_scored.tsv
" ::: NJU7044
