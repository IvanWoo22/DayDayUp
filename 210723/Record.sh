for i in NJU7105 NJU7106; do
  ../TrimGalore-0.6.6/trim_galore --cores 6 \
    --length 36 --paired --gzip \
    ${i}/${i}_*_R1_*.fastq.gz \
    ${i}/${i}_*_R2_*.fastq.gz \
    -o ${i}/trim_galoredata
done

for i in NJU7105 NJU7106; do
  umi_tools extract --extract-method=regex \
    -I data/${i}/trim_galoredata/${i}_*_R2_*_2.fq.gz \
    -p "^(?P<umi_1>.{6}).*" \
    -S data/${i}/processed.2.fq \
    --read2-in=data/${i}/trim_galoredata/${i}_*_R1_*_1.fq.gz \
    --bc-pattern2=".*(?P<discard_1>.{6})$" \
    --read2-out=data/${i}/processed.1.fq
done

for i in NJU7105 NJU7106; do
  time bowtie2 -p 12 -a -t --end-to-end \
    -D 20 -R 3 -N 0 -L 10 -i S,1,0.50 \
    --np 0 --xeq -x index/hsa_rrna \
    -1 data/${i}/processed.1.fq \
    -2 data/${i}/processed.2.fq \
    -S output/${i}_rrna.raw_umi.sam
done

for PREFIX in NJU7105 NJU7106; do
  rm log/${PREFIX}_mrna_alignment.log
  bsub -n 80 -q fat_768 -o log/${PREFIX}_mrna_alignment.log -J "${PREFIX}" "
time bowtie2 -p 20 -a -t --end-to-end -D 15 -R 2 -N 0 -L 22 -i S,1,0.50 --score-min L,-0.6,-0.2 --xeq -x index/hsa_basic_protein_coding -1 data/${PREFIX}/processed.1.fq.gz -2 data/${PREFIX}/processed.2.fq.gz | pigz > output/${PREFIX}/mrna.raw.sam.gz
"
done

for PREFIX in NJU7105 NJU7106; do
  bsub -n 11 -o log/${PREFIX}_mrna_filter.log -J "${PREFIX}" "
samtools view -@ 8 -h -f 97 -F 144 output/${PREFIX}/mrna.raw.sam.gz > output/${PREFIX}/mrna.temp.sam
samtools view -@ 8 -f 145 -F 96 output/${PREFIX}/mrna.raw.sam.gz >> output/${PREFIX}/mrna.temp.sam
samtools sort -m 2G -@ 8 -n output/${PREFIX}/mrna.temp.sam -o output/${PREFIX}/mrna.filter.sam
rm output/${PREFIX}/mrna.temp.sam
pigz output/${PREFIX}/mrna.filter.sam
"
done
