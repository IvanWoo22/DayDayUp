# ==> Rawdata
for SAMPLE in NJU6232 NJU6248 NJUzp6; do
	mkdir "${SAMPLE}"
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R1*.f*q.gz /home/ivan/fat/rRNA_identify/"${SAMPLE}"/R1.fq.gz
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R2*.f*q.gz /home/ivan/fat/rRNA_identify/"${SAMPLE}"/R2.fq.gz
done
for SAMPLE in NJU6256 NJU6268 NJUzp14; do
	mkdir "${SAMPLE}"
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R1*.f*q.gz /home/ivan/fat/rRNA_identify/"${SAMPLE}"/R1.fq.gz
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R2*.f*q.gz /home/ivan/fat/rRNA_identify/"${SAMPLE}"/R2.fq.gz
done
for SAMPLE in NJU6252 NJU6264; do
	mkdir "${SAMPLE}"
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R1*.f*q.gz /home/ivan/fat/rRNA_identify/"${SAMPLE}"/R1.fq.gz
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R2*.f*q.gz /home/ivan/fat/rRNA_identify/"${SAMPLE}"/R2.fq.gz
done
mkdir NJUlysly
ln -sf /home/ivan/cater_data/nlr_data/njuseq/Mp_Sly_1217/Sly_NC/*_1.fq.gz /home/ivan/fat/rRNA_identify/NJUlysly/R1.fq.gz
ln -sf /home/ivan/cater_data/nlr_data/njuseq/Mp_Sly_1217/Sly_NC/*_2.fq.gz /home/ivan/fat/rRNA_identify/NJUlysly/R2.fq.gz
for SAMPLE in NJU6260 NJU6272; do
	mkdir "${SAMPLE}"
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R1*.f*q.gz /home/ivan/fat/rRNA_identify/"${SAMPLE}"/R1.fq.gz
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R2*.f*q.gz /home/ivan/fat/rRNA_identify/"${SAMPLE}"/R2.fq.gz
done
mkdir NJUlymp
ln -sf /home/ivan/cater_data/nlr_data/njuseq/Mp_Sly_1217/Mp_NC/*_1.fq.gz /home/ivan/fat/rRNA_identify/NJUlymp/R1.fq.gz
ln -sf /home/ivan/cater_data/nlr_data/njuseq/Mp_Sly_1217/Mp_NC/*_2.fq.gz /home/ivan/fat/rRNA_identify/NJUlymp/R2.fq.gz


for SAMPLE in NJU6232 NJU6248 NJUzp6 NJU6256 NJU6268 NJUzp14 NJU6252 NJU6264 NJUlysly NJU6260 NJU6272; do
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAG \
		-m 15 -e 0.1 -o "${SAMPLE}"/R1_cleans.fq.gz -p "${SAMPLE}"/R2_cleans.fq.gz \
		-j 16 "${SAMPLE}"/R1.fq.gz "${SAMPLE}"/R2.fq.gz \
		>"${SAMPLE}"/cutadapts.log 2>&1
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAG \
		-m 25 -e 0.1 -o "${SAMPLE}"/R1_cleanl.fq.gz -p "${SAMPLE}"/R2_cleanl.fq.gz \
		-j 16 "${SAMPLE}"/R1.fq.gz "${SAMPLE}"/R2.fq.gz \
		>"${SAMPLE}"/cutadaptl.log 2>&1
done

# ==> A.tha
vim ath_rrna.fa
cat ~/NJU_seq/data/ath_rrna/*s.fa >ath_rrna.origin.fa
bowtie-build ath_rrna.fa index/ath_rrna
bowtie-build ath_rrna.origin.fa index/ath_rrna_origin

for SAMPLE in NJU6232 NJU6248 NJUzp6; do
	bowtie -v 2 -a --best --strata \
		-S -x index/ath_rrna_origin \
		"${SAMPLE}"/R1_cleans.fq.gz \
		>"${SAMPLE}"/rrna_origin.sam
	bowtie -v 2 -a --best --strata \
		-S -x index/ath_rrna \
		"${SAMPLE}"/R1_cleans.fq.gz \
		>"${SAMPLE}"/rrna.sam
done

# select in 2
for SAMPLE in NJU6232 NJU6248 NJUzp6; do
	rm "${SAMPLE}"/rrna_origin.sam
	samtools view -bS "${SAMPLE}"/rrna.sam \
		| samtools sort - \
		| samtools view -h - \
		| awk '$1 ~ /^@/ || ($5 >= 20 && $6 !~ /S/)' \
		| samtools view -b -F 16 - \
		| samtools sort -o "${SAMPLE}"/rrna.fwd.clean.bam
	samtools index "${SAMPLE}"/rrna.fwd.clean.bam
done

for SAMPLE in NJU6232 NJU6248 NJUzp6; do
	bcftools mpileup \
		-f ath_rrna.fa \
		-Q 20 -d 100000 \
		"${SAMPLE}"/rrna.fwd.clean.bam \
		| bcftools call -mv -Oz -o "${SAMPLE}"/rRNA.raw.vcf.gz
	bcftools norm -f ath_rrna.fa \
		-m -any "${SAMPLE}"/rRNA.raw.vcf.gz \
		-Oz -o "${SAMPLE}"/rRNA.norm.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.norm.vcf.gz
	bcftools filter \
		-i 'INFO/DP>=200 && (INFO/DP4[2] + INFO/DP4[3]) / (INFO/DP4[0] + INFO/DP4[1] + INFO/DP4[2] + INFO/DP4[3]) > 0.7 && INFO/DP4[2] > 150' \
		"${SAMPLE}"/rRNA.norm.vcf.gz -Oz -o "${SAMPLE}"/rRNA.filt.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.filt.vcf.gz
done
bcftools isec \
	-n+2 -w1 \
	NJU6232/rRNA.filt.vcf.gz \
	NJU6248/rRNA.filt.vcf.gz \
	NJUzp6/rRNA.filt.vcf.gz \
	-Oz -o ath_rRNA.isec.vcf.gz
tabix -p vcf ath_rRNA.isec.vcf.gz

bowtie2-build ath_rrna.fa index/ath_rrna
for SAMPLE in NJU6232 NJU6248 NJUzp6; do
	bowtie2 --end-to-end \
		--rdg 5,3 --rfg 5,3 \
		--score-min L,-3,-2 \
		-N 1 -L 20 \
		-x index/ath_rrna \
		-U "${SAMPLE}"/R1_cleanl.fq.gz \
		| samtools view -b -F 16 -q 20 \
		| samtools sort -o "${SAMPLE}"/rrna.indel.bam
done
for SAMPLE in NJU6232 NJU6248 NJUzp6; do
	bcftools mpileup \
		-f ath_rrna.fa \
		-Q 20 -d 100000 \
		"${SAMPLE}"/rrna.indel.bam \
		| bcftools call -mv -Oz -o "${SAMPLE}"/rRNA.indel.vcf.gz
	bcftools norm -f ath_rrna.fa \
		-m -any "${SAMPLE}"/rRNA.indel.vcf.gz \
		-Oz -o "${SAMPLE}"/rRNA.indel.norm.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.indel.norm.vcf.gz
	bcftools view \
		--types indels \
		"${SAMPLE}"/rRNA.indel.norm.vcf.gz \
		| bcftools filter \
			-i 'INFO/DP>=100 && INFO/IDV>=50 && INFO/IMF>=0.5' \
			-Oz -o "${SAMPLE}"/rRNA.indel.filt.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.indel.filt.vcf.gz
done
bcftools isec \
	-n+2 -w1 \
	NJU6232/rRNA.indel.filt.vcf.gz \
	NJU6248/rRNA.indel.filt.vcf.gz \
	NJUzp6/rRNA.indel.filt.vcf.gz \
	-Oz -o ath_rRNA.indel.isec.vcf.gz
tabix -p vcf ath_rRNA.indel.isec.vcf.gz

bcftools consensus -f ath_rrna.fa ath_rRNA.isec.vcf.gz ath_rRNA.indel.isec.vcf.gz >ath_rrna.updated.fa
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP4]\n' rRNA.isec.vcf.gz >rRNA.snp.tsv
#perl update_rrna_ref.pl ath_rrna.fa rRNA.snp.tsv >ath_rrna.updated.fa

bowtie2-build ath_rrna.updated.fa index/ath_rrna_updated
for SAMPLE in NJU6232 NJU6248 NJUzp6; do
	bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x index/ath_rrna \
		-1 "${SAMPLE}"/R1_cleans.fq.gz \
		-2 "${SAMPLE}"/R2_cleans.fq.gz \
		-S "${SAMPLE}"/rrna.bt2.sam \
		2>&1 | tee "${SAMPLE}"/rrna.bowtie2.log

	bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x index/ath_rrna_updated \
		-1 "${SAMPLE}"/R1_cleans.fq.gz \
		-2 "${SAMPLE}"/R2_cleans.fq.gz \
		-S "${SAMPLE}"/rrna.updated.sam \
		2>&1 | tee "${SAMPLE}"/rrna.updated.bowtie2.log
done

# ==> G.max
vim gma_rrna.fa
cat ~/NJU_seq/data/gma_rrna/*s.fa >gma_rrna.origin.fa
bowtie-build gma_rrna.fa index/gma_rrna
bowtie-build gma_rrna.origin.fa index/gma_rrna_origin

for SAMPLE in NJU6256 NJU6268 NJUzp14; do
	bowtie -v 2 -a --best --strata \
		-S -x index/gma_rrna_origin \
		"${SAMPLE}"/R1_cleans.fq.gz \
		>"${SAMPLE}"/rrna_origin.sam
	bowtie -v 2 -a --best --strata \
		-S -x index/gma_rrna \
		"${SAMPLE}"/R1_cleans.fq.gz \
		>"${SAMPLE}"/rrna.sam
done

## reads processed: 23773878
## reads with at least one alignment: 4857517 (20.43%)
## reads that failed to align: 18916361 (79.57%)
#Reported 4857972 alignments
## reads processed: 23773878
## reads with at least one alignment: 4845849 (20.38%)
## reads that failed to align: 18928029 (79.62%)
#Reported 4847290 alignments
## reads processed: 8991030
## reads with at least one alignment: 1398801 (15.56%)
## reads that failed to align: 7592229 (84.44%)
#Reported 1398827 alignments
## reads processed: 8991030
## reads with at least one alignment: 1390773 (15.47%)
## reads that failed to align: 7600257 (84.53%)
#Reported 1391189 alignments
## reads processed: 16308091
## reads with at least one alignment: 3592250 (22.03%)
## reads that failed to align: 12715841 (77.97%)
#Reported 3592382 alignments
## reads processed: 16308091
## reads with at least one alignment: 3584212 (21.98%)
## reads that failed to align: 12723879 (78.02%)
#Reported 3584434 alignments

# select in 2
for SAMPLE in NJU6256 NJU6268 NJUzp14; do
	rm "${SAMPLE}"/rrna.sam
	samtools view -bS "${SAMPLE}"/rrna_origin.sam \
		| samtools sort - \
		| samtools view -h - \
		| awk '$1 ~ /^@/ || ($5 >= 20 && $6 !~ /S/)' \
		| samtools view -b -F 16 - \
		| samtools sort -o "${SAMPLE}"/rrna.fwd.clean.bam
	samtools index "${SAMPLE}"/rrna.fwd.clean.bam
done

for SAMPLE in NJU6256 NJU6268 NJUzp14; do
	bcftools mpileup \
		-f gma_rrna.origin.fa \
		-Q 20 -d 100000 \
		"${SAMPLE}"/rrna.fwd.clean.bam \
		| bcftools call -mv -Oz -o "${SAMPLE}"/rRNA.raw.vcf.gz
	bcftools norm -f gma_rrna.origin.fa \
		-m -any "${SAMPLE}"/rRNA.raw.vcf.gz \
		-Oz -o "${SAMPLE}"/rRNA.norm.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.norm.vcf.gz
	bcftools filter \
		-i 'INFO/DP>=200 && (INFO/DP4[2] + INFO/DP4[3]) / (INFO/DP4[0] + INFO/DP4[1] + INFO/DP4[2] + INFO/DP4[3]) > 0.7 && INFO/DP4[2] > 150' \
		"${SAMPLE}"/rRNA.norm.vcf.gz -Oz -o "${SAMPLE}"/rRNA.filt.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.filt.vcf.gz
done
bcftools isec \
	-n+2 -w1 \
	NJU6256/rRNA.filt.vcf.gz \
	NJU6268/rRNA.filt.vcf.gz \
	NJUzp14/rRNA.filt.vcf.gz \
	-Oz -o gma_rRNA.isec.vcf.gz
tabix -p vcf gma_rRNA.isec.vcf.gz

bowtie2-build gma_rrna.origin.fa index/gma_rrna_origin
for SAMPLE in NJU6256 NJU6268 NJUzp14; do
	bowtie2 --end-to-end \
		--rdg 5,3 --rfg 5,3 \
		--score-min L,-3,-2 \
		-N 1 -L 20 \
		-x index/gma_rrna_origin \
		-U "${SAMPLE}"/R1_cleanl.fq.gz \
		| samtools view -b -F 16 -q 20 \
		| samtools sort -o "${SAMPLE}"/rrna.indel.bam
done
for SAMPLE in NJU6256 NJU6268 NJUzp14; do
	bcftools mpileup \
		-f gma_rrna.origin.fa \
		-Q 20 -d 100000 \
		"${SAMPLE}"/rrna.indel.bam \
		| bcftools call -mv -Oz -o "${SAMPLE}"/rRNA.indel.vcf.gz
	bcftools norm -f gma_rrna.origin.fa \
		-m -any "${SAMPLE}"/rRNA.indel.vcf.gz \
		-Oz -o "${SAMPLE}"/rRNA.indel.norm.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.indel.norm.vcf.gz
	bcftools view \
		--types indels \
		"${SAMPLE}"/rRNA.indel.norm.vcf.gz \
		| bcftools filter \
			-i 'INFO/DP>=100 && INFO/IDV>=50 && INFO/IMF>=0.5' \
			-Oz -o "${SAMPLE}"/rRNA.indel.filt.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.indel.filt.vcf.gz
done
bcftools isec \
	-n+2 -w1 \
	NJU6256/rRNA.indel.filt.vcf.gz \
	NJU6268/rRNA.indel.filt.vcf.gz \
	NJUzp14/rRNA.indel.filt.vcf.gz \
	-Oz -o gma_rRNA.indel.isec.vcf.gz
tabix -p vcf gma_rRNA.indel.isec.vcf.gz

bcftools consensus -f gma_rrna.origin.fa gma_rRNA.isec.vcf.gz gma_rRNA.indel.isec.vcf.gz >gma_rrna.updated.fa

bowtie2-build gma_rrna.updated.fa index/gma_rrna_updated
for SAMPLE in NJU6256 NJU6268 NJUzp14; do
	bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x index/gma_rrna_origin \
		-1 "${SAMPLE}"/R1_cleans.fq.gz \
		-2 "${SAMPLE}"/R2_cleans.fq.gz \
		-S "${SAMPLE}"/rrna.bt2.sam \
		2>&1 | tee "${SAMPLE}"/rrna.bowtie2.log

	bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x index/gma_rrna_updated \
		-1 "${SAMPLE}"/R1_cleans.fq.gz \
		-2 "${SAMPLE}"/R2_cleans.fq.gz \
		-S "${SAMPLE}"/rrna.updated.sam \
		2>&1 | tee "${SAMPLE}"/rrna.updated.bowtie2.log
done

# ==> O.sat
vim osa_rrna.fa
cat ~/NJU_seq/data/osa_rrna/*s.fa >osa_rrna.origin.fa
bowtie-build osa_rrna.fa index/osa_rrna
bowtie-build osa_rrna.origin.fa index/osa_rrna_origin

for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	bowtie -v 2 -a --best --strata \
		-S -x index/osa_rrna_origin \
		"${SAMPLE}"/R1_cleans.fq.gz \
		>"${SAMPLE}"/rrna_origin.sam
	bowtie -v 2 -a --best --strata \
		-S -x index/osa_rrna \
		"${SAMPLE}"/R1_cleans.fq.gz \
		>"${SAMPLE}"/rrna.sam
done
## reads processed: 24850734
## reads with at least one alignment: 7805899 (31.41%)
## reads that failed to align: 17044835 (68.59%)
#Reported 7806694 alignments
## reads processed: 24850734
## reads with at least one alignment: 7640761 (30.75%)
## reads that failed to align: 17209973 (69.25%)
#Reported 7641879 alignments
## reads processed: 9916419
## reads with at least one alignment: 2768829 (27.92%)
## reads that failed to align: 7147590 (72.08%)
#Reported 2768889 alignments
## reads processed: 9916419
## reads with at least one alignment: 2738257 (27.61%)
## reads that failed to align: 7178162 (72.39%)
#Reported 2738324 alignments

# select in 2
for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	rm "${SAMPLE}"/rrna.sam
	samtools view -bS "${SAMPLE}"/rrna_origin.sam \
		| samtools sort - \
		| samtools view -h - \
		| awk '$1 ~ /^@/ || ($5 >= 20 && $6 !~ /S/)' \
		| samtools view -b -F 16 - \
		| samtools sort -o "${SAMPLE}"/rrna.fwd.clean.bam
	samtools index "${SAMPLE}"/rrna.fwd.clean.bam
done

for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	bcftools mpileup \
		-f osa_rrna.origin.fa \
		-Q 20 -d 100000 \
		"${SAMPLE}"/rrna.fwd.clean.bam \
		| bcftools call -mv -Oz -o "${SAMPLE}"/rRNA.raw.vcf.gz
	bcftools norm -f osa_rrna.origin.fa \
		-m -any "${SAMPLE}"/rRNA.raw.vcf.gz \
		-Oz -o "${SAMPLE}"/rRNA.norm.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.norm.vcf.gz
	bcftools filter \
		-i 'INFO/DP>=200 && (INFO/DP4[2] + INFO/DP4[3]) / (INFO/DP4[0] + INFO/DP4[1] + INFO/DP4[2] + INFO/DP4[3]) > 0.7 && INFO/DP4[2] > 150' \
		"${SAMPLE}"/rRNA.norm.vcf.gz -Oz -o "${SAMPLE}"/rRNA.filt.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.filt.vcf.gz
done
bcftools isec \
	-n+2 -w1 \
	NJU6260/rRNA.filt.vcf.gz \
	NJU6272/rRNA.filt.vcf.gz \
	-Oz -o osa_rRNA.isec.vcf.gz
tabix -p vcf osa_rRNA.isec.vcf.gz

bowtie2-build osa_rrna.origin.fa index/osa_rrna_origin
for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	bowtie2 --end-to-end \
		--rdg 5,3 --rfg 5,3 \
		--score-min L,-3,-2 \
		-N 1 -L 20 \
		-x index/osa_rrna_origin \
		-U "${SAMPLE}"/R1_cleanl.fq.gz \
		| samtools view -b -F 16 -q 20 \
		| samtools sort -o "${SAMPLE}"/rrna.indel.bam
done
for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	bcftools mpileup \
		-f osa_rrna.origin.fa \
		-Q 20 -d 100000 \
		"${SAMPLE}"/rrna.indel.bam \
		| bcftools call -mv -Oz -o "${SAMPLE}"/rRNA.indel.vcf.gz
	bcftools norm -f osa_rrna.origin.fa \
		-m -any "${SAMPLE}"/rRNA.indel.vcf.gz \
		-Oz -o "${SAMPLE}"/rRNA.indel.norm.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.indel.norm.vcf.gz
	bcftools view \
		--types indels \
		"${SAMPLE}"/rRNA.indel.norm.vcf.gz \
		| bcftools filter \
			-i 'INFO/DP>=100 && INFO/IDV>=50 && INFO/IMF>=0.5' \
			-Oz -o "${SAMPLE}"/rRNA.indel.filt.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.indel.filt.vcf.gz
done
bcftools isec \
	-n+2 -w1 \
	NJU6260/rRNA.indel.filt.vcf.gz \
	NJU6272/rRNA.indel.filt.vcf.gz \
	-Oz -o osa_rRNA.indel.isec.vcf.gz
tabix -p vcf osa_rRNA.indel.isec.vcf.gz

bcftools consensus -f osa_rrna.origin.fa osa_rRNA.isec.vcf.gz osa_rRNA.indel.isec.vcf.gz >osa_rrna.updated.fa

bowtie2-build osa_rrna.updated.fa index/osa_rrna_updated
for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x index/osa_rrna_origin \
		-1 "${SAMPLE}"/R1_cleans.fq.gz \
		-2 "${SAMPLE}"/R2_cleans.fq.gz \
		-S "${SAMPLE}"/rrna.bt2.sam \
		2>&1 | tee "${SAMPLE}"/rrna.bowtie2.log

	bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x index/osa_rrna_updated \
		-1 "${SAMPLE}"/R1_cleans.fq.gz \
		-2 "${SAMPLE}"/R2_cleans.fq.gz \
		-S "${SAMPLE}"/rrna.updated.sam \
		2>&1 | tee "${SAMPLE}"/rrna.updated.bowtie2.log
done

# ==> M.po
vim mpo_rrna.fa
cat ~/NJU_seq/data/mpo_rrna/*s.fa >mpo_rrna.origin.fa
bowtie-build mpo_rrna.fa index/mpo_rrna
bowtie-build mpo_rrna.origin.fa index/mpo_rrna_origin

for SAMPLE in NJUlymp; do
	bowtie -v 2 -a --best --strata \
		-S -x index/mpo_rrna_origin \
		"${SAMPLE}"/R1_cleans.fq.gz \
		>"${SAMPLE}"/rrna_origin.sam
	bowtie -v 2 -a --best --strata \
		-S -x index/mpo_rrna \
		"${SAMPLE}"/R1_cleans.fq.gz \
		>"${SAMPLE}"/rrna.sam
done

# select in 2
for SAMPLE in NJUlymp; do
	rm "${SAMPLE}"/rrna_origin.sam
	samtools view -bS "${SAMPLE}"/rrna.sam \
		| samtools sort - \
		| samtools view -h - \
		| awk '$1 ~ /^@/ || ($5 >= 20 && $6 !~ /S/)' \
		| samtools view -b -F 16 - \
		| samtools sort -o "${SAMPLE}"/rrna.fwd.clean.bam
	samtools index "${SAMPLE}"/rrna.fwd.clean.bam
done

for SAMPLE in NJUlymp; do
	bcftools mpileup \
		-f mpo_rrna.fa \
		-Q 20 -d 100000 \
		"${SAMPLE}"/rrna.fwd.clean.bam \
		| bcftools call -mv -Oz -o "${SAMPLE}"/rRNA.raw.vcf.gz
	bcftools norm -f mpo_rrna.fa \
		-m -any "${SAMPLE}"/rRNA.raw.vcf.gz \
		-Oz -o "${SAMPLE}"/rRNA.norm.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.norm.vcf.gz
	bcftools filter \
		-i 'INFO/DP>=200 && (INFO/DP4[2] + INFO/DP4[3]) / (INFO/DP4[0] + INFO/DP4[1] + INFO/DP4[2] + INFO/DP4[3]) > 0.7 && INFO/DP4[2] > 150' \
		"${SAMPLE}"/rRNA.norm.vcf.gz -Oz -o "${SAMPLE}"/rRNA.filt.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.filt.vcf.gz
done

bowtie2-build mpo_rrna.fa index/mpo_rrna
for SAMPLE in NJUlymp; do
	bowtie2 --end-to-end \
		--rdg 5,3 --rfg 5,3 \
		--score-min L,-3,-2 \
		-N 1 -L 20 \
		-x index/mpo_rrna \
		-U "${SAMPLE}"/R1_cleanl.fq.gz \
		| samtools view -b -F 16 -q 20 \
		| samtools sort -o "${SAMPLE}"/rrna.indel.bam
done
for SAMPLE in NJUlymp; do
	bcftools mpileup \
		-f mpo_rrna.fa \
		-Q 20 -d 100000 \
		"${SAMPLE}"/rrna.indel.bam \
		| bcftools call -mv -Oz -o "${SAMPLE}"/rRNA.indel.vcf.gz
	bcftools norm -f mpo_rrna.fa \
		-m -any "${SAMPLE}"/rRNA.indel.vcf.gz \
		-Oz -o "${SAMPLE}"/rRNA.indel.norm.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.indel.norm.vcf.gz
	bcftools view \
		--types indels \
		"${SAMPLE}"/rRNA.indel.norm.vcf.gz \
		| bcftools filter \
			-i 'INFO/DP>=100 && INFO/IDV>=50 && INFO/IMF>=0.5' \
			-Oz -o "${SAMPLE}"/rRNA.indel.filt.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.indel.filt.vcf.gz
done

bcftools consensus -f mpo_rrna.fa NJUlymp/rRNA.filt.vcf.gz NJUlymp/rRNA.indel.filt.vcf.gz >mpo_rrna.updated.fa
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP4]\n' rRNA.isec.vcf.gz >rRNA.snp.tsv
#perl update_rrna_ref.pl mpo_rrna.fa rRNA.snp.tsv >mpo_rrna.updated.fa

bowtie2-build mpo_rrna.updated.fa index/mpo_rrna_updated
for SAMPLE in NJUlymp; do
	bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x index/mpo_rrna \
		-1 "${SAMPLE}"/R1_cleans.fq.gz \
		-2 "${SAMPLE}"/R2_cleans.fq.gz \
		-S "${SAMPLE}"/rrna.bt2.sam \
		2>&1 | tee "${SAMPLE}"/rrna.bowtie2.log

	bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x index/mpo_rrna_updated \
		-1 "${SAMPLE}"/R1_cleans.fq.gz \
		-2 "${SAMPLE}"/R2_cleans.fq.gz \
		-S "${SAMPLE}"/rrna.updated.sam \
		2>&1 | tee "${SAMPLE}"/rrna.updated.bowtie2.log
done

# ==> S.lyr
cat ~/NJU_seq/data/sly_rrna/*s.fa >sly_rrna.origin.fa
bowtie-build sly_rrna.origin.fa index/sly_rrna_origin

for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	bowtie -v 2 -a --best --strata \
		-S -x index/sly_rrna_origin \
		"${SAMPLE}"/R1_cleans.fq.gz \
		>"${SAMPLE}"/rrna_origin.sam
done

## reads processed: 22519825
## reads with at least one alignment: 6977777 (30.99%)
## reads that failed to align: 15542048 (69.01%)
#Reported 6978061 alignments
## reads processed: 8913190
## reads with at least one alignment: 1567139 (17.58%)
## reads that failed to align: 7346051 (82.42%)
#Reported 1567168 alignments
## reads processed: 33119181
## reads with at least one alignment: 20262353 (61.18%)
## reads that failed to align: 12856828 (38.82%)
#Reported 20262715 alignments

# select in 2
for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	samtools view -bS "${SAMPLE}"/rrna_origin.sam \
		| samtools sort - \
		| samtools view -h - \
		| awk '$1 ~ /^@/ || ($5 >= 20 && $6 !~ /S/)' \
		| samtools view -b -F 16 - \
		| samtools sort -o "${SAMPLE}"/rrna.fwd.clean.bam
	samtools index "${SAMPLE}"/rrna.fwd.clean.bam
done

for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	bcftools mpileup \
		-f sly_rrna.origin.fa \
		-Q 20 -d 100000 \
		"${SAMPLE}"/rrna.fwd.clean.bam \
		| bcftools call -mv -Oz -o "${SAMPLE}"/rRNA.raw.vcf.gz
	bcftools norm -f sly_rrna.origin.fa \
		-m -any "${SAMPLE}"/rRNA.raw.vcf.gz \
		-Oz -o "${SAMPLE}"/rRNA.norm.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.norm.vcf.gz
	bcftools filter \
		-i 'INFO/DP>=200 && (INFO/DP4[2] + INFO/DP4[3]) / (INFO/DP4[0] + INFO/DP4[1] + INFO/DP4[2] + INFO/DP4[3]) > 0.7 && INFO/DP4[2] > 150' \
		"${SAMPLE}"/rRNA.norm.vcf.gz -Oz -o "${SAMPLE}"/rRNA.filt.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.filt.vcf.gz
done
bcftools isec \
	-n+2 -w1 \
	NJU6252/rRNA.filt.vcf.gz \
	NJU6264/rRNA.filt.vcf.gz \
	NJUlysly/rRNA.filt.vcf.gz \
	-Oz -o sly_rRNA.isec.vcf.gz
tabix -p vcf sly_rRNA.isec.vcf.gz

bowtie2-build sly_rrna.origin.fa index/sly_rrna_origin
for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	bowtie2 --end-to-end \
		--rdg 5,3 --rfg 5,3 \
		--score-min L,-3,-2 \
		-N 1 -L 20 \
		-x index/sly_rrna_origin \
		-U "${SAMPLE}"/R1_cleanl.fq.gz \
		| samtools view -b -F 16 -q 20 \
		| samtools sort -o "${SAMPLE}"/rrna.indel.bam
done
for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	bcftools mpileup \
		-f sly_rrna.origin.fa \
		-Q 20 -d 100000 \
		"${SAMPLE}"/rrna.indel.bam \
		| bcftools call -mv -Oz -o "${SAMPLE}"/rRNA.indel.vcf.gz
	bcftools norm -f sly_rrna.origin.fa \
		-m -any "${SAMPLE}"/rRNA.indel.vcf.gz \
		-Oz -o "${SAMPLE}"/rRNA.indel.norm.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.indel.norm.vcf.gz
	bcftools view \
		--types indels \
		"${SAMPLE}"/rRNA.indel.norm.vcf.gz \
		| bcftools filter \
			-i 'INFO/DP>=100 && INFO/IDV>=50 && INFO/IMF>=0.5' \
			-Oz -o "${SAMPLE}"/rRNA.indel.filt.vcf.gz
	tabix -p vcf "${SAMPLE}"/rRNA.indel.filt.vcf.gz
done
bcftools isec \
	-n+2 -w1 \
	NJU6252/rRNA.indel.filt.vcf.gz \
	NJU6264/rRNA.indel.filt.vcf.gz \
	NJUlysly/rRNA.indel.filt.vcf.gz \
	-Oz -o sly_rRNA.indel.isec.vcf.gz
tabix -p vcf sly_rRNA.indel.isec.vcf.gz

bcftools consensus -f sly_rrna.origin.fa sly_rRNA.isec.vcf.gz sly_rRNA.indel.isec.vcf.gz >sly_rrna.updated.fa

bowtie2-build sly_rrna.updated.fa index/sly_rrna_updated
for SAMPLE in NJU6252 NJU6264 NJUlysly; do
	bowtie2 -p 16 -k 10 -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x index/sly_rrna_origin \
		-1 "${SAMPLE}"/R1_cleans.fq.gz \
		-2 "${SAMPLE}"/R2_cleans.fq.gz \
		-S "${SAMPLE}"/rrna.bt2.sam \
		2>&1 | tee "${SAMPLE}"/rrna.bowtie2.log

	bowtie2 -p 16 -k 10 -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x index/sly_rrna_updated \
		-1 "${SAMPLE}"/R1_cleans.fq.gz \
		-2 "${SAMPLE}"/R2_cleans.fq.gz \
		-S "${SAMPLE}"/rrna.updated.sam \
		2>&1 | tee "${SAMPLE}"/rrna.updated.bowtie2.log
done