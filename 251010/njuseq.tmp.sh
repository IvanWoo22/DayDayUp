for SAMPLE in NJU609{2..5} NJU61{08..11}; do
	mkdir "${SAMPLE}"
	ln -sf /home/ivan/cater_data/"${SAMPLE}"*R1*.fastq.gz /home/ivan/fat/NJU_seq/"${SAMPLE}"/R1.fq.gz
	ln -sf /home/ivan/cater_data/"${SAMPLE}"*R2*.fastq.gz /home/ivan/fat/NJU_seq/"${SAMPLE}"/R2.fq.gz
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAG \
		-O 6 -m 10 -e 0.1 --discard-untrimmed \
		-o "${SAMPLE}"/R1_clean.fq.gz -p "${SAMPLE}"/R2_clean.fq.gz \
		-j 16 "${SAMPLE}"/R1.fq.gz "${SAMPLE}"/R2.fq.gz \
		>"${SAMPLE}"/cutadapt.log 2>&1
done

for SAMPLE in NJU609{2..5}; do
	bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x temp_index/selected_miRNA_A549 \
		-1 "${SAMPLE}"/R1_clean.fq.gz \
		-2 "${SAMPLE}"/R2_clean.fq.gz \
		-S "${SAMPLE}"/raw.sam \
		2>&1 | tee "${SAMPLE}"/bowtie2.log
done

for SAMPLE in NJU61{08..11}; do
	bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x temp_index/selected_miRNA_piRNA_mmu \
		-1 "${SAMPLE}"/R1_clean.fq.gz \
		-2 "${SAMPLE}"/R2_clean.fq.gz \
		-S "${SAMPLE}"/raw.sam \
		2>&1 | tee "${SAMPLE}"/bowtie2.log
done

for SAMPLE in NJU609{2..5} NJU61{08..11}; do
	awk '$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}' "${SAMPLE}"/raw.sam \
		| perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl \
		| perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl \
			>"${SAMPLE}"/rna.out.tmp
done

for SAMPLE in NJU609{2..5}; do
	for RNA in $(ls A549_ref/*.fa | sed 's/A549_ref\///' | sed 's/.fa//'); do
		perl ~/NJU_seq/rrna_analysis/readend_count.pl \
			A549_ref/${RNA}.fa \
			${SAMPLE}/rna.out.tmp ${RNA} \
			>${SAMPLE}/${RNA}.tsv
	done
done

for SAMPLE in NJU61{08..11}; do
	for RNA in $(ls mmu_ref/*.fa | sed 's/mmu_ref\///' | sed 's/.fa//'); do
		perl ~/NJU_seq/rrna_analysis/readend_count.pl \
			mmu_ref/${RNA}.fa \
			${SAMPLE}/rna.out.tmp ${RNA} \
			>${SAMPLE}/${RNA}.tsv
	done
done

time parallel -j 3 "
	perl scoretmp.pl \
		NJU6092/{}.tsv \
		NJU6093/{}.tsv \
		NJU6094/{}.tsv \
		NJU6095/{}.tsv \
		>A549.rna.scored.Nm.{}.tsv
" :::: <(ls A549_ref/*.fa | sed 's/A549_ref\///' | sed 's/.fa//')

time parallel -j 3 "
	perl scoretmp.pl \
		NJU6108/{}.tsv \
		NJU6109/{}.tsv \
		NJU6110/{}.tsv \
		NJU6111/{}.tsv \
		>mmu.rna.scored.Nm.{}.tsv
" :::: <(ls mmu_ref/*.fa | sed 's/mmu_ref\///' | sed 's/.fa//')

for RNA in $(ls A549_ref/*.fa | sed 's/A549_ref\///' | sed 's/.fa//'); do
	awk -v RNA="${RNA}" '{print RNA "\t" $0}' A549.rna.scored.Nm."${RNA}".tsv >>A549.rna.scored.Nm.tsv
done

for RNA in $(ls mmu_ref/*.fa | sed 's/mmu_ref\///' | sed 's/.fa//'); do
	awk -v RNA="${RNA}" '{print RNA "\t" $0}' mmu.rna.scored.Nm."${RNA}".tsv >>mmu.rna.scored.Nm.tsv
done

for SAMPLE in NJU609{2..5}; do
	samtools view -h -F 4 "${SAMPLE}"/raw.sam \
		| awk '($0 ~ /^@/) || ($0 ~ /NM:i:0/)' \
		| samtools sort -o "${SAMPLE}"/raw.sorted.bam
	samtools index "${SAMPLE}"/raw.sorted.bam
done

for SAMPLE in NJU609{2..5}; do
	cp "${SAMPLE}"/raw.sorted.bam "${SAMPLE}".bam
	samtools index "${SAMPLE}".bam
done
