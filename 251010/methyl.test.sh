mkdir lambda_pUC19_index Col-0_index
mv lambda_pUC19.fa lambda_pUC19_index/.
cd lambda_pUC19_index || return
bismark_genome_preparation --bowtie2 --verbose --parallel 6 --genomic_composition ./

for SAMPLE in Osa_WGBS Osa_EM; do
	mkdir ${SAMPLE}
	trim_galore --paired --quality 20 \
		--illumina --fastqc --gzip -j 8 -a2 A\{15\} \
		-o ${SAMPLE}/. \
		--clip_R1 5 --clip_R2 15 \
		--three_prime_clip_R1 15 --three_prime_clip_R2 5 \
		/home/ivan/cater_data/nlr_data/*/*/${SAMPLE}*_1.fq.gz \
		/home/ivan/cater_data/nlr_data/*/*/${SAMPLE}*_2.fq.gz
	mv ${SAMPLE}/*val_1.fq.gz ${SAMPLE}/R1.fq.gz
	mv ${SAMPLE}/*val_2.fq.gz ${SAMPLE}/R2.fq.gz
done

for SAMPLE in Osa_ACE Col_ACE; do
	mkdir ${SAMPLE}
	ln -sf /home/ivan/cater_data/nlr_data/*/*/02.clean/${SAMPLE}*.1.fq.gz \
		/home/ivan/fat/methylseq/${SAMPLE}/R1.fq.gz
	ln -sf /home/ivan/cater_data/nlr_data/*/*/02.clean/${SAMPLE}*.2.fq.gz \
		/home/ivan/fat/methylseq/${SAMPLE}/R2.fq.gz
done

for SAMPLE in Osa_WGBS Osa_EM Osa_ACE Col_ACE; do
	bismark \
		-p 7 -X 1200 -I 50 --gzip -q \
		--genome lambda_index/ \
		-B ${SAMPLE}/lambda \
		-1 ${SAMPLE}/R1.fq.gz \
		-2 ${SAMPLE}/R2.fq.gz
	bismark \
		-p 7 -X 1200 -I 50 --gzip -q \
		--genome pUC19_index/ \
		-B ${SAMPLE}/pUC19 \
		-1 ${SAMPLE}/R1.fq.gz \
		-2 ${SAMPLE}/R2.fq.gz
done

for SAMPLE in Osa_WGBS Osa_EM Osa_ACE; do
	bismark \
		-p 7 -X 1200 -I 50 --gzip -q \
		--genome Osa_index/ \
		-B ${SAMPLE}/genome \
		-1 ${SAMPLE}/R1.fq.gz \
		-2 ${SAMPLE}/R2.fq.gz
	deduplicate_bismark \
		${SAMPLE}/genome_pe.bam \
		--output_dir ${SAMPLE}
	bismark_methylation_extractor \
		--bedGraph --cytosine_report --gzip \
		--multicore 7 --CX --buffer_size 20G --counts \
		--genome_folder /home/ivan/fat/methylseq/Osa_index/ \
		${SAMPLE}/genome_pe.deduplicated.bam \
		--output_dir ${SAMPLE}
	samtools view \
		-f 2 ${SAMPLE}/genome_pe.deduplicated.bam \
		| awk '{print sqrt(($9)^2)}' >${SAMPLE}/fragment_length.txt
done

for SAMPLE in Col_ACE; do
	bismark \
		-p 7 -X 1200 -I 50 --gzip -q \
		--genome Ath_Col-0_index/ \
		-B ${SAMPLE}/genome \
		-1 ${SAMPLE}/R1.fq.gz \
		-2 ${SAMPLE}/R2.fq.gz
	deduplicate_bismark \
		${SAMPLE}/genome_pe.bam \
		--output_dir ${SAMPLE}
	bismark_methylation_extractor \
		--bedGraph --cytosine_report --gzip \
		--multicore 7 --CX --buffer_size 20G --counts \
		--genome_folder /home/ivan/fat/methylseq/Ath_Col-0_index/ \
		${SAMPLE}/genome_pe.deduplicated.bam \
		--output_dir ${SAMPLE}
	samtools view \
		-f 2 ${SAMPLE}/genome_pe.deduplicated.bam \
		| awk '{print sqrt(($9)^2)}' >${SAMPLE}/fragment_length.txt
done

for SAMPLE in Col_ACE; do
	bismark \
		-p 7 -X 1200 -I 50 --gzip -q \
		--genome Ath_Tair10_index/ \
		-B ${SAMPLE}/Tair10 \
		-1 ${SAMPLE}/R1.fq.gz \
		-2 ${SAMPLE}/R2.fq.gz
	bismark \
		-p 7 -X 1200 -I 50 --gzip -q \
		--genome Ath_Col-PEK_index/ \
		-B ${SAMPLE}/ColPEK \
		-1 ${SAMPLE}/R1.fq.gz \
		-2 ${SAMPLE}/R2.fq.gz
done

for SAMPLE in Osa_WGBS Osa_EM Osa_ACE; do
	bismark_methylation_extractor --comprehensive --bedGraph \
		--multicore 7 --buffer_size 20G --counts \
		--genome_folder Osa_index --cytosine_report --gzip \
		Osa_${platform}_selected/bismark_pe.deduplicated.bam \
		--output_dir Osa_${platform}_selected
	samtools view -f 2 Osa_${platform}_selected/bismark_pe.deduplicated.bam \
		| awk '{print sqrt(($9)^2)}' >Osa_${platform}_selected/fragment_length.txt
done

pigz -dc Osativa_323_v7.0.gene_exons.gff3.gz | awk '$3=="gene"' \
	| awk '{ match($9,/ID=([^;]+)/,a); if(a[1]!="") print $1"\t"$4-1"\t"$5"\t"a[1]; }' \
	| sort -k1,1 -k2,2n >Osa.genes.bed
sed -i 's/\.MSUv7\.0//' Osa.genes.bed
grep -F -f <(grep -v "No_NLR" 1st_blast_classification.txt | grep "^LOC" | awk -F'.' '{print $1}') \
	Osa.genes.bed \
	| sort -k1,1 -k2,2n \
		>Osa.NLR.genes.bed
pigz -dc Osativa_323_v7.0.gene_exons.gff3.gz \
	| awk '
		BEGIN{OFS="\t"}
		$3=="gene"{
			chr=$1
			strand=$7
			if (strand=="+") {
				tss=$4
				start=tss-2000-1
				end=tss
			} else if (strand=="-") {
				tss=$5
				start=tss+1
				end=tss+2000
			}
			if (start<0) start=0
			match($9,/ID=([^;]+)/,a)
			print chr, start, end, a[1]
	}' | sort -k1,1 -k2,2n >Osa.prom.bed
sed -i 's/\.MSUv7\.0//' Osa.prom.bed
grep -F -f <(grep -v "No_NLR" 1st_blast_classification.txt | grep "^LOC" | awk -F'.' '{print $1}') \
	Osa.prom.bed \
	| sort -k1,1 -k2,2n \
		>Osa.NLR.prom.bed

for SAMPLE in Osa_ACE Col_ACE; do
	mkdir ${SAMPLE}
	ln -sf /home/ivan/cater_data/nlr_data/*/*/02.clean/${SAMPLE}*.1.fq.gz /home/ivan/fat/methylseq/${SAMPLE}/R1.fq.gz
	ln -sf /home/ivan/cater_data/nlr_data/*/*/02.clean/${SAMPLE}*.2.fq.gz /home/ivan/fat/methylseq/${SAMPLE}/R2.fq.gz
done

bismark --genome Osa_index/ \
	-p 7 -X 1000 -I 60 \
	-B Osa_ACE/bismark --gzip -q \
	-1 Osa_ACE/R1.fq.gz \
	-2 Osa_ACE/R2.fq.gz

deduplicate_bismark \
	Osa_ACE/bismark_pe.bam \
	--output_dir Osa_ACE
bismark_methylation_extractor --comprehensive --bedGraph \
	--multicore 7 --buffer_size 20G --counts \
	--genome_folder Osa_index --cytosine_report --gzip \
	Osa_ACE/bismark_pe.deduplicated.bam \
	--output_dir Osa_ACE
samtools view -f 2 Osa_ACE/bismark_pe.deduplicated.bam \
	| awk '{print sqrt(($9)^2)}' >Osa_ACE/fragment_length.txt
