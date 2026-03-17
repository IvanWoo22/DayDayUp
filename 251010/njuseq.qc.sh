NJU62{{32..35},{48..75}}
for dir in NJU62{{32..35},{48..75}} NJUzp{3..18}; do
	echo -n "$dir: "
	[ -d /home/ivan/cater_data/"$dir" ] && ls /home/ivan/cater_data/"$dir"/*.gz 2>/dev/null | wc -l || echo "文件夹不存在"
done

for SAMPLE in NJU62{{32..35},{48..75}} NJUzp{3..18}; do
	mkdir "${SAMPLE}"
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R1*.f*q.gz /home/ivan/fat/NJU_seq_fast/"${SAMPLE}"/R1.fq.gz
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R2*.f*q.gz /home/ivan/fat/NJU_seq_fast/"${SAMPLE}"/R2.fq.gz
done

for SAMPLE in Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3}; do
	mkdir "${SAMPLE}"
	ln -sf /home/ivan/cater_data/nlr_data/njuseq/*/"${SAMPLE}"/*_1.fq.gz /home/ivan/fat/NJU_seq_fast/"${SAMPLE}"/R1.fq.gz
	ln -sf /home/ivan/cater_data/nlr_data/njuseq/*/"${SAMPLE}"/*_2.fq.gz /home/ivan/fat/NJU_seq_fast/"${SAMPLE}"/R2.fq.gz
done

for SAMPLE in NJU62{{32..35},{48..75}} NJUzp{3..18} Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3}; do
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAG \
		-O 6 -m 10 -e 0.1 --discard-untrimmed \
		-o "${SAMPLE}"/R1_clean.fq.gz -p "${SAMPLE}"/R2_clean.fq.gz \
		-j 16 "${SAMPLE}"/R1.fq.gz "${SAMPLE}"/R2.fq.gz \
		>"${SAMPLE}"/cutadapt.log 2>&1
done

for SAMPLE in NJU62{{32..35},{48..75}} NJUzp{3..18} Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3}; do
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAG \
		-O 6 -m 30 -e 0.1 --discard-untrimmed \
		-o "${SAMPLE}"/R1_long.fq.gz -p "${SAMPLE}"/R2_long.fq.gz \
		-j 16 "${SAMPLE}"/R1.fq.gz "${SAMPLE}"/R2.fq.gz \
		>"${SAMPLE}"/cutadapt_long.log 2>&1
done

for SAMPLE in Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3} NJU62{{32..35},{48..75}} NJUzp{3..18}; do
	bsub -n 24 -q largemem -J "$SAMPLE" "
	/share/home/wangq/miniconda3/bin/kraken2 \
		--threads 20 --use-names \
		--gzip-compressed --paired \
		--report ${SAMPLE}/class.report \
		--classified-out ${SAMPLE}/temp#.fq \
		--db /share/home/wangq/wyf/kraken_plant \
		${SAMPLE}/R1_long.fq.gz \
		${SAMPLE}/R2_long.fq.gz \
		--output ${SAMPLE}/class.tsv
		"
done

awk '$3=="rRNA" {
  OFS="\t";
  print $1, $4-1, $5, $9, ".", $7
}' GCF_019754155.1_ASM1975415v1_genomic.gff >ecoli_rrna.bed
awk '$3=="tRNA" {
  OFS="\t";
  print $1, $4-1, $5, $9, ".", $7
}' GCF_019754155.1_ASM1975415v1_genomic.gff >ecoli_trna.bed
awk '($3=="ncRNA" || $3=="tmRNA" || $3=="RNase_P_RNA") {
  OFS="\t";
  print $1, $4-1, $5, $9, ".", $7
}' GCF_019754155.1_ASM1975415v1_genomic.gff >ecoli_ncrna.bed
bedtools getfasta \
	-fi GCF_019754155.1_ASM1975415v1_genomic.fna \
	-bed ecoli_rrna.bed -s -name \
	>ecoli_rrna.fa
bedtools getfasta \
	-fi GCF_019754155.1_ASM1975415v1_genomic.fna \
	-bed ecoli_trna.bed -s -name \
	>ecoli_trna.fa
bedtools getfasta \
	-fi GCF_019754155.1_ASM1975415v1_genomic.fna \
	-bed ecoli_ncrna.bed -s -name \
	>ecoli_ncrna.fa
seqkit stats ecoli_rrna.fa ecoli_trna.fa ecoli_ncrna.fa
#processed files:  3 / 3 [======================================] ETA: 0s. done
#file            format  type  num_seqs  sum_len  min_len  avg_len  max_len
#ecoli_rrna.fa   FASTA   DNA          8    5,143      115    642.9    2,906
#ecoli_trna.fa   FASTA   DNA         78    6,171       64     79.1      105
#ecoli_ncrna.fa  FASTA   DNA         11    2,029      132    184.5      377
cat ecoli_rrna.fa ecoli_trna.fa ecoli_ncrna.fa \
	| seqkit rmdup -s \
		>ecoli_rna.fa
rm ecoli_rrna.* ecoli_trna.* ecoli_ncrna.*
bowtie2-build ecoli_rna.fa ../index/ecoli_rna

zgrep -w "Enterobacteriaceae" tax_slv_ssu_138.2.txt.gz | cut -f1 | sort -u >entero_ssu.ids
zgrep -w "Enterobacteriaceae" tax_slv_lsu_138.2.txt.gz | cut -f1 | sort -u >entero_lsu.ids
zgrep -w "Lactobacillaceae" tax_slv_ssu_138.2.txt.gz | cut -f1 | sort -u >lacto_ssu.ids
zgrep -w "Lactobacillaceae" tax_slv_lsu_138.2.txt.gz | cut -f1 | sort -u >lacto_lsu.ids
zgrep -w "Streptococcaceae" tax_slv_ssu_138.2.txt.gz | cut -f1 | sort -u >strep_ssu.ids
zgrep -w "Streptococcaceae" tax_slv_lsu_138.2.txt.gz | cut -f1 | sort -u >strep_lsu.ids
seqkit grep \
	-n -r -f <(cat ./*_ssu.ids) \
	SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz \
	-o bac_ssu.fa.gz
seqkit grep \
	-n -r -f <(cat ./*_lsu.ids) \
	SILVA_138.2_LSURef_NR99_tax_silva_trunc.fasta.gz \
	-o bac_lsu.fa.gz
pigz -dc bac_ssu.fa.gz bac_lsu.fa.gz | sed '/^>/! s/U/T/g' >bac_rna.fa
rm ./*_*su.fa.gz ./*_*su.ids
bowtie2-build --threads 12 bac_rna.fa ../index/bac_rna

for SAMPLE in NJU62{56..59}; do
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R1*.f*q.gz /home/ivan/fat/NJU_seq_qc/"${SAMPLE}"/R1.fq.gz
	ln -sf /home/ivan/cater_data/"${SAMPLE}"/*_R2*.f*q.gz /home/ivan/fat/NJU_seq_qc/"${SAMPLE}"/R2.fq.gz
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAG \
		-O 6 -m 10 -e 0.1 --discard-untrimmed -o "${SAMPLE}"/R1.p1.fq.gz -p "${SAMPLE}"/R2.p1.fq.gz \
		-j 16 "${SAMPLE}"/R1.fq.gz "${SAMPLE}"/R2.fq.gz \
		>"${SAMPLE}"/cutadapt.log 2>&1
done

for SAMPLE in NJU62{56..59}; do
	bowtie2 --end-to-end -p 16 -k 20 -t \
		--no-unal --no-mixed --no-discordant \
		-D 20 -R 3 -N 0 -L 10 -i S,1,0.75 \
		--maxins 120 --score-min L,1.6,-0.8 \
		--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \
		--xeq -x ../index/ecoli_rna \
		-1 "${SAMPLE}"/R1.p1.fq.gz \
		-2 "${SAMPLE}"/R2.p1.fq.gz \
		-S "${SAMPLE}"/ecoli_rna.raw.sam \
		2>&1 | tee "${SAMPLE}"/ecoli_rna.bowtie2.log
done

parallel -j 8 --keep-order '
	samtools view -F 256 {1}/ecoli_rna.raw.bam \
		| awk '\''$6!="*"&&$7=="="&&$4==$8{print $1}'\'' \
		| perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl \
			>{1}/ecoli_rna.out.list
	perl ~/NJU_seq/tool/delete_fastq.pl -n {1}/ecoli_rna.out.list -i {1}/R1.p1.fq.gz -o {1}/R1.p2.fq.gz &
	perl ~/NJU_seq/tool/delete_fastq.pl -n {1}/ecoli_rna.out.list -i {1}/R2.p1.fq.gz -o {1}/R2.p2.fq.gz &
	wait
' ::: NJU62{56..59}

for SAMPLE in NJU62{56..59}; do
	time bowtie2 --end-to-end -p 16 -k 20 -t \
		--no-unal --no-mixed --no-discordant \
		-D 20 -R 3 -N 0 -L 10 -i S,1,0.75 \
		--maxins 120 --score-min L,1.6,-0.8 \
		--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \
		--xeq -x ../index/lacto_rna \
		-1 "${SAMPLE}"/R1.p2.fq.gz \
		-2 "${SAMPLE}"/R2.p2.fq.gz \
		-S "${SAMPLE}"/lacto_rna.raw.sam \
		2>&1 | tee "${SAMPLE}"/lacto_rna.bowtie2.log
done

parallel -j 8 --keep-order '
	samtools view -F 256 {1}/lacto_rna.raw.sam \
		| awk '\''$6!="*"&&$7=="="&&$4==$8{print $1}'\'' \
		| perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl \
			>{1}/lacto_rna.out.list
	perl ~/NJU_seq/tool/delete_fastq.pl -n {1}/lacto_rna.out.list -i {1}/R1.p2.fq.gz -o {1}/R1.p3.fq.gz &
	perl ~/NJU_seq/tool/delete_fastq.pl -n {1}/lacto_rna.out.list -i {1}/R2.p2.fq.gz -o {1}/R2.p3.fq.gz &
	wait
' ::: NJU62{56..59}

for SAMPLE in NJU62{56..59}; do
	seqkit seq -g -m 30 "${SAMPLE}"/R1.p3.fq.gz >"${SAMPLE}"/R1.tmp.fq
	seqkit seq -g -m 30 "${SAMPLE}"/R2.p3.fq.gz >"${SAMPLE}"/R2.tmp.fq
	seqkit pair -1 "${SAMPLE}"/R1.tmp.fq -2 "${SAMPLE}"/R2.tmp.fq
	pigz <"${SAMPLE}"/R1.tmp.paired.fq >"${SAMPLE}"/R1.clean_long.fq.gz
	pigz <"${SAMPLE}"/R2.tmp.paired.fq >"${SAMPLE}"/R2.clean_long.fq.gz
	rm "${SAMPLE}"/R*.tmp*.fq
done

##==> NJU6256/long
# 18.47	379542	0	P	35493	        Streptophyta
# 11.60	238405	5851	P	1239	        Bacillota
#  0.85	17464	884	P	1224	      Pseudomonadota
##==> NJU6256/p1
#  1.46	379654	0	P	35493	        Streptophyta
#  0.92	238441	5851	P	1239	        Bacillota
#  0.07	17477	885	P	1224	      Pseudomonadota
##==> NJU6256/p3
#  1.80	372638	0	P	35493	        Streptophyta
#  0.08	17418	414	P	1239	        Bacillota
#  0.08	15686	852	P	1224	      Pseudomonadota
##==> NJU6257/long
# 14.78	497252	8562	P	1224	      Pseudomonadota
#  4.84	162807	4183	P	1239	        Bacillota
#  7.97	268082	0	P	35493	        Streptophyta
#  0.91	497337	8566	P	1224	      Pseudomonadota
#  0.30	162824	4183	P	1239	        Bacillota
#  0.49	268293	0	P	35493	        Streptophyta
#  0.85	343018	2564	P	1224	      Pseudomonadota
#  0.03	14026	314	P	1239	        Bacillota
#  0.66	267519	0	P	35493	        Streptophyta
##==> NJU6258/long
# 15.45	568590	15026	P	1239	        Bacillota
# 13.66	502573	12318	P	1224	      Pseudomonadota
#  4.57	168138	0	P	35493	        Streptophyta
#  1.06	568740	15032	P	1239	        Bacillota
#  0.94	502722	12322	P	1224	      Pseudomonadota
#  0.31	168365	0	P	35493	        Streptophyta
#  0.94	336017	3709	P	1224	      Pseudomonadota
#  0.12	41967	980	P	1239	        Bacillota
#  0.47	166953	0	P	35493	        Streptophyta
##==> NJU6259/long
# 15.30	558794	11641	P	1224	      Pseudomonadota
# 14.25	520384	11410	P	1239	        Bacillota
#  4.78	174504	0	P	35493	        Streptophyta
#  0.86	558886	11649	P	1224	      Pseudomonadota
#  0.80	520443	11412	P	1239	        Bacillota
#  0.27	174628	0	P	35493	        Streptophyta
#  0.85	376016	3699	P	1224	      Pseudomonadota
#  0.08	34488	816	P	1239	        Bacillota
#  0.39	173582	0	P	35493	        Streptophyta

for SAMPLE in NJU62{{32..35},{48..75}} NJUzp{3..18} Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3}; do
	time bowtie2 --end-to-end -p 16 -a -t \
		--no-unal --no-mixed --no-discordant \
		-D 20 -R 3 -N 0 -L 10 -i S,1,0.75 \
		--maxins 120 --score-min L,-0.8,-0.8 \
		--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \
		--xeq -x ../index/ecoli_rna \
		-1 "${SAMPLE}"/R1_clean.fq.gz \
		-2 "${SAMPLE}"/R2_clean.fq.gz \
		-S "${SAMPLE}"/ecoli_rna.raw.sam \
		2>&1 | tee "${SAMPLE}"/ecoli_rna.bowtie2.log
done

parallel -j 10 --keep-order '
	samtools view -F 256 {1}/ecoli_rna.raw.sam \
		| awk '\''$6!="*"&&$7=="="&&$4==$8{print $1}'\'' \
		| perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl \
			>{1}/ecoli_rna.out.list
	rm {1}/ecoli_rna.raw.sam
	perl ~/NJU_seq/tool/delete_fastq.pl -n {1}/ecoli_rna.out.list -i {1}/R1_clean.fq.gz -o {1}/R1.p2.fq.gz &
	perl ~/NJU_seq/tool/delete_fastq.pl -n {1}/ecoli_rna.out.list -i {1}/R2_clean.fq.gz -o {1}/R2.p2.fq.gz &
	wait
' ::: NJU62{{32..35},{48..75}} NJUzp{3..18} Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3}

for SAMPLE in NJU62{{32..35},{48..75}} NJUzp{3..18} Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3}; do
	time bowtie2 --end-to-end -p 16 -k 20 -t \
		--no-unal --no-mixed --no-discordant \
		-D 20 -R 3 -N 0 -L 10 -i S,1,0.75 \
		--maxins 120 --score-min L,1.6,-0.8 \
		--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \
		--xeq -x ../index/lacto_rna \
		-1 "${SAMPLE}"/R1.p2.fq.gz \
		-2 "${SAMPLE}"/R2.p2.fq.gz \
		-S "${SAMPLE}"/lacto_rna.raw.sam \
		2>&1 | tee "${SAMPLE}"/lacto_rna.bowtie2.log
done

parallel -j 10 --keep-order '
	samtools view -F 256 {1}/lacto_rna.raw.sam \
		| awk '\''$6!="*"&&$7=="="&&$4==$8{print $1}'\'' \
		| perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl \
			>{1}/lacto_rna.out.list
	rm {1}/lacto_rna.raw.sam
	perl ~/NJU_seq/tool/delete_fastq.pl -n {1}/lacto_rna.out.list -i {1}/R1.p2.fq.gz -o {1}/R1.p3.fq.gz &
	perl ~/NJU_seq/tool/delete_fastq.pl -n {1}/lacto_rna.out.list -i {1}/R2.p2.fq.gz -o {1}/R2.p3.fq.gz &
	wait
' ::: NJU62{{32..35},{48..75}} NJUzp{3..18} Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3}

parallel -j 12 --keep-order '
	seqkit seq -g -m 30 ../NJU_seq_clean/{1}/R1.filter.fq.gz >{1}/R1.tmp.fq
	seqkit seq -g -m 30 ../NJU_seq_clean/{1}/R2.filter.fq.gz >{1}/R2.tmp.fq
	seqkit pair -1 {1}/R1.tmp.fq -2 {1}/R2.tmp.fq
	pigz <{1}/R1.tmp.paired.fq >{1}/R1.clean_long.fq.gz
	pigz <{1}/R2.tmp.paired.fq >{1}/R2.clean_long.fq.gz
	rm {}/R*.tmp*.fq
' ::: NJU62{{32..35},{48..75}} NJUzp{3..18} Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3} Mp_NC Mp_{1..3}

for SAMPLE in NJU62{{32..35},{48..75}} NJUzp{3..18} Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3} Mp_NC Mp_{1..3}; do
	mv ../NJU_seq_fast/"${SAMPLE}"/R1_clean.fq.gz "${SAMPLE}"/R1.origin.fq.gz
	mv ../NJU_seq_fast/"${SAMPLE}"/R2_clean.fq.gz "${SAMPLE}"/R2.origin.fq.gz
done

for SAMPLE in NJU62{{32..35},{48..75}} NJUzp{3..18} Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3} Mp_NC Mp_{1..3}; do
	mv ../NJU_seq_fast/"${SAMPLE}"/R1.p3.fq.gz "${SAMPLE}"/R1.filter.fq.gz
	mv ../NJU_seq_fast/"${SAMPLE}"/R2.p3.fq.gz "${SAMPLE}"/R2.filter.fq.gz
done

bsub -n 24 -q largemem -J kraken2 "
	for SAMPLE in Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3} NJU62{{32..35},{48..75}} NJUzp{3..18}; do
		/share/home/wangq/miniconda3/bin/kraken2 \\
			--threads 20 --use-names \\
			--gzip-compressed --paired \\
			--report \${SAMPLE}/class.clean_long.report \\
			--classified-out \${SAMPLE}/temp#.fq \\
			--db /share/home/wangq/wyf/kraken_plant \\
			\${SAMPLE}/R1.clean_long.fq.gz \\
			\${SAMPLE}/R2.clean_long.fq.gz \\
			--output \${SAMPLE}/class.tsv
	done
	rm ./*/temp* ./*/class*.tsv
"

for SAMPLE in NBG00{01..28} Osa_NC Osa_{1..3}; do
	mkdir -p ../NJU_seq_qc/"${SAMPLE}"
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAG \
		-O 6 -m 10 -e 0.1 --discard-untrimmed \
		-o "${SAMPLE}"/R1.p1.fq.gz -p "${SAMPLE}"/R2.p1.fq.gz \
		-j 16 "${SAMPLE}"/R1.fq.gz "${SAMPLE}"/R2.fq.gz \
		>"${SAMPLE}"/cutadapt.log 2>&1
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAG \
		-O 6 -m 30 -e 0.1 --discard-untrimmed \
		-o ../NJU_seq_qc/"${SAMPLE}"/R1.long.fq.gz \
		-p ../NJU_seq_qc/"${SAMPLE}"/R2.long.fq.gz \
		-j 16 "${SAMPLE}"/R1.fq.gz "${SAMPLE}"/R2.fq.gz \
		>../NJU_seq_qc/"${SAMPLE}"/cutadapt.long.log 2>&1
done

for SAMPLE in NBG00{01..28} Osa_NC Osa_{1..3}; do
	time bowtie2 --end-to-end -p 16 -a -t \
		--no-unal --no-mixed --no-discordant \
		-D 20 -R 3 -N 0 -L 10 -i S,1,0.75 \
		--maxins 120 --score-min L,-0.8,-0.8 \
		--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \
		--xeq -x ../index/ecoli_rna \
		-1 "${SAMPLE}"/R1_clean.fq.gz \
		-2 "${SAMPLE}"/R2_clean.fq.gz \
		-S "${SAMPLE}"/ecoli_rna.raw.sam \
		2>&1 | tee "${SAMPLE}"/ecoli_rna.bowtie2.log
done

parallel -j 16 --keep-order '
	samtools view -h -F 256 {1}/ecoli_rna.raw.sam \
		| samtools view -h -F 128 - | samtools view -F 16 - \
		| awk '\''$6!="*"&&$7=="="&&$4==$8{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" 10}'\'' \
		| perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl \
		| cut -f1 >{1}/ecoli_rna.out.list
	perl ~/NJU_seq/tool/delete_fastq.pl -n {1}/ecoli_rna.out.list -i {1}/R1.p1.fq.gz -o {1}/R1.p2.fq.gz &
	perl ~/NJU_seq/tool/delete_fastq.pl -n {1}/ecoli_rna.out.list -i {1}/R2.p1.fq.gz -o {1}/R2.p2.fq.gz &
	wait
' ::: NBG00{01..28} Osa_NC Osa_{1..3}

for SAMPLE in NJU62{{32..35},{48..75}} NJUzp{3..18} Col_SA_NC Col_SA_{1..3} Sly_NC Sly_{1..3}; do
	time bowtie2 --end-to-end -p 16 -k 20 -t \
		--no-unal --no-mixed --no-discordant \
		-D 20 -R 3 -N 0 -L 10 -i S,1,0.75 \
		--maxins 120 --score-min L,1.6,-0.8 \
		--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \
		--xeq -x ../index/lacto_rna \
		-1 "${SAMPLE}"/R1.p2.fq.gz \
		-2 "${SAMPLE}"/R2.p2.fq.gz \
		-S "${SAMPLE}"/lacto_rna.raw.sam \
		2>&1 | tee "${SAMPLE}"/lacto_rna.bowtie2.log
done
