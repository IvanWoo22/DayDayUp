for SAMPLE in Col_SA_NC Col_SA_1 Col_SA_2 Col_SA_3; do
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAG \
		-O 6 -m 10 -e 0.1 --discard-untrimmed -o ${SAMPLE}/R1_clean.fq.gz -p ${SAMPLE}/R2_clean.fq.gz \
		-j 16 ${SAMPLE}/R1.fq.gz ${SAMPLE}/R2.fq.gz \
		>${SAMPLE}/cutadapt.log 2>&1

	perl ~/NJU_seq/quality_control/pe_consistency.pl \
		${SAMPLE}/R1_clean.fq.gz \
		${SAMPLE}/R2_clean.fq.gz \
		${SAMPLE}/merged_temp.fq.gz

	time bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x ../index/ath_rrna \
		-1 ${SAMPLE}/R1_clean.fq.gz \
		-2 ${SAMPLE}/R2_clean.fq.gz \
		-S ${SAMPLE}/rrna.raw.sam \
		2>&1 | tee ${SAMPLE}/rrna.bowtie2.log

	time pigz -p 16 ${SAMPLE}/rrna.raw.sam
	time pigz -dcf ${SAMPLE}/rrna.raw.sam.gz \
		| parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j 12 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\'' |
    perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl' \
			>${SAMPLE}/rrna.out.tmp

	time bash ~/NJU_seq/tool/extract_fastq.sh \
		${SAMPLE}/rrna.out.tmp \
		${SAMPLE}/R1_clean.fq.gz ${SAMPLE}/R1.mrna.fq.gz \
		${SAMPLE}/R2_clean.fq.gz ${SAMPLE}/R2.mrna.fq.gz

	time bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 --score-min C,0,0 \
		--xeq -x ../index/ath_protein_coding \
		-1 ${SAMPLE}/R1.mrna.fq.gz -2 ${SAMPLE}/R2.mrna.fq.gz \
		-S ${SAMPLE}/mrna.raw.sam \
		2>&1 | tee ${SAMPLE}/mrna.bowtie2.log

	parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 12 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\'' |
    perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl
  ' <${SAMPLE}/mrna.raw.sam \
		| perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl \
			>${SAMPLE}/mrna.out.tmp

	parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 12 '
    /home/linuxbrew/.linuxbrew/bin/perl ~/NJU_seq/mrna_analysis/dedup.pl --refstr "Parent=transcript:" --transid "AT" --info ../data/ath_exon.info
  ' <${SAMPLE}/mrna.out.tmp \
		| /home/linuxbrew/.linuxbrew/bin/perl ~/NJU_seq/mrna_analysis/dedup.pl \
			--refstr "Parent=transcript:" \
			--transid "AT" \
			--info ../data/ath_exon.info \
			>${SAMPLE}/mrna.dedup.tmp

	time perl almostuniquematchneo.pl \
		${SAMPLE}/R1.mrna.fq.gz \
		${SAMPLE}/mrna.dedup.tmp \
		${SAMPLE}/mrna.almostuniquematchneo.tmp

	time perl countneo.pl \
		${SAMPLE}/mrna.almostuniquematchneo.tmp \
		>${SAMPLE}/mrna.count.tmp

	time gzip -dcf ../data/ath.gff3.gz \
		| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
		| perl ~/NJU_seq/mrna_analysis/merge.pl \
			--refstr "Parent=transcript:" \
			--geneid "AT" \
			--transid "AT" \
			-i ${SAMPLE}/mrna.count.tmp \
			-o ${SAMPLE}/mrna.tsv
	pigz ${SAMPLE}/mrna.raw.sam
done

for SAMPLE in Col_SA_NC Col_SA_1 Col_SA_2 Col_SA_3; do
	for RNA in 25s 18s 5-8s; do
		perl ~/NJU_seq/rrna_analysis/readend_count.pl \
			~/NJU_seq/data/ath_rrna/${RNA}.fa \
			${SAMPLE}/rrna.out.tmp ${RNA} \
			>${SAMPLE}/rrna_${RNA}.tsv
	done
done

time parallel -j 3 "
	perl ~/NJU_seq/rrna_analysis/score.pl \
		Col_SA_NC/rrna_{}.tsv \
		Col_SA_1/rrna_{}.tsv \
		Col_SA_2/rrna_{}.tsv \
		Col_SA_3/rrna_{}.tsv \
		>rrna.scored.Nm.{}.SA.tsv
" ::: 25s 18s 5-8s

perl /home/ivan/cater_data/ivan/fat/mscore.modified.pl \
	Col_SA_NC/mrna.tsv Col_SA_1/mrna.tsv Col_SA_2/mrna.tsv Col_SA_3/mrna.tsv \
	>mrna.scored.Nm.SA.tsv

for SAMPLE in 1 2 3; do
	perl /home/ivan/cater_data/ivan/fat/mscore.modified.pl \
		Col_SA_NC/mrna.tsv Col_SA_${SAMPLE}/mrna.tsv \
		>mrna.scored.Nm.SA_${SAMPLE}.tsv
done

perl ~/NJU_seq/mrna_analysis/motif_nm.pl \
	ath.fa mrna.scored.Nm.SA.tsv 30 30 \
	>mrna.scored.Nm.SA.ctx.tsv

perl ~/NJU_seq/presentation/signature_count.pl \
	mrna.scored.Nm.SA.tsv \
	mrna.scored.Nm.SA.signature.pdf
