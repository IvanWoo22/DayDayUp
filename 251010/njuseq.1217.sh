for SAMPLE in Mp_{{1..3},NC} Sly_{{1..3},NC}; do
	mkdir "${SAMPLE}"
	ln -sf /home/ivan/cater_data/nlr_data/njuseq/Mp_Sly_1217/"${SAMPLE}"/*_1.fq.gz /home/ivan/fat/NJU_seq/"${SAMPLE}"/R1.fq.gz
	ln -sf /home/ivan/cater_data/nlr_data/njuseq/Mp_Sly_1217/"${SAMPLE}"/*_2.fq.gz /home/ivan/fat/NJU_seq/"${SAMPLE}"/R2.fq.gz
done

for SAMPLE in Mp_{{1..3},NC} Sly_{{1..3},NC}; do
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAG \
		-O 6 -m 10 -e 0.1 --discard-untrimmed -o "${SAMPLE}"/R1_clean.fq.gz -p "${SAMPLE}"/R2_clean.fq.gz \
		-j 16 "${SAMPLE}"/R1.fq.gz "${SAMPLE}"/R2.fq.gz \
		>"${SAMPLE}"/cutadapt.log 2>&1

	perl ~/NJU_seq/quality_control/pe_consistency.pl \
		"${SAMPLE}"/R1_clean.fq.gz \
		"${SAMPLE}"/R2_clean.fq.gz \
		"${SAMPLE}"/merged_temp.fq.gz
done

for SAMPLE in Mp_{{1..3},NC} Sly_{{1..3},NC}; do
	kraken2 --threads 16 --use-names \
		--gzip-compressed --paired \
		--report "${SAMPLE}"/class.report \
		--classified-out "${SAMPLE}"/temp#.fq \
		--db /home/ivan/fat/data/kraken_plant/ \
		"${SAMPLE}"/R1_clean.fq.gz \
		"${SAMPLE}"/R2_clean.fq.gz \
		--output "${SAMPLE}"/class.tsv
done

cat ~/NJU_seq/data/sly_rrna/* >sly_rrna.fa
bowtie2-build sly_rrna.fa ../index/sly_rrna

for SAMPLE in Sly_{{1..3},NC}; do
	time bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x ../index/sly_rrna \
		-1 "${SAMPLE}"/R1_clean.fq.gz \
		-2 "${SAMPLE}"/R2_clean.fq.gz \
		-S "${SAMPLE}"/rrna.raw.sam \
		2>&1 | tee "${SAMPLE}"/rrna.bowtie2.log

	time pigz -p 16 "${SAMPLE}"/rrna.raw.sam
	time pigz -dcf "${SAMPLE}"/rrna.raw.sam.gz \
		| parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j 12 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\'' |
    perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl' \
			>"${SAMPLE}"/rrna.out.tmp

	time bash ~/NJU_seq/tool/extract_fastq.sh \
		"${SAMPLE}"/rrna.out.tmp \
		"${SAMPLE}"/R1_clean.fq.gz "${SAMPLE}"/R1.mrna.fq.gz \
		"${SAMPLE}"/R2_clean.fq.gz "${SAMPLE}"/R2.mrna.fq.gz
done

for SAMPLE in Sly_{{1..3},NC}; do
	for RNA in 25s 17s 5-8s; do
		perl ~/NJU_seq/rrna_analysis/readend_count.pl \
			~/NJU_seq/data/sly_rrna/${RNA}.fa \
			"${SAMPLE}"/rrna.out.tmp ${RNA} \
			>"${SAMPLE}"/rrna_${RNA}.tsv
	done
done

time parallel -j 3 "
	perl ~/NJU_seq/rrna_analysis/score.pl \
		Sly_NC/rrna_{}.tsv \
		Sly_1/rrna_{}.tsv \
		Sly_2/rrna_{}.tsv \
		Sly_3/rrna_{}.tsv \
		>rrna.scored.Nm.{}.Sly.tsv
" ::: 25s 17s 5-8s

for SAMPLE in Sly_{{1..3},NC}; do
	time bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 --score-min C,0,0 \
		--xeq -x ../index/sly_protein_coding \
		-1 "${SAMPLE}"/R1.mrna.fq.gz -2 "${SAMPLE}"/R2.mrna.fq.gz \
		-S "${SAMPLE}"/mrna.raw.sam \
		2>&1 | tee "${SAMPLE}"/mrna.bowtie2.log

	parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 12 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\'' |
    perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl
  ' <"${SAMPLE}"/mrna.raw.sam \
		| perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl \
			>"${SAMPLE}"/mrna.out.tmp
done

pigz -dc /home/ivan/fat/data/Slycopersicum_514_ITAG3.2.gene_exons.gff.gz \
	| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
	| sort -k1,1 -k2,2n >sly_exon.info

for SAMPLE in Sly_{{1..3},NC}; do
	parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 12 '
    /home/linuxbrew/.linuxbrew/bin/perl ~/NJU_seq/mrna_analysis/dedup.pl --refstr "Parent=" --transid "Solyc" --info sly_exon.info
  ' <"${SAMPLE}"/mrna.out.tmp \
		| /home/linuxbrew/.linuxbrew/bin/perl ~/NJU_seq/mrna_analysis/dedup.pl \
			--refstr "Parent=" \
			--transid "Solyc" \
			--info sly_exon.info \
			>"${SAMPLE}"/mrna.dedup.tmp

	time perl almostuniquematchneo.pl \
		"${SAMPLE}"/R1.mrna.fq.gz \
		"${SAMPLE}"/mrna.dedup.tmp \
		"${SAMPLE}"/mrna.almostuniquematchneo.tmp

	time perl countneo.pl \
		"${SAMPLE}"/mrna.almostuniquematchneo.tmp \
		| sort -k1,1 -k2,2n >"${SAMPLE}"/mrna.count.tmp

	perl ~/NJU_seq/mrna_analysis/merge.pl \
			--refstr "Parent=" \
			--geneid "Solyc" \
			--transid "Solyc" \
			-i "${SAMPLE}"/mrna.count.tmp \
			--stdout \
			<sly_exon.info \
		| sort -k1,1 -k2,2n >"${SAMPLE}"/mrna.tsv
done

for SAMPLE in Sly_{{1..3},NC}; do
	time perl ~/NJU_seq/mrna_analysis/almostuniquematch.pl \
		"${SAMPLE}"/R1.mrna.fq.gz \
		"${SAMPLE}"/mrna.dedup.tmp \
		"${SAMPLE}"/mrna.almostuniquematch.tmp

	time perl ~/NJU_seq/mrna_analysis/count.pl \
		"${SAMPLE}"/mrna.almostuniquematch.tmp \
		| sort -k1,1 -k2,2n >"${SAMPLE}"/mrna.count.uniq.tmp

	perl ~/NJU_seq/mrna_analysis/merge.pl \
			--refstr "Parent=" \
			--geneid "Solyc" \
			--transid "Solyc" \
			-i "${SAMPLE}"/mrna.count.uniq.tmp \
			--stdout \
			<sly_exon.info \
		| sort -k1,1 -k2,2n >"${SAMPLE}"/mrna.uniq.tsv
done

perl ~/NJU_seq/mrna_analysis/score.pl \
	Sly_NC/mrna.uniq.tsv Sly_1/mrna.uniq.tsv Sly_2/mrna.uniq.tsv Sly_3/mrna.uniq.tsv \
	>mrna.scored.Nm.Sly.uniq.tsv

for SAMPLE in 1 2 3; do
	perl ~/NJU_seq/mrna_analysis/score.pl \
		Sly_NC/mrna.uniq.tsv Sly_${SAMPLE}/mrna.uniq.tsv \
		>mrna.scored.Nm.Sly_${SAMPLE}.uniq.tsv
done

perl ~/NJU_seq/mrna_analysis/motif_nm.pl \
	sly.fa mrna.scored.Nm.Sly.tsv 20 20 \
	>mrna.scored.Nm.Sly.ctx.tmp
perl ~/NJU_seq/mrna_analysis/motif_nm.pl \
	sly.fa mrna.scored.Nm.Sly.ctx.tmp 50 50 \
	>mrna.scored.Nm.Sly.ctx.tsv

perl ~/NJU_seq/mrna_analysis/motif_nm.pl \
	sly.fa mrna.scored.Nm.Sly.uniq.tsv 20 20 \
	>mrna.scored.Nm.Sly.uniq.ctx.tmp
perl ~/NJU_seq/mrna_analysis/motif_nm.pl \
	sly.fa mrna.scored.Nm.Sly.uniq.ctx.tmp 50 50 \
	>mrna.scored.Nm.Sly.uniq.ctx.tsv

rm ./*.ctx.tmp

perl ~/NJU_seq/presentation/signature_count.pl \
	mrna.scored.Nm.SA.tsv \
	mrna.scored.Nm.SA.signature.pdf

for SAMPLE in Col_SA_{{1..3},NC}; do
	sed -n '32,35p' "${SAMPLE}"/cutadapt.log \
		| awk -vs="${SAMPLE}" \
			'BEGIN{print "base\t"s}
				{gsub(/^[ \t]+|[ \t]+$/,"");
				gsub(/[ \t]+/," ");
				if($0~/:/){split($0,a,":");
				gsub(/^[ \t]+|[ \t]+$/,"",a[1]);
				gsub(/^[ \t]+|[ \t]+$/,"",a[2]);
				print a[1]"\t"a[2]}}' \
		| datamash transpose \
			>"${SAMPLE}"/cutadapt.tsv
done
tsv-append -H ./*/cutadapt.tsv >cutadapt.tsv

for SAMPLE in Mp_{{1..3},NC}; do
	time bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x ../index/mpo_rrna \
		-1 "${SAMPLE}"/R1_clean.fq.gz \
		-2 "${SAMPLE}"/R2_clean.fq.gz \
		-S "${SAMPLE}"/rrna.raw.sam \
		2>&1 | tee "${SAMPLE}"/rrna.bowtie2.log

	time pigz -p 16 "${SAMPLE}"/rrna.raw.sam
	time pigz -dcf "${SAMPLE}"/rrna.raw.sam.gz \
		| parallel --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j 12 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\'' |
    perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl' \
			>"${SAMPLE}"/rrna.out.tmp

	time bash ~/NJU_seq/tool/extract_fastq.sh \
		"${SAMPLE}"/rrna.out.tmp \
		"${SAMPLE}"/R1_clean.fq.gz "${SAMPLE}"/R1.mrna.fq.gz \
		"${SAMPLE}"/R2_clean.fq.gz "${SAMPLE}"/R2.mrna.fq.gz
done

for SAMPLE in Mp_{{1..3},NC}; do
	for RNA in 25s 18s 5-8s; do
		perl ~/NJU_seq/rrna_analysis/readend_count.pl \
			~/NJU_seq/data/mpo_rrna/${RNA}.fa \
			"${SAMPLE}"/rrna.out.tmp ${RNA} \
			>"${SAMPLE}"/rrna_${RNA}.tsv
	done
done

for SAMPLE in Mp_{{1..3},NC}; do
	time bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 --score-min C,0,0 \
		--xeq -x ../index/mpo_protein_coding \
		-1 "${SAMPLE}"/R1.mrna.fq.gz -2 "${SAMPLE}"/R2.mrna.fq.gz \
		-S "${SAMPLE}"/mrna.raw.sam \
		2>&1 | tee "${SAMPLE}"/mrna.bowtie2.log

	parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 12 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\'' |
    perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl
  ' <"${SAMPLE}"/mrna.raw.sam \
		| perl ~/NJU_seq/mrna_analysis/multimatch_judge.pl \
			>"${SAMPLE}"/mrna.out.tmp
done

awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' MpTak1_v7.1.gff \
	| grep -v ".pre" | sort -k1,1 -k2,2n >mpo_exon.info

for SAMPLE in Mp_{{1..3},NC}; do
	parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 12 '
    /home/linuxbrew/.linuxbrew/bin/perl ~/NJU_seq/mrna_analysis/dedup.pl --refstr "Parent=" --transid "Mp" --info mpo_exon.info
  ' <"${SAMPLE}"/mrna.out.tmp \
		| /home/linuxbrew/.linuxbrew/bin/perl ~/NJU_seq/mrna_analysis/dedup.pl \
			--refstr "Parent=" \
			--transid "Mp" \
			--info mpo_exon.info \
			>"${SAMPLE}"/mrna.dedup.tmp

	time perl almostuniquematchneo.pl \
		"${SAMPLE}"/R1.mrna.fq.gz \
		"${SAMPLE}"/mrna.dedup.tmp \
		"${SAMPLE}"/mrna.almostuniquematchneo.tmp

	time perl countneo.pl \
		"${SAMPLE}"/mrna.almostuniquematchneo.tmp \
		| sort -k1,1 -k2,2n >"${SAMPLE}"/mrna.count.tmp

	awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' MpTak1_v7.1.gff \
		| grep -v ".pre" \
		| perl ~/NJU_seq/mrna_analysis/merge.pl \
			--refstr "Parent=" \
			--geneid "Mp" \
			--transid "Mp" \
			-i "${SAMPLE}"/mrna.count.tmp \
			--stdout \
		| sort -k1,1 -k2,2n >"${SAMPLE}"/mrna.tsv
done

perl ~/NJU_seq/mrna_analysis/score.pl \
	Mp_NC/mrna.tsv Mp_1/mrna.tsv Mp_2/mrna.tsv Mp_3/mrna.tsv \
	>mrna.scored.Nm.Mp.tsv

for SAMPLE in 1 2 3; do
	perl ~/NJU_seq/mrna_analysis/score.pl \
		Mp_NC/mrna.tsv Mp_${SAMPLE}/mrna.tsv \
		>mrna.scored.Nm.Mp_${SAMPLE}.tsv
done

for SAMPLE in Mp_{{1..3},NC}; do
	time perl ~/NJU_seq/mrna_analysis/almostuniquematch.pl \
		"${SAMPLE}"/R1.mrna.fq.gz \
		"${SAMPLE}"/mrna.dedup.tmp \
		"${SAMPLE}"/mrna.almostuniquematch.tmp

	time perl ~/NJU_seq/mrna_analysis/count.pl \
		"${SAMPLE}"/mrna.almostuniquematch.tmp \
		| sort -k1,1 -k2,2n >"${SAMPLE}"/mrna.count.uniq.tmp

	awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' MpTak1_v7.1.gff \
		| grep -v ".pre" \
		| perl ~/NJU_seq/mrna_analysis/merge.pl \
			--refstr "Parent=" \
			--geneid "Mp" \
			--transid "Mp" \
			-i "${SAMPLE}"/mrna.count.uniq.tmp \
			--stdout \
		| sort -k1,1 -k2,2n >"${SAMPLE}"/mrna.uniq.tsv
done

perl ~/NJU_seq/mrna_analysis/score.pl \
	Mp_NC/mrna.uniq.tsv Mp_1/mrna.uniq.tsv Mp_2/mrna.uniq.tsv Mp_3/mrna.uniq.tsv \
	>mrna.scored.Nm.Mp.uniq.tsv

for SAMPLE in 1 2 3; do
	perl ~/NJU_seq/mrna_analysis/score.pl \
		Mp_NC/mrna.uniq.tsv Mp_${SAMPLE}/mrna.uniq.tsv \
		>mrna.scored.Nm.Mp_${SAMPLE}.uniq.tsv
done

for SPECIES in Mp.uniq Mp Sly.uniq Sly; do
	perl ~/NJU_seq/presentation/signature_count.pl \
		mrna.scored.Nm.${SPECIES}.tsv \
		mrna.scored.Nm.${SPECIES}.signature.pdf
done

for SAMPLE in Col_SA_{{1..3},NC}; do
	sed -n '32,35p' "${SAMPLE}"/cutadapt.log \
		| awk -vs="${SAMPLE}" \
			'BEGIN{print "base\t"s}
				{gsub(/^[ \t]+|[ \t]+$/,"");
				gsub(/[ \t]+/," ");
				if($0~/:/){split($0,a,":");
				gsub(/^[ \t]+|[ \t]+$/,"",a[1]);
				gsub(/^[ \t]+|[ \t]+$/,"",a[2]);
				print a[1]"\t"a[2]}}' \
		| datamash transpose \
			>"${SAMPLE}"/cutadapt.tsv
done
tsv-append -H ./*/cutadapt.tsv >cutadapt.tsv
