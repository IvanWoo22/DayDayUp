for SAMPLE in Sly_{{1..3},NC}; do
	bsub -n 24 -J "$SAMPLE"_rrna "
	time bowtie2 -p 16 -a -t \
		--end-to-end -D 20 -R 3 \
		-N 0 -L 10 -i S,1,0.50 --np 0 \
		--xeq -x index/ath/ath_rrna_updated \
		-1 NJU_seq_reads/${SAMPLE}/R1_clean.fq.gz \
		-2 NJU_seq_reads/${SAMPLE}/R2_clean.fq.gz \
		-S NJU_seq_reads/${SAMPLE}/rrna.raw.sam \
		2>&1 | tee ${SAMPLE}/rrna.bowtie2.log
	"
done

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

for SAMPLE in NJU{1,2}{NC,1,2,3}; do
	bsub -n 24 -J "$SAMPLE"_rrna "
		bowtie2 -p 20 -a -t \\
			--no-unal --no-mixed --no-discordant \\
			--end-to-end -D 20 -R 3 \\
			-N 0 -L 10 -i S,1,0.75 --np 0 \\
			--xeq -x index/ath/ath_rrna_total \\
			-1 NJU_seq_reads/${SAMPLE}/R1.filter.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.filter.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/rrna.filter.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/rrna.filter.bowtie2.log
	"
done

bsub -n 24 -J sam_filter "
	parallel -j 20 --keep-order '
		samtools view -h -F 128 NJU_seq_output/{1}/rrna.filter.sam \\
			| samtools view -F 16 - \\
			| awk '\''\$6!=\"*\"&&\$7==\"=\"&&\$4==\$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \\
			| perl NJU_seq/rrna_analysis/matchquality_judge.pl \\
			| perl NJU_seq/rrna_analysis/multimatch_judge.pl \\
				>NJU_seq_output/{1}/rrna.filter.tmp
		cut -f1 NJU_seq_output/{1}/rrna.filter.tmp >NJU_seq_output/{1}/rrna.filter.list
	' ::: NJU{1,2}{NC,1,2,3}
	parallel -j 20 --keep-order '
		perl NJU_seq/tool/delete_fastq.pl \\
			-n NJU_seq_output/{1}/rrna.filter.list \\
			-i NJU_seq_reads/{1}/{2}.filter.fq.gz \\
			-o NJU_seq_reads/{1}/{2}.filter.mrna.fq.gz
	' ::: NJU{1,2}{NC,1,2,3} ::: R1 R2
"

for SAMPLE in NJU{1,2}{NC,1,2,3}; do
	bsub -n 24 -J "$SAMPLE"_mrna "
		bowtie2 --end-to-end -p 20 -a -t \\
			--no-unal --no-mixed --no-discordant \\
			 -D 20 -R 3 -i S,1,0.75 --maxins 120 \\
			-N 0 -L 10 --score-min L,0.8,-0.8 \\
			--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \\
			--xeq -x index/ath/ath_protein_coding \\
			-1 NJU_seq_reads/${SAMPLE}/R1.filter.mrna.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.filter.mrna.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/mrna.filter.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/mrna.filter.bowtie2.log
		"
done

bsub -n 24 -J sam_filter "
	parallel -j 20 --keep-order '
		samtools view -h -F 128 NJU_seq_output/{1}/mrna.filter.sam \\
			| samtools view -F 16 - \\
			| awk '\''\$6!=\"*\"&&\$7==\"=\"&&\$4==\$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \\
			| perl matchquality_judge.pl \\
			| perl NJU_seq/mrna_analysis/multimatch_judge.pl \\
				>NJU_seq_output/{1}/mrna.filter.tmp
	' ::: NJU{1,2}{NC,1,2,3}
"

for SAMPLE in NJU{1,2}{NC,1,2,3}; do
	bsub -n 1 -J "$SAMPLE"_mrna "
		perl NJU_seq/mrna_analysis/dedup.pl \\
			--refstr \"Parent=transcript:\" \\
			--info index/ath/exon.info \\
			<NJU_seq_output/${SAMPLE}/mrna.filter.tmp \\
			>NJU_seq_output/${SAMPLE}/mrna.dedup.filter.tmp
	"
done

bsub -n 24 -J dedup_count "
	parallel -j 10 --keep-order '
		perl NJU_seq/mrna_analysis/almostuniquematch.pl \\
			NJU_seq_reads/{1}/R1.filter.mrna.fq.gz \\
			NJU_seq_output/{1}/mrna.dedup.filter.tmp \\
			NJU_seq_output/{1}/mrna.almostuniquematch.filter.tmp
		perl NJU_seq/mrna_analysis/count.pl \\
			NJU_seq_output/{1}/mrna.almostuniquematch.filter.tmp \\
			| sort -k1,1 -k2,2n \\
			>NJU_seq_output/{1}/mrna.count.filter.tmp
	' ::: NJU{1,2}{NC,1,2,3}
"

bsub -n 24 -J sam_filter "
	parallel -j 20 --keep-order '
		perl NJU_seq/mrna_analysis/merge.pl \\
			--refstr \"Parent=transcript:\" \\
			--geneid \"AT\" \\
			--transid \"AT\" \\
			-i NJU_seq_output/{1}/mrna.count.filter.tmp \\
			--stdout \\
			<index/ath/exon.info \\
			| sort -k1,1 -k2,2n \\
			  >NJU_seq_output/{1}/mrna.filter.tsv
	' ::: NJU{1,2}{NC,1,2,3}
"

perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJU1{NC,1,2,3}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.NJU1.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJU2{NC,1,2,3}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.NJU2.tsv

for SPECIES in NJU{1,2}; do
	perl NJU_seq/presentation/signature_count.pl \
		NJU_seq_output/mrna.filter."${SPECIES}".tsv \
		NJU_seq_output/mrna.filter."${SPECIES}".signature.pdf
done

for SPECIES in NJU{1,2}; do
	perl NJU_seq/mrna_analysis/motif_nm.pl \
		index/ath/ath.fa <(sed '1d' NJU_seq_output/mrna.filter."${SPECIES}".tsv) 20 20 \
		>NJU_seq_output/mrna.filter."${SPECIES}".ctx.tsv
done

for SAMPLE in NJU62{52..55} Sly_{4,5,6}{TR,NC} Sly_{NC,{1..3}} NJUSl{{1..4},NC} Sly_8_1NC Sly_8_2NC Sly_8_1TR; do
	bsub -n 24 -J "$SAMPLE"_bac "
		bowtie2 --end-to-end -p 16 -k 20 -t \\
			--no-unal --no-mixed --no-discordant \\
			-D 20 -R 3 -N 0 -L 10 -i S,1,0.75 \\
			--maxins 120 --score-min L,-0.8,-0.8 \\
			--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 \\
			--ignore-quals --xeq -x index/bac_rna \\
			-1 NJU_seq_reads/${SAMPLE}/R1.origin.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.origin.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/bac_rna.raw.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/bac_rna.bowtie2.log
		"
done

bsub -n 24 -J sam_filter "
	parallel -j 20 --keep-order '
		samtools view -h -F 128 \\
			NJU_seq_output/{1}/bac_rna.raw.sam \\
			| samtools view -F 16 - \\
			| awk '\''\$6!=\"*\"&&\$7==\"=\"&&\$4==\$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \\
			| perl NJU_seq/rrna_analysis/multimatch_judge.pl \\
			| cut -f1 >NJU_seq_output/{1}/bac_rna.out.list
	' ::: NJU62{52..55} Sly_{4,5,6}{TR,NC} Sly_{NC,{1..3}} NJUSl{{1..4},NC} Sly_8_1NC Sly_8_2NC Sly_8_1TR
"

bsub -n 24 -J fastq_filter "
	parallel -j 20 --keep-order '
		perl NJU_seq/tool/delete_fastq.pl \\
			-n NJU_seq_output/{1}/bac_rna.out.list \\
			-i NJU_seq_reads/{1}/{2}.origin.fq.gz \\
			-o NJU_seq_reads/{1}/{2}.filter.fq.gz
	' ::: NJU62{52..55} Sly_{4,5,6}{TR,NC} Sly_{NC,{1..3}} NJUSl{{1..4},NC} Sly_8_1NC Sly_8_2NC Sly_8_1TR ::: R1 R2
"

for SAMPLE in NJU62{52..55} Sly_{4,5,6}{TR,NC} Sly_{NC,{1..3}} NJUSl{{1..4},NC} Sly_8_1NC Sly_8_2NC Sly_8_1TR; do
	bsub -n 24 -J "$SAMPLE"_rrna "
		bowtie2 -p 20 -a -t \\
			--no-unal --no-mixed --no-discordant \\
			--end-to-end -D 20 -R 3 \\
			-N 0 -L 10 -i S,1,0.75 --np 0 \\
			--xeq -x index/mct/mct_rrna_total \\
			-1 NJU_seq_reads/${SAMPLE}/R1.filter.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.filter.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/rrna.filter.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/rrna.filter.bowtie2.log
	"
done

parallel -j 10 --keep-order '
		samtools view -h -F 128 NJU_seq_fast/{1}/rrna.filter.sam \
			| samtools view -F 16 - \
			| awk '\''$6!="*"&&$7=="="&&$4==$8{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\'' \
			| perl ~/NJU_seq/rrna_analysis/matchquality_judge.pl \
			| perl ~/NJU_seq/rrna_analysis/multimatch_judge.pl \
				>NJU_seq_fast/{1}/rrna.filter.tmp
		cut -f1 NJU_seq_fast/{1}/rrna.filter.tmp >NJU_seq_fast/{1}/rrna.filter.list
		perl ~/NJU_seq/tool/delete_fastq.pl \
			-n NJU_seq_fast/{1}/rrna.filter.list \
			-i NJU_seq_clean/{1}/{2}.filter.fq.gz \
			-o NJU_seq_clean/{1}/{2}.filter.mrna.fq.gz
	' ::: Sly_{4,5,6}{TR,NC} ::: R1 R2

bsub -n 4 -J stats_length '
	for SAMPLE in NJU62{52..55} Sly_{4,5,6}{TR,NC} Sly_{NC,{1..3}} NJUSl{{1..4},NC} Sly_8_1NC Sly_8_2NC Sly_8_1TR; do
		for SOURCE in rrna.filter mrna.filter bac_rna.raw; do
			samtools sort -o NJU_seq_output/${SAMPLE}/${SOURCE}.bam \
				NJU_seq_output/${SAMPLE}/${SOURCE}.sam
			samtools index NJU_seq_output/${SAMPLE}/${SOURCE}.bam
			samtools stats -d NJU_seq_output/${SAMPLE}/${SOURCE}.bam \
				>NJU_seq_output/${SAMPLE}/${SOURCE}.sorted.stats.txt 2>/dev/null
		done
	done
'

for SAMPLE in NJU62{52..55} Sly_{4,5,6}{TR,NC} Sly_{NC,{1..3}} NJUSl{{1..4},NC} Sly_8_1NC Sly_8_2NC Sly_8_1TR; do
	for SOURCE in rrna.filter bac_rna.raw mrna.filter; do
		grep ^IS NJU_seq_output/"${SAMPLE}"/"${SOURCE}".sorted.stats.txt \
			>NJU_seq_output/"${SAMPLE}"."${SOURCE}".insert_size_distribution.txt
	done
done

for SAMPLE in Osa_{NC,{1..3}}; do
	bsub -n 24 -J "$SAMPLE"_bac "
		bowtie2 --end-to-end -p 16 -k 20 -t \\
			--no-unal --no-mixed --no-discordant \\
			-D 20 -R 3 -N 0 -L 10 -i S,1,0.75 \\
			--maxins 120 --score-min L,-0.8,-0.8 \\
			--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 \\
			--ignore-quals --xeq -x index/bac_rna \\
			-1 NJU_seq_reads/${SAMPLE}/R1.origin.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.origin.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/bac_rna.raw.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/bac_rna.bowtie2.log
	"
done

bsub -n 24 -J sam_filter "
	parallel -j 20 --keep-order '
		samtools view -h -F 128 \\
			NJU_seq_output/{1}/bac_rna.raw.sam \\
			| samtools view -F 16 - \\
			| awk '\''\$6!=\"*\"&&\$7==\"=\"&&\$4==\$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \\
			| perl NJU_seq/rrna_analysis/multimatch_judge.pl \\
			| cut -f1 >NJU_seq_output/{1}/bac_rna.out.list
	' ::: Osa_{NC,{1..3}}
"

bsub -n 24 -J fastq_filter "
	parallel -j 20 --keep-order '
		perl NJU_seq/tool/delete_fastq.pl \\
			-n NJU_seq_output/{1}/bac_rna.out.list \\
			-i NJU_seq_reads/{1}/{2}.origin.fq.gz \\
			-o NJU_seq_reads/{1}/{2}.filter.fq.gz
	' ::: Osa_{NC,{1..3}} ::: R1 R2
"

for SAMPLE in Osa_{NC,{1..3}}; do
	bsub -n 24 -J "$SAMPLE"_rrna "
		bowtie2 -p 20 -a -t \\
			--no-unal --no-mixed --no-discordant \\
			--end-to-end -D 20 -R 3 \\
			-N 0 -L 10 -i S,1,0.75 --np 0 \\
			--xeq -x index/osa/osa_rrna_total \\
			-1 NJU_seq_reads/${SAMPLE}/R1.filter.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.filter.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/rrna.filter.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/rrna.filter.bowtie2.log
	"
done

bsub -n 24 -J sam_filter "
	parallel -j 20 --keep-order '
		samtools view -h -F 128 NJU_seq_output/{1}/rrna.filter.sam \\
			| samtools view -F 16 - \\
			| awk '\''\$6!=\"*\"&&\$7==\"=\"&&\$4==\$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \\
			| perl NJU_seq/rrna_analysis/matchquality_judge.pl \\
			| perl NJU_seq/rrna_analysis/multimatch_judge.pl \\
				>NJU_seq_output/{1}/rrna.filter.tmp
		cut -f1 NJU_seq_output/{1}/rrna.filter.tmp >NJU_seq_output/{1}/rrna.filter.list
	' ::: Osa_{NC,{1..3}}
	parallel -j 20 --keep-order '
		perl NJU_seq/tool/delete_fastq.pl \\
			-n NJU_seq_output/{1}/rrna.filter.list \\
			-i NJU_seq_reads/{1}/{2}.filter.fq.gz \\
			-o NJU_seq_reads/{1}/{2}.filter.mrna.fq.gz
	' ::: Osa_{NC,{1..3}} ::: R1 R2
"

for SAMPLE in NJU62{{60..63},{72..75}} Osa_{NC,{1..3}}; do
	bsub -n 24 -J "$SAMPLE"_mrna "
		bowtie2 --end-to-end -p 20 -a -t \\
			--no-unal --no-mixed --no-discordant \\
			 -D 20 -R 3 -i S,1,0.75 --maxins 120 \\
			-N 0 -L 10 --score-min L,0.8,-0.8 \\
			--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \\
			--xeq -x index/osa/osa_protein_coding \\
			-1 NJU_seq_reads/${SAMPLE}/R1.filter.mrna.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.filter.mrna.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/mrna.filter.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/mrna.filter.bowtie2.log
	"
done

bsub -n 24 -J sam_filter "
	parallel -j 20 --keep-order '
		cut -f1 NJU_seq_output/{1}/mrna.filter.tmp >NJU_seq_output/{1}/mrna.filter.list
		perl NJU_seq/tool/delete_fastq.pl \\
			-n NJU_seq_output/{1}/mrna.filter.list \\
			-i NJU_seq_reads/{1}/{2}.filter.mrna.fq.gz \\
			-o NJU_seq_reads/{1}/{2}.others.fq.gz
	' ::: NJU62{52..55} Sly_{4,5,6}{TR,NC} Sly_{NC,{1..3}} NJUSl{{1..4},NC} Sly_8_1NC Sly_8_2NC Sly_8_1TR ::: R1 R2
"