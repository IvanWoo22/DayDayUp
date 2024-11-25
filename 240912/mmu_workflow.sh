for PREFIX in brain heart kidney lung muscle spleen; do
	pigz -dc data/mmu.gff3.gz | awk '$3=="gene"' \
		| perl ~/NJU_seq/tool/add_gene_name.pl \
			--id "gene_id=" --name "gene_name=" --col "8" \
			--file "NJU_seq/mrna.scored.Nm.mmu_${PREFIX}.tsv" \
			>NJU_seq/mrna.scored.Nm.mmu_${PREFIX}.tmp1
done

for PREFIX in brain heart kidney lung muscle spleen; do
	pigz -dc data/mmu_ens.gff3.gz | awk '$3=="gene"' \
		| perl ~/NJU_seq/tool/add_gene_name.pl \
			--id "gene_id=" --name "description=" --col "9" \
			--file "NJU_seq/mrna.scored.Nm.mmu_${PREFIX}.tmp1" \
			>NJU_seq/mrna.scored.Nm.mmu_${PREFIX}.tmp2
done

for PREFIX in brain heart kidney lung muscle spleen; do
	perl ~/NJU_seq/mrna_analysis/motif_nm.pl \
		data/GRCm38.primary_assembly.genome.fa \
		NJU_seq/mrna.scored.Nm.mmu_${PREFIX}.tmp2 30 30 \
		>NJU_seq/mrna.scored.Nm.mmu_${PREFIX}.tmp3
done

for PREFIX in brain heart kidney lung muscle spleen; do
	perl ~/NJU_seq/mrna_analysis/main_transcript_3.pl \
		data/mmu_basic_transcript_region.tsv \
		NJU_seq/mrna.scored.Nm.mmu_${PREFIX}.tmp3 \
		NJU_seq/mrna.scored.Nm.mmu_${PREFIX}.detailed.tsv
done

Total: 793929
Consistency: 744950
Proportion: 93.83%
Total: 3063928
Consistency: 3032450
Proportion: 98.97%

for PREFIX in NJU-BMS-05 NJU-BMS-06; do
	bsub -n 24 -J "count" "
		time bowtie2 -p 12 -a -t \
			--end-to-end -D 20 -R 3 \
			-N 0 -L 10 -i S,1,0.50 --np 0 \
			--xeq -x index/mmu_rrna \
			-1 data/${PREFIX}/R1_clean.fq.gz -2 data/${PREFIX}/R2_clean.fq.gz \
			-S output/${PREFIX}/rrna.raw.sam \
			2>&1 \
			| tee output/${PREFIX}/rrna.bowtie2.log
"
done

for PREFIX in NJU-BMS-05 NJU-BMS-06; do
	parallel <output/"${PREFIX}"/rrna.raw.sam --pipe --block 10M --no-run-if-empty --linebuffer --keep-order -j 4 '
    perl -nla -F"\t" -e '\''
      $F[5] ne qq(*) or next;
      $F[6] eq qq(=) or next;
      print join qq(\t), $F[0], $F[2], $F[3], $F[5], $F[9];
    '\'' |
    perl NJU_seq/rrna_analysis/matchquality_judge.pl |
    perl NJU_seq/rrna_analysis/multimatch_judge.pl
  ' \
		>output/"${PREFIX}"/rrna.out.tmp
done

for PREFIX in NJU-BMS-05 NJU-BMS-06; do
	time parallel -j 3 "
  perl NJU_seq/rrna_analysis/readend_count.pl \\
    NJU_seq/data/mmu_rrna/{}.fa output/${PREFIX}/rrna.out.tmp {} \\
    >output/${PREFIX}/rrna_{}.tsv
  " ::: 28s 18s 5-8s
done

time parallel -j 3 "
  perl NJU_seq/rrna_analysis/score.pl \\
    output/NJU-BMS-05/rrna_{}.tsv \\
    output/NJU-BMS-06/rrna_{}.tsv \\
      >output/NJU-BMS_rrna_{}_scored.tsv
  " ::: 28s 18s 5-8s

for PREFIX in NJU-BMS-05 NJU-BMS-06; do
	time bash NJU_seq/tool/extract_fastq.sh \
		output/"${PREFIX}"/rrna.out.tmp \
		data/"${PREFIX}"/R1_clean.fq.gz data/"${PREFIX}"/R1.mrna.fq.gz \
		data/"${PREFIX}"/R2_clean.fq.gz data/"${PREFIX}"/R2.mrna.fq.gz
done

for PREFIX in NJU-BMS-05 NJU-BMS-06; do
	bsub -n 24 -J "count" "
		time bowtie2 -p 20 -a -t \
			--end-to-end -D 20 -R 3 \
			-N 0 -L 10 --score-min C,0,0 \
			--xeq -x index/mmu_basic_protein_coding \
			-1 data/${PREFIX}/R1.mrna.fq.gz -2 data/${PREFIX}/R2.mrna.fq.gz \
			-S output/${PREFIX}/mrna.raw.sam \
			2>&1 \
			| tee output/${PREFIX}/mrna.bowtie2.log
			"
done

for PREFIX in NJU-BMS-05 NJU-BMS-06; do
	bsub -n 24 -J "count" "
		bash flow.sh ${PREFIX}
	"
done

PREFIX=${1}

parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 20 '
    awk '\''$6!="*"&&$7=="="{print $1 "\t" $3 "\t" $4 "\t" $6 "\t" $10}'\'' |
    perl NJU_seq/mrna_analysis/multimatch_judge.pl
  ' <output/${PREFIX}/mrna.raw.sam | perl NJU_seq/mrna_analysis/multimatch_judge.pl \
	>output/${PREFIX}/mrna.out.tmp

parallel --pipe --block 100M --no-run-if-empty --linebuffer --keep-order -j 20 '
    perl NJU_seq/mrna_analysis/dedup.pl --refstr "Parent=" --transid "ENSMUST" --info data/mmu_exon.info
  ' <output/"${PREFIX}"/mrna.out.tmp | perl NJU_seq/mrna_analysis/dedup.pl \
	--refstr "Parent=" \
	--transid "ENSMUST" \
	--info data/mmu_exon.info \
	>output/"${PREFIX}"/mrna.dedup.tmp

for PREFIX in NJU-BMS-05 NJU-BMS-06; do
	time bash NJU_seq/mrna_analysis/almostunique.sh \
		output/"${PREFIX}"/mrna.dedup.tmp \
		data/"${PREFIX}"/R1.mrna.fq.gz \
		output/"${PREFIX}" \
		output/"${PREFIX}"/mrna.almostunique.tmp
done

for PREFIX in NJU-BMS-05 NJU-BMS-06; do
	time perl NJU_seq/mrna_analysis/count.pl \
		output/"${PREFIX}"/mrna.almostunique.tmp \
		>output/"${PREFIX}"/mrna.count.tmp

	time gzip -dcf data/mmu.gff3.gz \
		| awk '$3=="exon" {print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' \
		| perl NJU_seq/mrna_analysis/merge.pl \
			--refstr "Parent=" \
			--geneid "ENSMUSG" \
			--transid "ENSMUST" \
			-i output/"${PREFIX}"/mrna.count.tmp \
			-o output/"${PREFIX}"/mrna.tsv
done
