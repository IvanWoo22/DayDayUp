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
