barrnap --kingdom euk --threads 12 \
	--outseq data/ath_rrna.fa --reject 0.05 \
	nlr_eco/Col-PEK/chr.fasta \
	>data/ath_rrna.gff3
/home/ivan/miniconda3/envs/perl_old/bin/perl \
	rRNA_identify/rnammer/rnammer \
	-d -S euk -m lsu,ssu,tsu -multi \
	-gff ath.rrna.gff -xml ath.rrna.xml \
	-h ath.rrna.hmm -f ath.rrna.fasta \
	<nlr_eco/Col-PEK/chr.fasta
awk '{if($1!~"^#") print $1"\t"($4-1)"\t"$5"\t"$9"\t"$6"\t"$7}' ath.rrna.gff >ath.rnammer.bed
grep -v "^#" data/ath_rrna.gff3 | awk 'OFS="\t" {print $1, $4-1, $5, $9, $6, $7}' >ath.barrnap.bed
cat ath.rnammer.bed ath.barrnap.bed | sort -k1,1 -k2,2n >all_rDNA_sorted.bed
bedtools getfasta -fi nlr_eco/Col-PEK/chr.fasta -bed merged_rDNA.bed -fo merged_rDNA_sequences.fasta
cat rRNA_identify/ath_rrna.updated.fa merged_rDNA_sequences.fasta >rRNA_identify/ath_rrna.total.fa
bowtie2-build rRNA_identify/ath_rrna.total.fa rRNA_identify/index/ath_rrna_total

barrnap --kingdom euk --threads 12 \
	--outseq data/osa_rrna.fa --reject 0.05 \
	NJU_seq/osa.fa \
	>data/osa_rrna.gff3
/home/ivan/miniconda3/envs/perl_old/bin/perl \
	rRNA_identify/rnammer/rnammer \
	-d -S euk -m lsu,ssu,tsu -multi \
	-gff osa.rrna.gff -xml osa.rrna.xml \
	-h osa.rrna.hmm -f osa.rrna.fasta \
	<NJU_seq/osa.fa
awk '{if($1!~"^#") print $1"\t"($4-1)"\t"$5"\t"$9"\t"$6"\t"$7}' osa.rrna.gff >osa.rnammer.bed
grep -v "^#" data/osa_rrna.gff3 | awk 'OFS="\t" {print $1, $4-1, $5, $9, $6, $7}' >osa.barrnap.bed
cat osa.rnammer.bed osa.barrnap.bed | sort -k1,1 -k2,2n >all_rDNA_sorted.bed
bedtools merge -i all_rDNA_sorted.bed >merged_rDNA.bed
bedtools getfasta -fi NJU_seq/osa.fa -bed merged_rDNA.bed -fo merged_rDNA_sequences.fasta
cat rRNA_identify/osa_rrna.updated.fa merged_rDNA_sequences.fasta >rRNA_identify/osa_rrna.total.fa
bowtie2-build rRNA_identify/osa_rrna.total.fa rRNA_identify/index/osa_rrna_total

barrnap --kingdom euk --threads 12 \
	--outseq data/gma_rrna.fa --reject 0.05 \
	NJU_seq/gma.fa \
	>data/gma_rrna.gff3
/home/ivan/miniconda3/envs/perl_old/bin/perl \
	rRNA_identify/rnammer/rnammer \
	-d -S euk -m lsu,ssu,tsu -multi \
	-gff gma.rrna.gff -xml gma.rrna.xml \
	-h gma.rrna.hmm -f gma.rrna.fasta \
	<NJU_seq/gma.fa
awk '{if($1!~"^#") print $1"\t"($4-1)"\t"$5"\t"$9"\t"$6"\t"$7}' gma.rrna.gff >gma.rnammer.bed
grep -v "^#" data/gma_rrna.gff3 | awk 'OFS="\t" {print $1, $4-1, $5, $9, $6, $7}' >gma.barrnap.bed
cat gma.rnammer.bed gma.barrnap.bed | sort -k1,1 -k2,2n >all_rDNA_sorted.bed
bedtools merge -i all_rDNA_sorted.bed >merged_rDNA.bed
bedtools getfasta -fi NJU_seq/gma.fa -bed merged_rDNA.bed -fo merged_rDNA_sequences.fasta
cat rRNA_identify/gma_rrna.updated.fa merged_rDNA_sequences.fasta >rRNA_identify/gma_rrna.total.fa
bowtie2-build rRNA_identify/gma_rrna.total.fa rRNA_identify/index/gma_rrna_total

barrnap --kingdom euk --threads 12 \
	--outseq data/sly_rrna.fa --reject 0.05 \
	data/Microtom_genome/microTom.genome.fa \
	>data/sly_rrna.gff3
/home/ivan/miniconda3/envs/perl_old/bin/perl \
	rRNA_identify/rnammer/rnammer \
	-d -S euk -m lsu,ssu,tsu -multi \
	-gff sly.rrna.gff -xml sly.rrna.xml \
	-h sly.rrna.hmm -f sly.rrna.fasta \
	<data/Microtom_genome/microTom.genome.fa
awk '{if($1!~"^#") print $1"\t"($4-1)"\t"$5"\t"$9"\t"$6"\t"$7}' sly.rrna.gff >sly.rnammer.bed
grep -v "^#" data/sly_rrna.gff3 | awk 'OFS="\t" {print $1, $4-1, $5, $9, $6, $7}' >sly.barrnap.bed
cat sly.rnammer.bed sly.barrnap.bed | sort -k1,1 -k2,2n >all_rDNA_sorted.bed
bedtools merge -i all_rDNA_sorted.bed >merged_rDNA.bed
bedtools getfasta -fi NJU_seq/sly.fa -bed merged_rDNA.bed -fo merged_rDNA_sequences.fasta
cat rRNA_identify/sly_rrna.updated.fa merged_rDNA_sequences.fasta >rRNA_identify/mct_rrna.total.fa
bowtie2-build rRNA_identify/mct_rrna.total.fa rRNA_identify/index/mct_rrna_total

for SAMPLE in NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}; do
	bsub -n 24 -J "$SAMPLE"_rrna "
		bowtie2 -p 20 -a -t \\
			--no-unal --no-mixed --no-discordant \\
			--end-to-end -D 20 -R 3 \\
			-N 0 -L 10 -i S,1,0.75 --np 0 \\
			--xeq -x index/ath/ath_rrna_all \\
			-1 NJU_seq_reads/${SAMPLE}/R1.filter.mrna.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.filter.mrna.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/rrna.all.filter.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/rrna.all.filter.bowtie2.log
	"
done

bsub -n 24 -q largemem -J sam_filter "
	parallel -j 20 --keep-order '
		samtools view -h -F 128 NJU_seq_output/{1}/rrna.all.filter.sam \\
			| samtools view -F 16 - \\
			| awk '\''\$6!=\"*\"&&\$7==\"=\"&&\$4==\$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \\
			| perl NJU_seq/rrna_analysis/matchquality_judge.pl \\
			| perl NJU_seq/rrna_analysis/multimatch_judge.pl \\
				>NJU_seq_output/{1}/rrna.all.filter.tmp
		cut -f1 NJU_seq_output/{1}/rrna.all.filter.tmp >NJU_seq_output/{1}/rrna.all.filter.list
	' ::: NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}
"

parallel -j 20 --keep-order '
	perl NJU_seq/tool/delete_fastq.pl \\
		-n NJU_seq_output/{1}/rrna.filter.list \\
		-i NJU_seq_reads/{1}/{2}.filter.fq.gz \\
		-o NJU_seq_reads/{1}/{2}.filter.mrna.fq.gz
' ::: NJU62{{32..35},{48..75}} NBG00{01..28} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}} ::: R1 R2

for SAMPLE in NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}; do
	bsub -n 24 -J "$SAMPLE"_rrna "
		bowtie2 -p 20 -a -t --no-hd --end-to-end \\
			--no-unal --no-mixed --no-discordant \\
			-D 20 -R 3 -i S,1,0.75 --maxins 120 \\
			-N 0 -L 10 --score-min L,-2.4,-0.8 \\
			--rdg 8,4 --rfg 8,4 --np 12 --mp 12,12 --ignore-quals \\
			--xeq -x index/ath/ath_rrna_total \\
			-1 NJU_seq_reads/${SAMPLE}/R1.filter.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.filter.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/rrna.total.filter.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/rrna.total.filter.bowtie2.log
	"
done

bsub -n 24 -q largemem -J sam_filter "
	parallel -j 20 --keep-order '
		samtools view -h -F 128 NJU_seq_output/{1}/rrna.total.filter.sam \\
			| samtools view -F 16 - \\
			| awk '\''\$6!=\"*\"&&\$7==\"=\"&&\$4==\$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \\
			| perl NJU_seq/rrna_analysis/multimatch_judge.pl \\
				>NJU_seq_output/{1}/rrna.total.filter.tmp
		cut -f1 NJU_seq_output/{1}/rrna.total.filter.tmp \\
			| uniq | sort | uniq >NJU_seq_output/{1}/rrna.total.filter.list
	' ::: NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}
"

bsub -n 128 -q amd_milan -J sam_filter "
	parallel -j 60 --keep-order '
		perl NJU_seq/tool/delete_fastq.pl \\
			-n NJU_seq_output/{1}/rrna.total.filter.list \\
			-i NJU_seq_reads/{1}/{2}.filter.fq.gz \\
			-o NJU_seq_reads/{1}/{2}.filter.mrnas.fq.gz
	' ::: NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC} ::: R1 R2
"

for SAMPLE in NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}; do
	bsub -n 24 -J "$SAMPLE"_mrna "
		bowtie2 --end-to-end -p 20 -a -t \\
			--no-unal --no-mixed --no-discordant \\
			 -D 20 -R 3 -i S,1,0.75 --maxins 120 \\
			-N 0 -L 10 --score-min L,0.8,-0.8 \\
			--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \\
			--xeq -x index/ath/ath_protein_coding \\
			-1 NJU_seq_reads/${SAMPLE}/R1.filter.mrnas.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.filter.mrnas.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/mrnas.filter.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/mrnas.filter.bowtie2.log
		"
done

bsub -n 24 -q largemem -J sam_filter "
	parallel -j 20 --keep-order '
		samtools view -h -F 128 NJU_seq_output/{1}/mrnas.filter.sam \\
			| samtools view -F 16 - \\
			| awk '\''\$6!=\"*\"&&\$7==\"=\"&&\$4==\$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \\
			| perl matchquality_judge.pl \\
			| perl NJU_seq/mrna_analysis/multimatch_judge.pl \\
				>NJU_seq_output/{1}/mrnas.filter.tmp
	' ::: NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}
"

for SAMPLE in NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}; do
	bsub -n 1 -J "$SAMPLE"_mrna "
		perl NJU_seq/mrna_analysis/dedup.pl \\
			--refstr \"Parent=transcript:\" \\
			--info index/ath/exon.info \\
			<NJU_seq_output/${SAMPLE}/mrnas.filter.tmp \\
			>NJU_seq_output/${SAMPLE}/mrnas.dedup.filter.tmp
	"
done

for SAMPLE in NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}; do
	bsub -n 1 -J "$SAMPLE"_mrna "
		perl NJU_seq/mrna_analysis/almostuniquematch.pl \\
			NJU_seq_reads/${SAMPLE}/R1.filter.mrnas.fq.gz \\
			NJU_seq_output/${SAMPLE}/mrnas.dedup.filter.tmp \\
			NJU_seq_output/${SAMPLE}/mrnas.almostuniquematch.filter.tmp
		perl NJU_seq/mrna_analysis/count.pl \\
			NJU_seq_output/${SAMPLE}/mrnas.almostuniquematch.filter.tmp \\
			| sort -k1,1 -k2,2n \\
			>NJU_seq_output/${SAMPLE}/mrnas.count.filter.tmp
	"
done

for SAMPLE in NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}; do
	bsub -n 1 -J "$SAMPLE"_mrna "
		perl NJU_seq/mrna_analysis/merge.pl \\
			--refstr \"Parent=transcript:\" \\
			--geneid \"AT\" \\
			--transid \"AT\" \\
			-i NJU_seq_output/${SAMPLE}/mrnas.count.filter.tmp \\
			--stdout \\
			<index/ath/exon.info \\
			| sort -k1,1 -k2,2n \\
			  >NJU_seq_output/${SAMPLE}/mrnas.filter.tsv
	"
done

bsub -n 24 -q largemem -J sam_filter "
	parallel -j 20 --keep-order '
		perl new.dedup.pl --stdout \\
			--refstr \"Parent=transcript:\" \\
			--info index/ath/exon.info \\
			--geneid \"AT\" --transid \"AT\" \\
			<NJU_seq_output/{1}/mrnas.filter.tmp \\
			>NJU_seq_output/{1}/mrnas.dedup.filter.tmp
	' ::: NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}
"

bsub -n 24 -q largemem -J sam_filter "
	parallel -j 20 --keep-order '
		perl new.count.pl --dedup \\
			<NJU_seq_output/{1}/mrnas.dedup.filter.tmp \\
			| sort -k1,1 -k2,2n \\
			>NJU_seq_output/{1}/mrnas.dedup.filter.tsv
		perl new.count.pl \\
			<NJU_seq_output/{1}/mrnas.dedup.filter.tmp \\
			| sort -k1,1 -k2,2n \\
			>NJU_seq_output/{1}/mrnas.filter.tsv
	' ::: NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}
"

perl new.score.pl \
	NJU_seq_output/NJU62{32..35}/mrnas.filter.tsv \
	>NJU_seq_output/mrnas.filter.Ath_leaf_RF.tsv
perl new.score.pl \
	NJU_seq_output/NJU62{48..51}/mrnas.filter.tsv \
	>NJU_seq_output/mrnas.filter.Ath_leaf.tsv
perl new.score.pl \
	NJU_seq_output/NJUzp{06,{03..05}}/mrnas.filter.tsv \
	>NJU_seq_output/mrnas.filter.Ath_zp.tsv
perl new.score.pl \
	NJU_seq_output/NJUzp{10,{07..09}}/mrnas.filter.tsv \
	>NJU_seq_output/mrnas.filter.Ath_zpdc.tsv
perl new.score.pl \
	NJU_seq_output/Col_SA_{NC,{1..3}}/mrnas.filter.tsv \
	>NJU_seq_output/mrnas.filter.Ath_SA.tsv
perl new.score.pl \
	NJU_seq_output/NJU62{32..35}/mrnas.dedup.filter.tsv \
	>NJU_seq_output/mrnas.dedup.filter.Ath_leaf_RF.tsv
perl new.score.pl \
	NJU_seq_output/NJU62{48..51}/mrnas.dedup.filter.tsv \
	>NJU_seq_output/mrnas.dedup.filter.Ath_leaf.tsv
perl new.score.pl \
	NJU_seq_output/NJUzp{06,{03..05}}/mrnas.dedup.filter.tsv \
	>NJU_seq_output/mrnas.dedup.filter.Ath_zp.tsv
perl new.score.pl \
	NJU_seq_output/NJUzp{10,{07..09}}/mrnas.dedup.filter.tsv \
	>NJU_seq_output/mrnas.dedup.filter.Ath_zpdc.tsv
perl new.score.pl \
	NJU_seq_output/Col_SA_{NC,{1..3}}/mrnas.dedup.filter.tsv \
	>NJU_seq_output/mrnas.dedup.filter.Ath_SA.tsv

for SPECIES in Ath_leaf_RF Ath_leaf Ath_zp Ath_zpdc Ath_SA; do
	perl new.merge.pl \
		NJU_seq_output/mrnas.dedup.filter.${SPECIES}.tsv \
		NJU_seq_output/mrnas.filter.${SPECIES}.tsv \
		>NJU_seq_output/mrna.moreinfo.${SPECIES}.tsv
	cut -f 1-7,14,21-23,30-32,39-46 \
		NJU_seq_output/mrna.moreinfo.${SPECIES}.tsv \
		>NJU_seq_output/mrna.Nm.${SPECIES}.tsv
done

for SPECIES in Ath_leaf_RF Ath_leaf Ath_zp Ath_zpdc Ath_SA; do
	for base in A G C T; do
		perl NJU_seq/mrna_analysis/motif_nm.pl \
			index/ath/ath.fa \
			<(sed '1d' NJU_seq_output/mrna.merge.${SPECIES}.tsv | awk -v base=${base} '$NF!="dup"&&$5==base') \
			20 20 | awk -v OFS='' '
				{header=$1; header=header":" $2; header=header"|" $4$5$6; header=header" " $7; seq=$NF; print ">" header; print seq}
				' >NJU_seq_output/motif/mrna.dedup.${SPECIES}.20nt.${base}.fa
	done
done

for SPECIES in Ath_leaf_RF Ath_leaf Ath_zp Ath_zpdc Ath_SA; do
	awk '$' NJU_seq_output/mrna.moreinfo.${SPECIES}.tsv \
		>NJU_seq_output/mrna.highconf.${SPECIES}.tsv
done
