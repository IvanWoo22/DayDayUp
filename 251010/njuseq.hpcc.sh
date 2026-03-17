for SAMPLE in NJU62{{32..35},{48..75}} NBG00{01..28} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}}; do
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
	' ::: NJU62{32..35} NJU62{48..75} NBG00{01..28} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}} NJUOs{NC,{1..4}}
"

bsub -n 128 -q amd_milan -J fastq_filter1 "
	parallel -j 100 --keep-order '
		perl NJU_seq/tool/delete_fastq.pl \\
			-n NJU_seq_output/{1}/bac_rna.out.list \\
			-i NJU_seq_reads/{1}/{2}.origin.fq.gz \\
			-o NJU_seq_reads/{1}/{2}.filter.fq.gz
	' ::: NJU62{48..75} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}} ::: R1 R2
"

bsub -n 128 -q amd_milan -J fastq_filter2 "
	parallel -j 100 --keep-order '
		perl NJU_seq/tool/delete_fastq.pl \\
			-n NJU_seq_output/{1}/bac_rna.out.list \\
			-i NJU_seq_reads/{1}/{2}.origin.fq.gz \\
			-o NJU_seq_reads/{1}/{2}.filter.fq.gz
	' ::: NBG00{01..28} NJU62{32..35} NJUOs{NC,{1..4}} ::: R1 R2
"

parallel -j 12 --keep-order '
	seqkit seq -g -m 30 ../NJU_seq_reads/{1}/R1.filter.fq.gz >{1}/R1.tmp.fq
	seqkit seq -g -m 30 ../NJU_seq_reads/{1}/R2.filter.fq.gz >{1}/R2.tmp.fq
	seqkit pair -1 {1}/R1.tmp.fq -2 {1}/R2.tmp.fq
	pigz <{1}/R1.tmp.paired.fq >{1}/R1.clean_long.fq.gz
	pigz <{1}/R2.tmp.paired.fq >{1}/R2.clean_long.fq.gz
	rm {1}/R*.tmp*.fq
' ::: NJU62{48..75} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}} NBG00{01..28} NJU62{32..35} NJUOs{NC,{1..4}}

for SAMPLE in NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}; do
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

for SAMPLE in NJU62{{56..59},{68..71}} NJUzp{11..18}; do
	bsub -n 24 -J "$SAMPLE"_rrna "
		bowtie2 -p 20 -a -t \\
			--no-unal --no-mixed --no-discordant \\
			--end-to-end -D 20 -R 3 \\
			-N 0 -L 10 -i S,1,0.75 --np 0 \\
			--xeq -x index/gma/gma_rrna_total \\
			-1 NJU_seq_reads/${SAMPLE}/R1.filter.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.filter.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/rrna.filter.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/rrna.filter.bowtie2.log
	"
done

for SAMPLE in NJU62{{60..63},{72..75}} Osa_{NC,{1..3}} NJUOs{{1..4},NC}; do
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

#for SAMPLE in NJU62{{52..55},{64..67}}; do
#	bsub -n 24 -J "$SAMPLE"_rrna "
#		bowtie2 -p 20 -a -t \\
#			--no-unal --no-mixed --no-discordant \\
#			--end-to-end -D 20 -R 3 \\
#			-N 0 -L 10 -i S,1,0.75 --np 0 \\
#			--xeq -x index/sly/sly_rrna_updated \\
#			-1 NJU_seq_reads/${SAMPLE}/R1.filter.fq.gz \\
#			-2 NJU_seq_reads/${SAMPLE}/R2.filter.fq.gz \\
#			-S NJU_seq_output/${SAMPLE}/rrna.filter.sam \\
#			2>&1 | tee NJU_seq_output/${SAMPLE}/rrna.filter.bowtie2.log
#	"
#done

for SAMPLE in NJU62{{52..55},{64..67}} Sly_{{1..3},NC}; do
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

bsub -n 128 -q amd_milan -J sam_filter "
	parallel -j 100 --keep-order '
		samtools view -h -F 128 NJU_seq_output/{1}/rrna.filter.sam \\
			| samtools view -F 16 - \\
			| awk '\''\$6!=\"*\"&&\$7==\"=\"&&\$4==\$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \\
			| perl NJU_seq/rrna_analysis/matchquality_judge.pl \\
			| perl NJU_seq/rrna_analysis/multimatch_judge.pl \\
				>NJU_seq_output/{1}/rrna.filter.tmp
		cut -f1 NJU_seq_output/{1}/rrna.filter.tmp >NJU_seq_output/{1}/rrna.filter.list
	' ::: NJU62{{32..35},{48..75}} NBG00{01..28} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}} NJUT{1,2,3,NC} NJUZ{1,2,3,NC} NJUSl{{1..4},NC} NJUOs{{1..4},NC}
	parallel -j 100 --keep-order '
		perl NJU_seq/tool/delete_fastq.pl \\
			-n NJU_seq_output/{1}/rrna.filter.list \\
			-i NJU_seq_reads/{1}/{2}.filter.fq.gz \\
			-o NJU_seq_reads/{1}/{2}.filter.mrna.fq.gz
	' ::: NJU62{{32..35},{48..75}} NBG00{01..28} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}} NJUT{1,2,3,NC} NJUZ{1,2,3,NC} NJUSl{{1..4},NC} NJUOs{{1..4},NC} ::: R1 R2
"

bsub -n 24 -J rrna_score "
	parallel -j 20 --keep-order '
		perl NJU_seq/rrna_analysis/readend_count.pl \
			index/ath/{2}.fa \
			NJU_seq_output/{1}/rrna.filter.tmp {2} \
			>NJU_seq_output/{1}/rrna.filter.{2}.tsv
	' ::: NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC} \\
	::: 25s 18s 5-8s 5s
	parallel -j 20 --keep-order '
		perl NJU_seq/rrna_analysis/readend_count.pl \
			index/gma/{2}.fa \
			NJU_seq_output/{1}/rrna.filter.tmp {2} \
			>NJU_seq_output/{1}/rrna.filter.{2}.tsv
	' ::: NJU62{{56..59},{68..71}} NJUzp{11..18} NBG00{01..28} \\
	::: 25s 18s 5-8s 5s
	parallel -j 20 --keep-order '
		perl NJU_seq/rrna_analysis/readend_count.pl \
			index/osa/{2}.fa \
			NJU_seq_output/{1}/rrna.filter.tmp {2} \
			>NJU_seq_output/{1}/rrna.filter.{2}.tsv
	' ::: NJU62{{60..63},{72..75}} Osa_{NC,{1..3}} \\
	::: 25s 18s 5-8s 5s
#	parallel -j 20 --keep-order '
#		perl NJU_seq/rrna_analysis/readend_count.pl \
#			index/sly/{2}.fa \
#			NJU_seq_output/{1}/rrna.filter.tmp {2} \
#			>NJU_seq_output/{1}/rrna.filter.{2}.tsv
#	' ::: NJU62{{52..55},{64..67}} \\
#	::: 25s 18s 5-8s 5s
	parallel -j 20 --keep-order '
		perl NJU_seq/rrna_analysis/readend_count.pl \
			index/mct/{2}.fa \
			NJU_seq_output/{1}/rrna.filter.tmp {2} \
			>NJU_seq_output/{1}/rrna.filter.{2}.tsv
	' ::: Sly_{{1..3},NC} NJU62{{52..55},{64..67}} \\
	::: 25s 18s 5-8s 5s
"

for RNA in 25s 18s 5-8s 5s; do
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJU62{32..35}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Ath_leaf_RF.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJU62{48..51}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Ath_leaf.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJU62{52..55}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Sly_leaf_RF.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJU62{56..59}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Gma_leaf_RF.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJU62{60..63}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Osa_leaf_RF.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJU62{64..67}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Sly_leaf.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJU62{68..71}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Gma_leaf.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJU62{72..75}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Osa_leaf.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJUzp{06,{03..05}}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Ath_zp.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJUzp{10,{07..09}}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Ath_zpdc.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJUzp{14,{11..13}}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Gma_zp.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NJUzp{18,{15..17}}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Gma_zpsmv.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/Sly_{NC,{1..3}}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Sly_ly.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/Col_SA_{NC,{1..3}}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Ath_SA.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/Osa_{NC,{1..3}}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Osa_ly.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NBG00{01..04}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Gma_wy_G1.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NBG00{05..08}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Gma_wy_G2.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NBG00{09..12}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Gma_wy_G3.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NBG00{13..16}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Gma_wy_G4.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NBG00{17..20}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Gma_wy_G5.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NBG00{21..24}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Gma_wy_G6.tsv &
	perl rrna_analysis_score.pl \
		NJU_seq_output/NBG00{25..28}/rrna.filter.${RNA}.tsv \
		>NJU_seq_output/rrna.filter.${RNA}.Gma_wy_G7.tsv &
	wait
done

for SPECIES in Ath_{leaf_RF,leaf,zp,zpdc,SA} Gma_{leaf_RF,leaf,zp,zpsmv,wy_G{1..7}} Osa_{leaf_RF,leaf,ly} Sly_{leaf_RF,leaf,ly}; do
	awk '
		FNR==1 && $1 ~ /^#Position/ {
			if (!printed_header++) {
				sub(/^#Position/, "Position", $1)
				print "#rRNA\t"$0
			}
			next
		}
		{
			rrna = (FILENAME ~ /18s/) ? "18s" :
						 (FILENAME ~ /25s/) ? "25s" :
						 (FILENAME ~ /5-8s/) ? "5-8s" :
						 (FILENAME ~ /5s/) ? "5s" : "NA"
			print rrna "\t" $0
		}
	' NJU_seq_output/rrna.filter.{25,18,5-8,5}s."${SPECIES}".tsv \
		>NJU_seq_output/rrna.filter."${SPECIES}".tsv
	awk '$12=="meet"&&$14=="meet"{print $1"\t" $2"\t"$3}' NJU_seq_output/rrna.filter."${SPECIES}".tsv >NJU_seq_output/rrna.Nm."${SPECIES}".tsv
	wc -l NJU_seq_output/rrna.Nm."${SPECIES}".tsv
done

#35 NJU_seq_output/rrna.Nm.Ath_leaf_RF.tsv
#58 NJU_seq_output/rrna.Nm.Ath_leaf.tsv
#93 NJU_seq_output/rrna.Nm.Ath_zp.tsv
#71 NJU_seq_output/rrna.Nm.Ath_zpdc.tsv
#67 NJU_seq_output/rrna.Nm.Ath_SA.tsv
#53 NJU_seq_output/rrna.Nm.Gma_leaf_RF.tsv
#55 NJU_seq_output/rrna.Nm.Gma_leaf.tsv
#60 NJU_seq_output/rrna.Nm.Gma_zp.tsv
#61 NJU_seq_output/rrna.Nm.Gma_zpsmv.tsv
#63 NJU_seq_output/rrna.Nm.Gma_wy_G1.tsv
#41 NJU_seq_output/rrna.Nm.Gma_wy_G2.tsv
#46 NJU_seq_output/rrna.Nm.Gma_wy_G3.tsv
#19 NJU_seq_output/rrna.Nm.Gma_wy_G4.tsv
#42 NJU_seq_output/rrna.Nm.Gma_wy_G5.tsv
#44 NJU_seq_output/rrna.Nm.Gma_wy_G6.tsv
#46 NJU_seq_output/rrna.Nm.Gma_wy_G7.tsv
#48 NJU_seq_output/rrna.Nm.Osa_leaf_RF.tsv
#59 NJU_seq_output/rrna.Nm.Osa_leaf.tsv
#71 NJU_seq_output/rrna.Nm.Osa_ly.tsv
#52 NJU_seq_output/rrna.Nm.Sly_leaf_RF.tsv
#60 NJU_seq_output/rrna.Nm.Sly_leaf.tsv
#53 NJU_seq_output/rrna.Nm.Sly_ly.tsv

for SAMPLE in NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}; do
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

for SAMPLE in NJU62{{56..59},{68..71}} NJUzp{11..18} NBG00{01..28}; do
	bsub -n 24 -J "$SAMPLE"_mrna "
		bowtie2 --end-to-end -p 20 -a -t \\
			--no-unal --no-mixed --no-discordant \\
			 -D 20 -R 3 -i S,1,0.75 --maxins 120 \\
			-N 0 -L 10 --score-min L,0.8,-0.8 \\
			--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \\
			--xeq -x index/gma/gma_protein_coding \\
			-1 NJU_seq_reads/${SAMPLE}/R1.filter.mrna.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.filter.mrna.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/mrna.filter.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/mrna.filter.bowtie2.log
	"
done

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

#for SAMPLE in NJU62{{52..55},{64..67}}; do
#	bsub -n 24 -J "$SAMPLE"_mrna "
#		bowtie2 --end-to-end -p 20 -a -t \\
#			--no-unal --no-mixed --no-discordant \\
#			 -D 20 -R 3 -i S,1,0.75 --maxins 120 \\
#			-N 0 -L 10 --score-min L,0.8,-0.8 \\
#			--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \\
#			--xeq -x index/sly/sly_protein_coding \\
#			-1 NJU_seq_reads/${SAMPLE}/R1.filter.mrna.fq.gz \\
#			-2 NJU_seq_reads/${SAMPLE}/R2.filter.mrna.fq.gz \\
#			-S NJU_seq_output/${SAMPLE}/mrna.filter.sam \\
#			2>&1 | tee NJU_seq_output/${SAMPLE}/mrna.filter.bowtie2.log
#	"
#done

for SAMPLE in NJU62{{52..55},{64..67}} Sly_{{1..3},NC}; do
	bsub -n 24 -J "$SAMPLE"_mrna "
		bowtie2 --end-to-end -p 20 -a -t \\
			--no-unal --no-mixed --no-discordant \\
			 -D 20 -R 3 -i S,1,0.75 --maxins 120 \\
			-N 0 -L 10 --score-min L,0.8,-0.8 \\
			--rdg 8,8 --rfg 8,8 --np 12 --mp 12,12 --ignore-quals \\
			--xeq -x index/mct/mct_protein_coding \\
			-1 NJU_seq_reads/${SAMPLE}/R1.filter.mrna.fq.gz \\
			-2 NJU_seq_reads/${SAMPLE}/R2.filter.mrna.fq.gz \\
			-S NJU_seq_output/${SAMPLE}/mrna.filter.sam \\
			2>&1 | tee NJU_seq_output/${SAMPLE}/mrna.filter.bowtie2.log
	"
done

bsub -n 128 -q amd_milan -J sam_filter "
	parallel -j 100 --keep-order '
		samtools view -h -F 128 NJU_seq_output/{1}/mrna.filter.sam \\
			| samtools view -F 16 - \\
			| awk '\''\$6!=\"*\"&&\$7==\"=\"&&\$4==\$8{print \$1 \"\\t\" \$3 \"\\t\" \$4 \"\\t\" \$6 \"\\t\" \$10}'\'' \\
			| perl matchquality_judge.pl \\
			| perl NJU_seq/mrna_analysis/multimatch_judge.pl \\
				>NJU_seq_output/{1}/mrna.filter.tmp
	' ::: NJU62{{32..35},{48..75}} NBG00{01..28} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}}
"

for SAMPLE in NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}; do
	bsub -n 1 -J "$SAMPLE"_mrna "
		perl NJU_seq/mrna_analysis/dedup.pl \\
			--refstr \"Parent=transcript:\" \\
			--info index/ath/exon.info \\
			<NJU_seq_output/${SAMPLE}/mrna.filter.tmp \\
			>NJU_seq_output/${SAMPLE}/mrna.dedup.filter.tmp
	"
done

#for SAMPLE in NJU62{{56..59},{68..71}} NJUzp{11..18} NBG00{01..28}; do
#	bsub -n 1 -J "$SAMPLE"_mrna "
#		perl NJU_seq/mrna_analysis/dedup.pl \\
#			--refstr \"Parent=\" \\
#			--transid \"Glyma.\" \\
#			--info index/gma/exon.info \\
#			<NJU_seq_output/${SAMPLE}/mrna.filter.tmp \\
#			>NJU_seq_output/${SAMPLE}/mrna.dedup.filter.tmp
#	"
#done

bsub -n 24 -J dedup "
	parallel -j 20 --keep-order '
		perl NJU_seq/mrna_analysis/dedup.pl \\
			--refstr \"Parent=\" \\
			--transid \"Glyma.\" \\
			--info index/gma/exon.info \\
			<NJU_seq_output/{1}/mrna.filter.tmp \\
			>NJU_seq_output/{1}/mrna.dedup.filter.tmp
	' ::: NJU62{{56..59},{68..71}} NJUzp{11..18} NBG00{01..28}
"

for SAMPLE in NJU62{{60..63},{72..75}} Osa_{NC,{1..3}}; do
	sed -i 's/\.fgenesh\.mRNA/_fgenesh_mRNA/g' NJU_seq_output/"${SAMPLE}"/mrna.filter.tmp
done

for SAMPLE in NJU62{{60..63},{72..75}} Osa_{NC,{1..3}}; do
	bsub -n 1 -J "$SAMPLE"_mrna "
		perl NJU_seq/mrna_analysis/dedup.pl \\
			--refstr \"Parent=\" \\
			--info index/osa/exon.info \\
			<NJU_seq_output/${SAMPLE}/mrna.filter.tmp \\
			>NJU_seq_output/${SAMPLE}/mrna.dedup.filter.tmp
	"
done

#for SAMPLE in NJU62{{52..55},{64..67}}; do
#	bsub -n 1 -J "$SAMPLE"_mrna "
#		perl NJU_seq/mrna_analysis/dedup.pl \\
#			--refstr \"Parent=\" \\
#			--transid \"Solyc\" \\
#			--info index/sly/exon.info \\
#			<NJU_seq_output/${SAMPLE}/mrna.filter.tmp \\
#			>NJU_seq_output/${SAMPLE}/mrna.dedup.filter.tmp
#	"
#done

for SAMPLE in NJU62{{52..55},{64..67}} Sly_{{1..3},NC}; do
	bsub -n 1 -J "$SAMPLE"_mrna "
		perl NJU_seq/mrna_analysis/dedup.pl \\
			--refstr \"Parent=\" \\
			--transid \"SolyMT\" \\
			--info index/mct/exon.info \\
			<NJU_seq_output/${SAMPLE}/mrna.filter.tmp \\
			>NJU_seq_output/${SAMPLE}/mrna.dedup.filter.tmp
	"
done

#for SAMPLE in NJU62{{32..35},{48..75}} NBG00{01..28} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}}; do
#	bsub -n 1 -J "$SAMPLE"_mrna "
#		perl NJU_seq/mrna_analysis/almostuniquematch.pl \\
#			NJU_seq_reads/${SAMPLE}/R1.filter.mrna.fq.gz \\
#			NJU_seq_output/${SAMPLE}/mrna.dedup.filter.tmp \\
#			NJU_seq_output/${SAMPLE}/mrna.almostuniquematch.filter.tmp
#		perl NJU_seq/mrna_analysis/count.pl \\
#			NJU_seq_output/${SAMPLE}/mrna.almostuniquematch.filter.tmp \\
#			| sort -k1,1 -k2,2n \\
#			>NJU_seq_output/${SAMPLE}/mrna.count.filter.tmp
#	"
#done

bsub -n 128 -q amd_milan -J dedup_count "
	parallel -j 100 --keep-order '
		perl NJU_seq/mrna_analysis/almostuniquematch.pl \\
			NJU_seq_reads/{1}/R1.filter.mrna.fq.gz \\
			NJU_seq_output/{1}/mrna.dedup.filter.tmp \\
			NJU_seq_output/{1}/mrna.almostuniquematch.filter.tmp
		perl NJU_seq/mrna_analysis/count.pl \\
			NJU_seq_output/{1}/mrna.almostuniquematch.filter.tmp \\
			| sort -k1,1 -k2,2n \\
			>NJU_seq_output/{1}/mrna.count.filter.tmp
	' ::: NJU62{{32..35},{48..75}} NBG00{01..28} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}}
"

#for SAMPLE in NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}; do
#	bsub -n 1 -J "$SAMPLE"_mrna "
#		perl NJU_seq/mrna_analysis/merge.pl \\
#			--refstr \"Parent=transcript:\" \\
#			--geneid \"AT\" \\
#			--transid \"AT\" \\
#			-i NJU_seq_output/${SAMPLE}/mrna.count.filter.tmp \\
#			--stdout \\
#			<index/ath/exon.info \\
#			| sort -k1,1 -k2,2n \\
#			  >NJU_seq_output/${SAMPLE}/mrna.filter.tsv
#	"
#done

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
	' ::: NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}
"

#for SAMPLE in NJU62{{56..59},{68..71}} NJUzp{11..18} NBG00{01..28}; do
#	bsub -n 1 -J "$SAMPLE"_mrna "
#		perl NJU_seq/mrna_analysis/merge.pl \\
#			--refstr \"Parent=\" \\
#			--geneid \"Glyma.\" \\
#			--transid \"Glyma.\" \\
#			-i NJU_seq_output/${SAMPLE}/mrna.count.filter.tmp \\
#			--stdout <index/gma/exon.info \\
#			| sort -k1,1 -k2,2n \\
#			  >NJU_seq_output/${SAMPLE}/mrna.filter.tsv
#	"
#done

bsub -n 24 -J sam_filter "
	parallel -j 20 --keep-order '
		perl NJU_seq/mrna_analysis/merge.pl \\
			--refstr \"Parent=\" \\
			--geneid \"Glyma.\" \\
			--transid \"Glyma.\" \\
			-i NJU_seq_output/{1}/mrna.count.filter.tmp \\
			--stdout <index/gma/exon.info \\
			| sort -k1,1 -k2,2n \\
			  >NJU_seq_output/{1}/mrna.filter.tsv
	' ::: NJU62{{56..59},{68..71}} NJUzp{11..18} NBG00{01..28}
"

#for SAMPLE in NJU62{{60..63},{72..75}} Osa_{NC,{1..3}}; do
#	bsub -n 1 -J "$SAMPLE"_mrna "
#		perl NJU_seq/mrna_analysis/merge.pl \\
#			--refstr \"Parent=\" \\
#			-i NJU_seq_output/${SAMPLE}/mrna.count.filter.tmp \\
#			--stdout <index/osa/exon.info \\
#			| sort -k1,1 -k2,2n \\
#			  >NJU_seq_output/${SAMPLE}/mrna.filter.tsv
#	"
#done

bsub -n 24 -J sam_filter "
	parallel -j 20 --keep-order '
		perl NJU_seq/mrna_analysis/merge.pl \\
			--refstr \"Parent=\" \\
			-i NJU_seq_output/{1}/mrna.count.filter.tmp \\
			--stdout <index/osa/exon.info \\
			| sort -k1,1 -k2,2n \\
			  >NJU_seq_output/{1}/mrna.filter.tsv
	' ::: NJU62{{60..63},{72..75}} Osa_{NC,{1..3}}
"

#for SAMPLE in NJU62{{52..55},{64..67}}; do
#	bsub -n 1 -J "$SAMPLE"_mrna "
#		perl NJU_seq/mrna_analysis/merge.pl \\
#			--refstr \"Parent=\" \\
#			--geneid \"Solyc\" \\
#			--transid \"Solyc\" \\
#			-i NJU_seq_output/${SAMPLE}/mrna.count.filter.tmp \\
#			--stdout <index/sly/exon.info \\
#			| sort -k1,1 -k2,2n \\
#			  >NJU_seq_output/${SAMPLE}/mrna.filter.tsv
#	"
#done

#for SAMPLE in NJU62{{52..55},{64..67}} Sly_{{1..3},NC}; do
#	bsub -n 1 -J "$SAMPLE"_mrna "
#		perl NJU_seq/mrna_analysis/merge.pl \\
#			--refstr \"Parent=\" \\
#			--geneid \"SolyMT\" \\
#			--transid \"SolyMT\" \\
#			-i NJU_seq_output/${SAMPLE}/mrna.count.filter.tmp \\
#			--stdout <index/mct/exon.info \\
#			| sort -k1,1 -k2,2n \\
#			  >NJU_seq_output/${SAMPLE}/mrna.filter.tsv
#	"
#done

bsub -n 24 -J sam_filter "
	parallel -j 20 --keep-order '
		perl NJU_seq/mrna_analysis/merge.pl \\
			--refstr \"Parent=\" \\
			--geneid \"SolyMT\" \\
			--transid \"SolyMT\" \\
			-i NJU_seq_output/{1}/mrna.count.filter.tmp \\
			--stdout <index/mct/exon.info \\
			| sort -k1,1 -k2,2n \\
			  >NJU_seq_output/{1}/mrna.filter.tsv
	' ::: NJU62{{52..55},{64..67}} Sly_{{1..3},NC}
"

perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJU62{32..35}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Ath_leaf_RF.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJU62{48..51}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Ath_leaf.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJU62{52..55}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Sly_leaf_RF.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJU62{56..59}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Gma_leaf_RF.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJU62{60..63}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Osa_leaf_RF.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJU62{64..67}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Sly_leaf.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJU62{68..71}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Gma_leaf.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJU62{72..75}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Osa_leaf.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJUzp{06,{03..05}}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Ath_zp.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJUzp{10,{07..09}}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Ath_zpdc.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJUzp{14,{11..13}}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Gma_zp.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NJUzp{18,{15..17}}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Gma_zpsmv.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/Sly_{NC,{1..3}}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Sly_ly.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/Col_SA_{NC,{1..3}}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Ath_SA.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/Osa_{NC,{1..3}}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Osa_ly.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NBG00{01..04}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Gma_wy_G1.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NBG00{05..08}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Gma_wy_G2.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NBG00{09..12}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Gma_wy_G3.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NBG00{13..16}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Gma_wy_G4.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NBG00{17..20}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Gma_wy_G5.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NBG00{21..24}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Gma_wy_G6.tsv
perl NJU_seq/mrna_analysis/score.pl \
	NJU_seq_output/NBG00{25..28}/mrna.filter.tsv \
	>NJU_seq_output/mrna.filter.Gma_wy_G7.tsv

for SPECIES in Ath_{leaf_RF,zp,zpdc,SA} Gma_{leaf_RF,zp,zpsmv,wy_G{1..7}} Osa_{leaf_RF,ly} Sly_{leaf_RF,ly}; do
	perl NJU_seq/presentation/signature_count.pl \
		NJU_seq_output/mrna.filter."${SPECIES}".tsv \
		NJU_seq_output/mrna.filter."${SPECIES}".signature.pdf
done

#for SPECIES in Sly_leaf_RF; do
#	for length in 20 100; do
#		for base in A G C T; do
#			perl NJU_seq/mrna_analysis/motif_nm.pl \
#				index/sly/sly.fa \
#				<(sed '1d' NJU_seq_output/mrna.filter.${SPECIES}.tsv | awk -v base=${base} '$5==base') \
#				${length} ${length} | awk -v OFS='' '
#				{header=$1; header=header":" $2; header=header"|" $4$5$6; header=header" " $7; seq=$NF; print ">" header; print seq}
#				' >NJU_seq_output/motif/mrna.filter.${SPECIES}.${length}nt.${base}.fa
#		done
#	done
#done
for SPECIES in Sly_leaf_RF Sly_ly; do
	for base in A G C T; do
		perl NJU_seq/mrna_analysis/motif_nm.pl \
			index/mct/mct.fa \
			<(sed '1d' NJU_seq_output/mrna.filter.${SPECIES}.tsv | awk -v base=${base} '$5==base') \
			20 20 | awk -v OFS='' '
				{header=$1; header=header":" $2; header=header"|" $4$5$6; header=header" " $7; seq=$NF; print ">" header; print seq}
				' >NJU_seq_output/motif/mrna.filter.${SPECIES}.20nt.${base}.fa
	done
done
for SPECIES in Ath_leaf_RF Ath_zp Ath_zpdc Ath_SA; do
	for base in A G C T; do
		perl NJU_seq/mrna_analysis/motif_nm.pl \
			index/ath/ath.fa \
			<(sed '1d' NJU_seq_output/mrna.filter.${SPECIES}.tsv | awk -v base=${base} '$5==base') \
			20 20 | awk -v OFS='' '
				{header=$1; header=header":" $2; header=header"|" $4$5$6; header=header" " $7; seq=$NF; print ">" header; print seq}
				' >NJU_seq_output/motif/mrna.filter.${SPECIES}.20nt.${base}.fa
	done
done
for SPECIES in Gma_leaf_RF Gma_zp Gma_zpsmv; do
	for base in A G C T; do
		perl NJU_seq/mrna_analysis/motif_nm.pl \
			index/gma/gma.fa \
			<(sed '1d' NJU_seq_output/mrna.filter."${SPECIES}".tsv | awk -v base=${base} '$5==base') \
			20 20 | awk -v OFS='' '
				{header=$1; header=header":" $2; header=header"|" $4$5$6; header=header" " $7; seq=$NF; print ">" header; print seq}
				' >NJU_seq_output/motif/mrna.filter."${SPECIES}".20nt.${base}.fa
	done
done
for SPECIES in Osa_leaf_RF Osa_ly; do
	for base in A G C T; do
		perl NJU_seq/mrna_analysis/motif_nm.pl \
			index/osa/osa.fa \
			<(sed '1d' NJU_seq_output/mrna.filter.${SPECIES}.tsv | awk -v base=${base} '$5==base') \
			20 20 | awk -v OFS='' '
				{header=$1; header=header":" $2; header=header"|" $4$5$6; header=header" " $7; seq=$NF; print ">" header; print seq}
				' >NJU_seq_output/motif/mrna.filter.${SPECIES}.20nt.${base}.fa
	done
done

for SPECIES in Sly_leaf_RF Sly_ly Ath_leaf_RF Ath_zp Ath_zpdc; do
	for base in A G C T; do
		./meme-5.5.9/src/meme -rna -mod zoops -objfun classic \
			-time 14400 -p 12 -minw 6 -maxw 20 \
			-markov_order 0 -evt 0.05 -nmotifs 30 \
			mrna.filter.${SPECIES}.20nt.${base}.fa \
			-oc ${SPECIES}_20nt_${base}
	done
done
