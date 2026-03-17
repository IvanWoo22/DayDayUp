i=1
for SAMPLE in NJU62{33..35}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJU6232/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Ath_leaf_RF_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NJU62{49..51}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJU6248/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Ath_leaf_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NJU62{53..55}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJU6252/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Sly_leaf_RF_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NJU62{57..59}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJU6256/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Gma_leaf_RF_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NJU62{61..63}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJU6260/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Osa_leaf_RF_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NJU62{65..67}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJU6264/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Sly_leaf_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NJU62{69..71}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJU6268/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Gma_leaf_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NJU62{73..75}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJU6272/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Osa_leaf_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NJUzp{03..05}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJUzp06/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Ath_zp_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NJUzp{07..09}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJUzp10/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Ath_zpdc_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NJUzp{11..13}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJUzp14/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Gma_zp_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NJUzp{15..17}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NJUzp18/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Gma_zpsmv_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NBG00{02..04}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NBG0001/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Gma_wy_G1_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NBG00{06..08}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NBG0005/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Gma_wy_G2_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NBG00{10..12}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NBG0009/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Gma_wy_G3_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NBG00{14..16}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NBG0013/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Gma_wy_G4_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NBG00{18..20}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NBG0017/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Gma_wy_G5_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NBG00{22..24}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NBG0021/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Gma_wy_G6_${i}.tsv
	((i++))
done
i=1
for SAMPLE in NBG00{26..28}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/NBG0025/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Gma_wy_G7_${i}.tsv
	((i++))
done
i=1
for SAMPLE in Osa_{1..3}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/Osa_NC/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Osa_ly_${i}.tsv
	((i++))
done
i=1
for SAMPLE in Col_SA_{1..3}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/Col_SA_NC/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Ath_SA_${i}.tsv
	((i++))
done
i=1
for SAMPLE in Sly_{1..3}; do
	perl NJU_seq/mrna_analysis/score.pl \
		NJU_seq_output/Sly_NC/mrna.filter.tsv \
		NJU_seq_output/"${SAMPLE}"/mrna.filter.tsv \
		| cut -f 1-3 >NJU_seq_output/consistency/mrna.filter.Sly_ly_${i}.tsv
	((i++))
done

for SPECIES in Ath_{leaf_RF,leaf,zp,zpdc,SA} Gma_{leaf_RF,leaf,zp,zpsmv,wy_G{1..7}} Osa_{leaf_RF,leaf,ly} Sly_{leaf_RF,leaf,ly}; do
	Rscript Comparison.R consistency/mrna.filter."${SPECIES}".*.tsv --out consistency/mrna.filter."${SPECIES}".comparison.pdf
done

for SAMPLE in NJU62{{32..35},{48..75}} NBG00{01..28} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}}; do
	bsub -n 1 -J "$SAMPLE"_mrna "

	"
done

bsub -n 128 -q amd_milan -J sam_position "
	parallel -j 100 --keep-order '
		perl readend.pl \\
			NJU_seq_output/{1}/mrna.almostuniquematch.filter.tmp \\
			| sort -k1,1 -k2,2n \\
			>NJU_seq_output/{1}/mrna.endposition.filter.tmp
	' ::: NJU62{{32..35},{48..75}} NBG00{01..28} NJUzp{03..18} {Col_SA,Sly,Osa}_{NC,{1..3}}
"

for SAMPLE in NJU62{{32..35},{48..51}} NJUzp{03..10} Col_SA_{{1..3},NC}; do
	bsub -n 1 -J "$SAMPLE"_mrna "
		perl read_position.pl \\
			--refstr \"Parent=transcript:\" \\
			--geneid \"AT\" \\
			--transid \"AT\" \\
			-i NJU_seq_output/${SAMPLE}/mrna.endposition.filter.tmp \\
			--stdout \\
			<index/ath/exon.info \\
			| sort -k6,6 -k7,7n \\
			  >NJU_seq_output/${SAMPLE}/mrna.endposition.filter.tsv
	"
done

for SAMPLE in NJU62{{56..59},{68..71}} NJUzp{11..18} NBG00{01..28}; do
	bsub -n 1 -J "$SAMPLE"_mrna "
		perl read_position.pl \\
			--refstr \"Parent=\" \\
			--geneid \"Glyma.\" \\
			--transid \"Glyma.\" \\
			-i NJU_seq_output/${SAMPLE}/mrna.endposition.filter.tmp \\
			--stdout <index/gma/exon.info \\
			| sort -k6,6 -k7,7n \\
			  >NJU_seq_output/${SAMPLE}/mrna.endposition.filter.tsv
	"
done

for SAMPLE in NJU62{{60..63},{72..75}} Osa_{NC,{1..3}}; do
	bsub -n 1 -J "$SAMPLE"_mrna "
		perl read_position.pl \\
			--refstr \"Parent=\" \\
			-i NJU_seq_output/${SAMPLE}/mrna.endposition.filter.tmp \\
			--stdout <index/osa/exon.info \\
			| sort -k6,6 -k7,7n \\
			  >NJU_seq_output/${SAMPLE}/mrna.endposition.filter.tsv
	"
done

for SAMPLE in NJU62{{52..55},{64..67}} Sly_{{1..3},NC}; do
	bsub -n 1 -J "$SAMPLE"_mrna "
		perl read_position.pl \\
			--refstr \"Parent=\" \\
			--geneid \"SolyMT\" \\
			--transid \"SolyMT\" \\
			-i NJU_seq_output/${SAMPLE}/mrna.endposition.filter.tmp \\
			--stdout <index/mct/exon.info \\
			| sort -k6,6 -k7,7n \\
			  >NJU_seq_output/${SAMPLE}/mrna.endposition.filter.tsv
	"
done

for SPECIE in Ath Osa Gma Sly; do
	{
		grep "10-\[1\]-20" </home/ivan/fat/NJU_seq_fast/motif/${SPECIE}_*_20nt_A/mast/mast.txt \
			| grep -v "DIAGRAM" | awk '{print $1}' \
			| awk -F"[:|]" '{print $1"\t"$2}' | sort | uniq \
			| sort -k1,1 -k2,2n
		grep "11-\[1\]-20" </home/ivan/fat/NJU_seq_fast/motif/${SPECIE}_*_20nt_G/mast/mast.txt \
			| grep -v "DIAGRAM" | awk '{print $1}' \
			| awk -F"[:|]" '{print $1"\t"$2}' | sort | uniq \
			| sort -k1,1 -k2,2n
		grep "12-\[1\]-19" </home/ivan/fat/NJU_seq_fast/motif/${SPECIE}_*_20nt_C/mast/mast.txt \
			| grep -v "DIAGRAM" | awk '{print $1}' \
			| awk -F"[:|]" '{print $1"\t"$2}' | sort | uniq \
			| sort -k1,1 -k2,2n
		grep "10-\[1\]-19" </home/ivan/fat/NJU_seq_fast/motif/${SPECIE}_*_20nt_T/mastT1/mast.txt \
			| grep -v "DIAGRAM" | awk '{print $1}' \
			| awk -F"[:|]" '{print $1"\t"$2}' | sort | uniq \
			| sort -k1,1 -k2,2n
		grep "12-\[1\]-19" </home/ivan/fat/NJU_seq_fast/motif/${SPECIE}_*_20nt_T/mastT2/mast.txt \
			| grep -v "DIAGRAM" | awk '{print $1}' \
			| awk -F"[:|]" '{print $1"\t"$2}' | sort | uniq \
			| sort -k1,1 -k2,2n
		grep "12-\[1\]-19" </home/ivan/fat/NJU_seq_fast/motif/${SPECIE}_*_20nt_T/mastT3/mast.txt \
			| grep -v "DIAGRAM" | awk '{print $1}' \
			| awk -F"[:|]" '{print $1"\t"$2}' | sort | uniq \
			| sort -k1,1 -k2,2n
	} >${SPECIE}.motif.sites
done
