mkdir output
while IFS= read -r ecotype; do
	mkdir -p output/"${ecotype}"
	pigz -dc upload_genome/"${ecotype}".fasta.gz \
		>output/"${ecotype}"/chr.fasta

	gffread -x output/"${ecotype}"/cds.fasta \
		-g output/"${ecotype}"/chr.fasta \
		upload_gene_gff/"${ecotype}".*.protein_coding_genes.gff

	transeq -frame 1 \
		-sequence output/"${ecotype}"/cds.fasta \
		-outseq output/"${ecotype}"/proteins.faa

	sed -i 's/_1$//' output/"${ecotype}"/proteins.faa
done <ecotype.list

parallel -j 6 "
	/home/ivan/miniconda3/bin/perl nbs_hmm_blast.pl \
		--hmmquery NB-ARC.hmm \
		--blastp_database angiosperm_190aa_land_plant_ref_R_genes_classification \
		--pfam_database /home/ivan/fat/index/pfamA/Pfam-A.hmm \
		--output_dir output/{1}/ \
		--input_file output/{1}/proteins.faa

	java -jar ../NLR-Annotator/NLR-Annotator-v2.1b.jar \
		-t 6 -n 1000 -i output/{1}/chr.fasta \
		-x ../NLR-Annotator/src/mot.txt \
		-y ../NLR-Annotator/src/store.txt \
		-o output/{1}/NLR.anno.txt -g output/{1}/NLR.anno.gff \
		-b output/{1}/NLR.anno.bed -a output/{}/NLR.anno.fa
" :::: ecotype.list

echo -e "Ecotype\tCNL\tRNL\tTNL\tSum_BH\tSum_Anno" >ecotype_NLR_counts.tsv

while IFS= read -r ecotype; do
	awk '$3=="gene"' upload_gene_gff/"${ecotype}".*.protein_coding_genes.gff \
		| awk '{ match($9,/ID=([^;]+)/,a); if(a[1]!="") print $1"\t"$4"\t"$5"\t"a[1]; }' \
		| sort -k1,1 -k2,2n \
			>output/"${ecotype}"/genes.bed

	grep -F -f <(cut -f 1 output/"${ecotype}"/1st_blast_classification.txt | awk -F'.' '{print $1}') \
		output/"${ecotype}"/genes.bed \
		| sort -k1,1 -k2,2n \
			>output/"${ecotype}"/NLR.genes.bed

	perl calc_distance.pl output/"${ecotype}"/NLR.genes.bed >output/"${ecotype}"/NLR.gene_pairs_distance.tsv
	Rscript pairs_classify.R \
		output/"${ecotype}"/NLR.gene_pairs_distance.tsv \
		output/"${ecotype}"/NLR.gene_pairs_classified.tsv \
		output/"${ecotype}"/NLR.genes.bed \
		"${ecotype}" \
		output/"${ecotype}"/target_genes_classified.tsv

	perl make_clusters.pl \
		--bed output/"${ecotype}"/genes.bed \
		--target <(cut -f 1 output/"${ecotype}"/1st_blast_classification.txt | awk -F'.' '{print $1}') \
		--output output/"${ecotype}"/NLR.clusters.tsv \
		--max_gap 200000 \
		--max_non 1 \
		--h2h_thr 15000

	awk -vb="${ecotype}" '
		NR>1{
		split($6, a, ",");
		for(i in a){
			print $1"\t"$2"\t"$3"\t"$4"\t"a[i]"\t"$5"\t"$7"\t"b
		}}' output/"${ecotype}"/NLR.clusters.tsv \
		>output/"${ecotype}"/NLR.per_gene_classified.tsv

	sed -i '1iClusterID\tChr\tStart\tEnd\tGene\tSize\tCategory\tResource' \
		output/"${ecotype}"/NLR.per_gene_classified.tsv

	counts=$(grep -v "No_NLR" "output/${ecotype}/1st_blast_classification.txt" \
		| awk -F'_' '{print $1}' \
		| sort | uniq | cut -f 2 \
		| sort | uniq -c)
	cnl=0
	rnl=0
	tnl=0
	while read -r count type; do
		case "$type" in
		CNL) cnl=$count ;;
		RNL) rnl=$count ;;
		TNL) tnl=$count ;;
		esac
	done < <(echo "$counts")
	sum1=$(grep -v "No_NLR" "output/${ecotype}/1st_blast_classification.txt" \
		| awk -F'_' '{print $1}' \
		| sort | uniq | wc -l)
	sum2=$(wc -l <"output/${ecotype}/NLR.anno.txt")
	echo -e "${ecotype}\t${cnl}\t${rnl}\t${tnl}\t${sum1}\t${sum2}" >>ecotype_NLR_counts.tsv
done <ecotype.list
rm output/*/proteins.faa.p* \
	output/*/hmmtemp.txt \
	output/*/blastresult.txt \
	output/*/hmmresult.txt \
	output/*/*_out.fas \
	output/*/cds.fasta

tsv-join -H --data-fields 1 --key-fields 1 -a 2-6 \
	-f ecotype_NLR_counts.tsv seedexist.tsv \
	>ecotype_NLR_counts_land.tsv

tsv-join -H --data-fields 1 --key-fields 1 -a 2-11 \
	-f ecotype_NLR_counts_land.tsv ecotype.ordered.list \
	>ecotype_NLR_counts_land.ordered.tsv

tsv-append -H ./output/*/NLR.clusters.tsv \
	>All_clusters.tsv
tsv-append -H ./output/*/NLR.per_gene_classified.tsv \
	>All_target_genes_classified.tsv

gffread -x Col-PEK/cds.fasta \
	-g Col-PEK/file1.Col-PEK1.5_Chr1-5_20220523.fasta \
	Col-PEK/file2.Col-PEK1.5_Araport11-codinggene-27445_liftoff-cas_27561_revise.gff
transeq -frame 1 \
	-sequence Col-PEK/cds.fasta \
	-outseq proteins/Col-PEK.tmp.fa
sed -i 's/_1$//' proteins/Col-PEK.tmp.fa
seqkit fx2tab proteins/Col-PEK.tmp.fa \
	| awk '{split($1,a,"."); gene=a[1]; len=length($2); if(len>max[gene]){seq[gene]=$0; max[gene]=len}} END{for(g in seq) print seq[g]}' \
	| seqkit tab2fx >proteins/Col-PEK.fa
rm proteins/Col-PEK.tmp.fa

while IFS= read -r ecotype; do
	cp output/"${ecotype}"/proteins.faa proteins/"${ecotype}".fa
done <ecotype.list


orthofinder -f proteins/ -a 2 -S diamond_ultra_sens
