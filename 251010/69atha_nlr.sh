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

parallel -j 4 "
	mkdir -p output/{1}
	cp {1}/proteins.faa output/{1}/proteins.faa
" ::: Atha Alyr Crub Col-PEK

seqkit stats output/*/proteins.faa >proteins.stats.tsv

parallel -j 6 "
	/home/ivan/miniconda3/bin/perl nbs_hmm_blast.pl \
		--domain_id NB-ARC --hmmquery Pfam_index/PF00931.hmm \
		--blastp_database angiosperm_190aa_land_plant_ref_R_genes_classification \
		--pfam_database /home/ivan/fat/nlr_eco/Pfam_index/Pfam-A-neo.hmm \
		--output_dir output/{1}/ \
		--input_file output/{1}/proteins.faa

	java -jar ../NLR-Annotator/NLR-Annotator-v2.1b.jar \
		-t 6 -n 1000 -i output/{1}/chr.fasta \
		-x ../NLR-Annotator/src/mot.txt \
		-y ../NLR-Annotator/src/store.txt \
		-o output/{1}/NLR.anno.txt -g output/{1}/NLR.anno.gff \
		-b output/{1}/NLR.anno.bed -a output/{}/NLR.anno.fa
" :::: <(
	cat ecotype.list
	printf '%s\n' Atha Alyr Crub Col-PEK
)

while IFS= read -r ecotype; do
	awk '$3=="gene"' upload_gene_gff/"${ecotype}".*.protein_coding_genes.gff \
		| awk '{ match($9,/ID=([^;]+)/,a); if(a[1]!="") print $1"\t"$4"\t"$5"\t"a[1]; }' \
		| sort -k1,1 -k2,2n \
			>output/"${ecotype}"/genes.bed

	grep -F -f <(grep -v "No_NLR" output/"${ecotype}"/1st_blast_classification.txt | awk -F'.' '{print $1}') \
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
		--target <(grep -v "No_NLR" output/"${ecotype}"/1st_blast_classification.txt | awk -F'.' '{print $1}') \
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
done <ecotype.list

while IFS= read -r ecotype; do
	awk '$3=="gene"' "${ecotype}"/*gene_exons.gff3 \
		| awk '{ match($9,/ID=([^;]+)/,a); if(a[1]!="") print $1"\t"$4"\t"$5"\t"a[1]; }' \
		| sort -k1,1 -k2,2n \
			>output/"${ecotype}"/genes.bed

	grep -F -f <(grep -v "No_NLR" output/"${ecotype}"/1st_blast_classification.txt | awk -F'.' '{print $1}') \
		output/"${ecotype}"/genes.bed \
		| sort -k1,1 -k2,2n \
			>output/"${ecotype}"/NLR.genes.bed
done < <(
	printf '%s\n' Atha Alyr Col-PEK
)

sed -i 's/\.v2\.1$//' output/Alyr/NLR.genes.bed
sed -i 's/\.v2\.1$//' output/Alyr/genes.bed

while IFS= read -r ecotype; do
	awk '$3=="gene"' "${ecotype}"/*gene_exons.gff3 \
		| awk '{ match($9,/ID=([^;]+)/,a); if(a[1]!="") print $1"\t"$4"\t"$5"\t"a[1]; }' \
		| sed 's/\.v1\.1$//' | sort -k1,1 -k2,2n \
		>output/"${ecotype}"/genes.bed

	grep -F -f <(grep -v "No_NLR" output/"${ecotype}"/1st_blast_classification.txt | awk -F'.' '{print $2}') \
		output/"${ecotype}"/genes.bed \
		| sort -k1,1 -k2,2n \
			>output/"${ecotype}"/NLR.genes.bed
done < <(
	printf '%s\n' Crub
)

rm output/*/proteins.faa.p* \
	output/*/hmmtemp.txt \
	output/*/blastresult.txt \
	output/*/hmmresult.txt \
	output/*/*_out.fas \
	output/*/cds.fasta

echo -e "Ecotype\tCNL\tRNL\tTNL\tSum_BH\tSum_Anno" >ecotype_NLR_counts.tsv
while IFS= read -r ecotype; do
	counts=$(grep -v "No_NLR" output/"${ecotype}"/1st_blast_classification.txt \
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
	sum1=$(grep -v "No_NLR" output/"${ecotype}"/1st_blast_classification.txt \
		| awk -F'_' '{print $1}' \
		| sort | uniq | wc -l)
	sum2=$(wc -l <output/"${ecotype}"/NLR.anno.txt)
	echo -e "${ecotype}\t${cnl}\t${rnl}\t${tnl}\t${sum1}\t${sum2}" >>ecotype_NLR_counts.tsv
done < <(
	cat ecotype.list
	printf '%s\n' Atha Alyr Crub Col-PEK
)

tsv-join -H --data-fields 1 --key-fields 1 -a 2-11 \
	-f <(tsv-join -H --data-fields 1 --key-fields 1 -a 2-6 \
		-f seedexist.tsv ecotype_NLR_counts.tsv -w NA) \
	ecotype.ordered.list \
	>ecotype_NLR_counts_land.ordered.tsv

tsv-append -H ./output/*/NLR.clusters.tsv \
	>All_clusters.tsv
tsv-append -H ./output/*/NLR.per_gene_classified.tsv \
	>All_target_genes_classified.tsv

while IFS= read -r ecotype; do
	faops some \
		output/"${ecotype}"/proteins.faa \
		<(grep -v "No_NLR" output/"${ecotype}"/1st_blast_classification.txt | cut -f 1) \
		output/"${ecotype}"/NLR.proteins.faa
done < <(
	cat ecotype.list
	printf '%s\n' Atha Alyr Crub Col-PEK
)

seqkit stats output/*/NLR.proteins.faa >NLR.proteins.stats.tsv

while IFS= read -r ecotype; do
	hmmsearch \
		--tblout output/"${ecotype}"/nb_results.tbl \
		Pfam_index/PF00931.hmm \
		output/"${ecotype}"/NLR.proteins.faa
done < <(
	cat ecotype.list
	printf '%s\n' Atha Alyr Crub Col-PEK
)

# ==> orthofinder
while IFS= read -r ecotype; do
	cp output/"${ecotype}"/proteins.faa proteins/"${ecotype}".fa
done <ecotype.list
sed -i 's/^>\([^ ]*\).*/>\1/' proteins/Alyr.fa
sed -i 's/^>\([^ ]*\).*/>\1/' proteins/Crub.fa
sed -i 's/^>\([^|]*\).*/>\1/' proteins/Atha.fa

orthofinder -f proteins/ -a 2 -S diamond_ultra_sens

awk -F"\t" -v thresh=0.85 'NR==1{N=NF-1; next}{
  present=0; ok=1;
  for(i=2;i<=NF;i++){
    if($i!="") present++;
    if($i ~ /,/) ok=0;
  }
  if(ok && present>=thresh*N) print $1
}' proteins/OrthoFinder/Results_Nov14/Orthogroups/Orthogroups.tsv >singlecopy_ogs_85pct.txt

cut -f 4 output/*/NLR.genes.bed \
	| sort | uniq >all.nlr.genes.list
grep -f all.nlr.genes.list \
	<(cat /home/ivan/fat/nlr_eco/proteins/OrthoFinder/Results_Nov14/Orthogroups/Orthogroups.tsv \
		/home/ivan/fat/nlr_eco/proteins/OrthoFinder/Results_Nov14/Orthogroups/Orthogroups_UnassignedGenes.tsv) \
	| cut -f 1 | sort | uniq >all.nlr.ogs.list
grep -f all.nlr.ogs.list \
	<(cat /home/ivan/fat/nlr_eco/proteins/OrthoFinder/Results_Nov14/Orthogroups/Orthogroups.tsv \
		/home/ivan/fat/nlr_eco/proteins/OrthoFinder/Results_Nov14/Orthogroups/Orthogroups_UnassignedGenes.tsv) \
	| awk -F "\t" '{for(i=2;i<=NF;i++) print $i}' \
	| awk -F ", " '{for(i=1;i<=NF;i++) print $i}' \
	| sort | uniq >all.nlr.genes_byogs.list

while IFS= read -r ecotype; do
	mkdir -p output/"${ecotype}"/chrs
	faops split-name output/"${ecotype}"/chr.fasta output/"${ecotype}"/chrs/.
	sed -i "s/Chr/${ecotype}#1#Chr/g" output/"${ecotype}"/chrs/Chr*.fa
	for chr in Chr{1..5}; do
		cat output/"${ecotype}"/chrs/"${chr}".fa >pggb."${chr}".fa
	done
done < <(
	printf '%s\n' Col-PEK
	cat ecotype.selected.list
)

for chr in Chr{1..5}; do
	bgzip pggb."${chr}".fa
	samtools faidx pggb."${chr}".fa.gz
done

for chr in Chr{1..5}; do
	mkdir -p pggb_out/"${chr}"
	pggb \
		-i pggb."${chr}".fa.gz \
		-p 95 -s 2000 \
		-G 5000 -n 26 -m \
		-t 16 -k 311 -O 0.03 \
		-o pggb_out/"${chr}"/.
done
