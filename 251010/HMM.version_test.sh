# ==> Test HMMER3/f [3.3 | Nov 2019]
parallel -j 6 "
	mkidr -p output/{1}-neo/
	/home/ivan/miniconda3/bin/perl nbs_hmm_blast.pl \
		--domain_id NB-ARC --hmmquery Pfam_index/PF00931.hmm \
		--blastp_database angiosperm_190aa_land_plant_ref_R_genes_classification \
		--pfam_database /home/ivan/fat/nlr_eco/Pfam_index/Pfam-A-neo.hmm \
		--output_dir output/{1}-neo/ \
		--input_file output/{1}/proteins.faa
" :::: ecotype.list

echo -e "Ecotype\tCNL\tRNL\tTNL\tSum_BH" >ecotype_NLR_counts.neo.tsv
while IFS= read -r ecotype; do
	counts=$(grep -v "No_NLR" "output/${ecotype}-neo/1st_blast_classification.txt" \
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
	sum1=$(grep -v "No_NLR" "output/${ecotype}-neo/1st_blast_classification.txt" \
		| awk -F'_' '{print $1}' \
		| sort | uniq | wc -l)
	echo -e "${ecotype}\t${cnl}\t${rnl}\t${tnl}\t${sum1}" >>ecotype_NLR_counts.neo.tsv
done <ecotype.list

