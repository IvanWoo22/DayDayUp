perl split_lines.pl ../revis/RMasked/Results/Atha/Atha.links.tsv Atha_links2.txt
grep -v "Mt" Atha_links2.txt >Atha_links2.tsv

makeblastdb -in Atha.pep -dbtype prot -parse_seqids -out ATdb
for i in Asue Alyr Risl Esyr Cvio Ghir; do
	blastp -query ${i}.pep -db ATdb -out Atha_${i}.blast -evalue 10 -num_threads 18 -outfmt 6 -max_target_seqs 5
done

for i in Asue Alyr Risl Esyr Cvio Ghir; do
	cat Atha.gff ${i}.gff >Atha_${i}.gff
	../MCScanX/MCScanX -s 3 -m 2 -w 0 -a Atha_${i}
done
for i in Asue Alyr; do
	perl ../colline2yml_AL.pl Atha.gff Atha_${i}.collinearity >Atha_${i}.yml
done
for i in Risl Esyr Cvio Ghir; do
	perl ../colline2yml.pl Atha.gff Atha_${i}.collinearity >Atha_${i}.yml
done

perl ../data/classify_links_neo.pl Atha_links2.tsv \
	Atha_Asue.yml Atha_Alyr.yml Atha_Risl.yml Atha_Esyr.yml Atha_Cvio.yml Atha_Ghir.yml \
	0.5 >Atha_links2_timespot.tsv

cut -f 3,4,5,6,7,8 Atha_links2_timespot.tsv | sort | uniq -c \
	| awk -va=1 '{
    count0 = 0;
    count1 = 0;
    count2 = 0;
    count3 = 0;
    for (i = (a+1); i <= NF; i++) {
        if ($i == 0) {
            count0++;
        } else if ($i == 1) {
            count1++;
        } else if ($i == 2) {
            count2++;
        } else if ($i == 3) {
            count3++;
        }
    }
    if ((count0 >= 0 && count1 >= 0 && count2 == 0 && count3 == 0) || (count0 >= 0 && count1 == 0 && count2 >= 0 && count3 == 0)) {
        print;
    }
}' | awk '{print $2 $3 $4 $5 $6 $7 "\t" $1}' >time1.txt

for i in {2..7}; do
	cut -f 3,4,5,6,7,8 Atha_links2_timespot.tsv | sort | uniq -c \
		| awk -va="$i" '{
    count0 = 0;
    count1 = 0;
    count2 = 0;
    count3 = 0;
    for (i = (a+1); i <= NF; i++) {
        if ($i == 0) {
            count0++;
        } else if ($i == 1) {
            count1++;
        } else if ($i == 2) {
            count2++;
        } else if ($i == 3) {
            count3++;
        }
    }
    if ((count0 >= 0 && count1 >= 0 && count2 == 0 && count3 == 0) || (count0 >= 0 && count1 == 0 && count2 >= 0 && count3 == 0)) {
        print;
    }
}' | awk -va="$i" '$a==3{print $2 $3 $4 $5 $6 $7 "\t" $1}' >time"${i}".txt
done

for i in {1..7}; do
	perl fetch_timespot.pl \
		time"${i}".txt Atha_links2_timespot.tsv >Atha_links2_timespot_"${i}".tsv
done

for i in {1..7}; do
	perl ../promotor_intsec.pl \
		struc_promoter.bed Atha_links2_timespot_"${i}".tsv >Atha_links2_timespot_"${i}"_prom.tsv
	perl ../data/arg_meth_link.pl \
		../Memory/AT.beta.1.tsv Atha_links2_timespot_"${i}"_prom.tsv Atha_links2_timespot_"${i}"_prom_beta.tsv 3
done
