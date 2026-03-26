for i in U1 U2 U3 U4 U5 U11 U14 U15; do
	mkdir ${i}
	ln -sf /home/ivan/cater_data/Ara_Nm_terminal/${i}/*_1.fq.gz /home/ivan/fat/miR_Nm/${i}/R1.fq.gz
	ln -sf /home/ivan/cater_data/Ara_Nm_terminal/${i}/*_2.fq.gz /home/ivan/fat/miR_Nm/${i}/R2.fq.gz
done

for i in U1 U2 U3 U4 U5 U11 U14 U15; do
	cutadapt -a AGATCGGAAGAGCACA -A GATCGTCGGACTGTAGAA -O 5 -m 17 -e 0.18 -j 8 -o ${i}/clean_R1.fq.gz -p ${i}/clean_R2.fq.gz ${i}/R1.fq.gz ${i}/R2.fq.gz --discard-untrimmed >${i}/cutadapt.log 2>&1
done

parallel --keep-order --xapply -j 6 '
	pear -j 4 -f {}/clean_R1.fq.gz -r {}/clean_R2.fq.gz -o {} -n 10
	perl fastq2count.pl <{}.assembled.fastq >{}/reads_cutadapt.count
	perl mirna_count.pl {}/reads_cutadapt.count ath_mirna_neo.fa >{}/mir_cutadapt.count.tsv
	rm {}.assembled.fastq {}.discarded.fastq {}.unassembled.forward.fastq {}.unassembled.reverse.fastq
' ::: U1 U2 U3 U4 U5 U11 U14 U15

for i in U1 U2 U3 U4 U5 U11 U14 U15; do
	sed -i "1i\\Name\\t${i}\\t${i}_ployA" ${i}/mir_cutadapt.count.tsv
	sed -i 's/^>//g' ${i}/mir_cutadapt.count.tsv
done

cp U1/mir_cutadapt.count.tsv 0325.ath.count.tmp

for i in U2 U3 U4 U5 U11 U14 U15; do
	tsv-join -H --data-fields "Name" --key-fields "Name" -f ${i}/mir_cutadapt.count.tsv 0325.ath.count.tmp -a 2,3 >0325.merged.tsv
	mv 0325.merged.tsv 0325.ath.count.tmp
done

mv 0325.ath.count.tmp 0325.ath.merged.tsv
