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
