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

for chr in Chr{1..5}; do
        bgzip pggb."${chr}".fa
        samtools faidx pggb."${chr}".fa.gz
done


odgi extract -i pggb.Chr1.*.smooth.final.og \
	-r "Col-PEK#1#Chr1:23540998-23549864" \
	-o chr1.sub.og -t 20 --full-range -d 30000 -c 10000
odgi sort -i chr1.sub.og -o chr1.sub.sorted.og -p Ygsr -O -t 20
odgi layout -i chr1.sub.sorted.og -o chr1.sub.2d -t 20
odgi draw -i chr1.sub.sorted.og -c chr1.sub.2d -s chr1_sub_rainbow.svg -R 0.04 -B 1000 -w 20 -S 80
odgi viz -i chr1.sub.sorted.og -o chr1_sub_rainbow.png

