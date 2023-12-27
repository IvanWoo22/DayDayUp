mkdir -p Ensembl/Atha
aria2c -x 12 https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz
pigz -dc Arabidopsis_thaliana.TAIR10.dna_sm.toplevel.fa.gz >Ensembl/Atha/toplevel.fa
egaz prepseq Ensembl/Atha

faops count Ensembl/Atha/toplevel.fa \
	| perl -nla -e 'next if $F[0] eq '\''total'\'';
    print $F[0] if $F[1] > 50000;
    print $F[0] if $F[1] > 5000 and $F[6]/$F[1] < 0.05;
    ' | uniq >Ensembl/Atha/listFile

faops some Ensembl/Atha/toplevel.fa Ensembl/Atha/listFile stdout \
	| faops filter -N stdin stdout \
	| faops split-name stdin Ensembl/Atha/.
rm Ensembl/Atha/toplevel.fa Ensembl/Atha/listFile
mv Ensembl/Atha/Mt.fa Ensembl/Atha/Mt.fa.skip
mv Ensembl/Atha/Pt.fa Ensembl/Atha/Pt.fa.skip
faops masked Ensembl/Atha/*.fa | spanr cover stdin | spanr stat --all Ensembl/Atha/chr.sizes stdin
#chrLength,size,coverage
#119146348,38123558,0.3200
#V52
#chrLength,size,coverage
#119146348,28109149,0.2359
#RepeatMask
#chrLength,size,coverage
#119146348,21267278,0.1785
egaz template Ensembl/Atha/ --self -o arabidopsis --taxon ../data/ensembl_taxon.csv --circos --parallel 16 -v

cd arabidopsis || exit
mkdir -p Pairwise
egaz lastz \
	--isself --set set01 -C 0 \
	--parallel 16 --verbose \
	/home/ivan/fat/MASED/revis/Ensembl/Atha /home/ivan/fat/MASED/revis/Ensembl/Atha \
	-o Pairwise/AthavsSelf

egaz lpcnam \
	--parallel 16 --verbose \
	/home/ivan/fat/MASED/revis/Ensembl/Atha /home/ivan/fat/MASED/revis/Ensembl/Atha \
	Pairwise/AthavsSelf

mkdir -p Results/Atha
mkdir -p Processing/Atha
ln -s /home/ivan/fat/MASED/revis/Ensembl/Atha/chr.fasta Processing/Atha/genome.fa
cp -f /home/ivan/fat/MASED/revis/Ensembl/Atha/chr.sizes Processing/Atha/chr.sizes

cd Processing/Atha || exit
fasops axt2fas \
	../../Pairwise/AthavsSelf/axtNet/*.axt.gz \
	-l 1000 -s chr.sizes -o stdout >axt.fas
fasops separate axt.fas -o . --nodash -s .sep.fasta
egaz exactmatch target.sep.fasta genome.fa --length 500 --discard 50 -o replace.target.tsv
fasops replace axt.fas replace.target.tsv -o axt.target.fas
egaz exactmatch query.sep.fasta genome.fa --length 500 --discard 50 -o replace.query.tsv
fasops replace axt.target.fas replace.query.tsv -o axt.correct.fas
fasops covers axt.correct.fas -o axt.correct.yml
spanr split axt.correct.yml -s .temp.yml -o .
spanr compare --op union target.temp.yml query.temp.yml -o axt.union.yml
spanr stat chr.sizes axt.union.yml -o union.csv
fasops links axt.correct.fas -o stdout | perl -nl -e "s/(target|query)\.//g; print;" >links.lastz.tsv
fasops separate axt.correct.fas --nodash --rc -o stdout \
	| perl -nl -e "/^>/ and s/^>(target|query)\./\>/; print;" \
	| faops filter -u stdin stdout \
	| faops filter -n 250 stdin stdout \
		>axt.gl.fasta

egaz blastn axt.gl.fasta genome.fa -o axt.bg.blast --parallel 8
egaz blastmatch axt.bg.blast -c 0.95 -o axt.bg.region --parallel 8
samtools faidx genome.fa -r axt.bg.region --continue \
	| perl -p -e "/^>/ and s/:/(+):/" >axt.bg.fasta
cat axt.gl.fasta axt.bg.fasta | faops filter -u stdin stdout \
	| faops filter -n 250 stdin stdout >axt.all.fasta
egaz blastn axt.all.fasta axt.all.fasta -o axt.all.blast --parallel 8
egaz blastlink axt.all.blast -c 0.95 -o links.blast.tsv --parallel 8

linkr sort links.lastz.tsv links.blast.tsv -o links.sort.tsv
linkr clean links.sort.tsv -o links.sort.clean.tsv
rgr merge links.sort.clean.tsv -c 0.95 -o links.merge.tsv

linkr clean links.sort.clean.tsv -r links.merge.tsv --bundle 500 -o links.clean.tsv
linkr connect links.clean.tsv -r 0.9 -o links.connect.tsv
linkr filter links.connect.tsv -r 0.8 -o links.filter.tsv

fasops create links.filter.tsv -g genome.fa -o multi.temp.fas
fasops refine multi.temp.fas --msa mafft -p 16 --chop 10 -o multi.refine.fas
fasops links multi.refine.fas -o stdout | linkr sort stdin -o stdout | linkr filter stdin -n 2-50 -o links.refine.tsv
fasops links multi.refine.fas -o stdout --best | linkr sort stdin -o links.best.tsv
fasops create links.best.tsv -g genome.fa --name Atha -o pair.temp.fas
fasops refine pair.temp.fas --msa mafft -p 16 -o pair.refine.fas

perl -nla -F"\t" -e "print for @F" <links.refine.tsv | spanr cover stdin -o cover.yml
echo "key,count" >links.count.csv
for n in 2 3 4-50; do
	linkr filter links.refine.tsv -n ${n} -o stdout \
		>links.copy${n}.tsv
	perl -nla -F"\t" -e "print for @F" <links.copy${n}.tsv | spanr cover stdin -o copy${n}.temp.yml
	wc -l links.copy${n}.tsv \
		| perl -nl -e "
            @fields = grep {/\S+/} split /\s+/;
            next unless @fields == 2;
            next unless \$fields[1] =~ /links\.([\w-]+)\.tsv/;
            printf qq{%s,%s\n}, \$1, \$fields[0];
        " \
			>>links.count.csv
	rm links.copy${n}.tsv
done
spanr merge copy2.temp.yml copy3.temp.yml copy4-50.temp.yml -o copy.yml
spanr stat chr.sizes copy.yml --all -o links.copy.csv
fasops mergecsv links.copy.csv links.count.csv --concat -o copy.csv
spanr stat chr.sizes cover.yml -o cover.yml.csv

cp cover.yml ../../Results/Atha/Atha.cover.yml
cp copy.yml ../../Results/Atha/Atha.copy.yml
mv cover.yml.csv ../../Results/Atha/Atha.cover.csv
mv copy.csv ../../Results/Atha/Atha.copy.csv
cp links.refine.tsv ../../Results/Atha/Atha.links.tsv
mv multi.refine.fas ../../Results/Atha/Atha.multi.fas
mv pair.refine.fas ../../Results/Atha/Atha.pair.fas
