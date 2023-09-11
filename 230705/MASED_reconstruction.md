# MASED

**MASED** - **M**ethylation **A**nalysis of **SE**gmental **D**uplications.

## System Information

#### Hardware

#### Software

```shell
cd ~
git clone https://github.com/IvanWoo22/perbool.git

mkdir ucsc_tools
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToTwoBit ./faToTwoBit

brew tap wang-q/tap
brew install faops sparsemem samtools

wget https://github.com/wang-q/intspan/archive/refs/tags/v0.7.0.tar.gz
tar -zvxf v0.7.0.tar.gz
cd intspan-0.7.0
cargo install --force --path .

singularity pull docker://wangq/egaz:master
echo "alias egaz='singularity run \$HOME/egaz_master.sif egaz'" >>.zshrc
source .zshrc
```

## Self alignment for segment duplication of Arabidopsis thaliana

#### Genome preparation

| Short Organism Name | Genome/Metagenome Name         | filename                                                    | JGI Grouping ID | directory/path                                                          |
|---------------------|--------------------------------|-------------------------------------------------------------|-----------------|-------------------------------------------------------------------------|
| Alyrata             | Arabidopsis lyrata v2.1        | Alyrata_384_v2.1.protein_primaryTranscriptOnly.fa.gz        | Phytozome-384   | Phytozome/PhytozomeV12/Alyrata/annotation                               |
| Alyrata             | Arabidopsis lyrata v2.1        | Alyrata_384_v2.1.gene.gff3.gz                               | Phytozome-384   | Phytozome/PhytozomeV12/Alyrata/annotation                               |
| Athaliana           | Arabidopsis thaliana Araport11 | Athaliana_447_TAIR10.fa.gz                                  | Phytozome-447   | Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/assembly   |
| Athaliana           | Arabidopsis thaliana Araport11 | Athaliana_447_Araport11.gene.gff3.gz                        | Phytozome-447   | Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/annotation |
| Athaliana           | Arabidopsis thaliana Araport11 | Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa.gz | Phytozome-447   | Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/annotation |
| BrapaFPsc           | Brassica rapa FPsc v1.3        | BrapaFPsc_277_v1.3.gene.gff3.gz                             | Phytozome-277   | Phytozome/PhytozomeV10/BrapaFPsc/annotation                             |
| BrapaFPsc           | Brassica rapa FPsc v1.3        | BrapaFPsc_277_v1.3.protein_primaryTranscriptOnly.fa.gz      | Phytozome-277   | Phytozome/PhytozomeV10/BrapaFPsc/annotation                             |
| Cpapaya             | Carica papaya ASGPBv0.4        | Cpapaya_113_ASGPBv0.4.protein_primaryTranscriptOnly.fa.gz   | Phytozome-113   | Phytozome/PhytozomeV11/Cpapaya/annotation                               |
| Cpapaya             | Carica papaya ASGPBv0.4        | Cpapaya_113_ASGPBv0.4.gene.gff3.gz                          | Phytozome-113   | Phytozome/PhytozomeV11/Cpapaya/annotation                               |
| Crubella            | Capsella rubella v1.1          | Crubella_474_v1.1.gene.gff3.gz                              | Phytozome-474   | Phytozome                                                               |
| Crubella            | Capsella rubella v1.1          | Crubella_474_v1.1.protein_primaryTranscriptOnly.fa.gz       | Phytozome-474   | Phytozome                                                               |
| Tcacao              | Theobroma cacao v2.1           | Tcacao_523_v2.1.gene.gff3.gz                                | Phytozome-523   | Phytozome/PhytozomeV13/Tcacao/v2.1/annotation                           |
| Tcacao              | Theobroma cacao v2.1           | Tcacao_523_v2.1.protein_primaryTranscriptOnly.fa.gz         | Phytozome-523   | Phytozome/PhytozomeV13/Tcacao/v2.1/annotation                           |

We obtained all genome files from JGI.

```shell
git clone https://github.com/IvanWoo22/MASED.git && cd ./MASED && rm -rf .git

# Take a look of chromosome names.
pigz -dc ./data/Atha/Atha.fa.gz | grep "^>"
#>Chr1 CHROMOSOME dumped from ADB: Jun/20/09 14:53; last updated: 2009-02-02
#>Chr2 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02
#>Chr3 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02
#>Chr4 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02
#>Chr5 CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2009-02-02
#>ChrC CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2005-06-03
#>ChrM CHROMOSOME dumped from ADB: Jun/20/09 14:54; last updated: 2005-06-03

pigz -dc ./data/Atha/Atha.fa.gz | perl ~/perbool/fetch_fasta.pl -s 'Chr[1-5]' --stdin >./data/Atha_raw_genome.fa
parallel --keep-order --xapply -j 1 '
  sed -i "s/^>{1}.*/>{2}/" ./data/Atha_raw_genome.fa
' ::: Chr{1..5} ::: chr{1..5}

# Check the names
grep "^>" ./data/Atha_raw_genome.fa
#>Chr1
#>Chr2
#>Chr3
#>Chr4
#>Chr5

mkdir Atha
faops filter -N -s ./data/Atha_raw_genome.fa stdout | faops split-name stdin Atha/.
egaz repeatmasker ./Atha/*.fa -o ./Atha/. --gff --parallel 12
faops size ./Atha/*.fa >./Atha/chr.sizes
faops filter -U stdin ./Atha/chr.fasta <./Atha/*.fa
~/ucsc_tools/faToTwoBit ./Atha/*.fa ./Atha/chr.2bit

grep -Pv "ChrC|ChrM" ./data/Atha/Atha.gff3 >./data/Atha_raw_genome.gff
parallel --keep-order --xapply -j 1 '
  sed -i "s/{1}/{2}/" ./data/Atha_raw_genome.gff
' ::: Chr{1..5} ::: chr{1..5}
runlist gff --tag CDS --remove ./data/Atha_raw_genome.gff -o cds.yml
runlist gff --remove ./Atha/*.rm.gff -o repeat.yml
runlist merge repeat.yml cds.yml -o ./Atha/anno.yml
rm repeat.yml cds.yml ./Atha/*.rm.gff ./Atha/*.rm.out
```

#### Self-align with lastz

```shell
#This create scripts including code below.
egaz template \
  --parallel 12 --verbose \
  ./Atha/. --self --circos \
  --taxon ./ensembl_taxon.csv \
  -o ATFSD

egaz lastz \
  --isself --set set01 -C 0 \
  --parallel 12 --verbose \
  ./Atha ./Atha \
  -o ATFSD/Pairwise/AthavsSelf

egaz lpcnam \
  --parallel 12 --verbose \
  ./Atha ./Atha \
  -o ATFSD/Pairwise/AthavsSelf
```

#### Self-align with blastn

```shell
mkdir ATFSD/Processing
ln -s ${PWD}/Atha/chr.fasta ATFSD/Processing/genome.fa
cp -f ${PWD}/Atha/chr.sizes ATFSD/Processing/chr.sizes
cd ATFSD/Processing || exit

# axt2fas
fasops axt2fas \
  ../Pairwise/AthavsSelf/axtNet/*.axt.gz \
  -l 1000 -s chr.sizes -o stdout >axt.fas
fasops separate axt.fas -o . --nodash -s .sep.fasta

# Target positions
egaz exactmatch target.sep.fasta genome.fa \
  --length 500 --discard 50 -o replace.target.tsv
fasops replace axt.fas replace.target.tsv -o axt.target.fas

# Query positions
egaz exactmatch query.sep.fasta genome.fa \
  --length 500 --discard 50 -o replace.query.tsv
fasops replace axt.target.fas replace.query.tsv -o axt.correct.fas

# Coverage stats
fasops covers axt.correct.fas -o axt.correct.yml
spanr split axt.correct.yml -s .temp.yml -o .
spanr compare --op union target.temp.yml query.temp.yml -o axt.union.yml
spanr stat chr.sizes axt.union.yml -o union.csv

# links by lastz-chain
fasops links axt.correct.fas -o stdout |
  perl -nl -e "s/(target|query)\.//g; print;" \
    >links.lastz.tsv

# remove species names
# remove duplicated sequences
# remove sequences with more than 250 Ns
fasops separate axt.correct.fas --nodash --rc -o stdout |
  perl -nl -e "/^>/ and s/^>(target|query)\./\>/; print;" |
  faops filter -u stdin stdout |
  faops filter -n 250 stdin stdout \
    >axt.gl.fasta

# Get more paralogs
egaz blastn \
  --parallel 12 \
  axt.gl.fasta genome.fa \
  -o axt.bg.blast 
egaz blastmatch \
  --parallel 12 -c 0.95 \
  axt.bg.blast \
  -o axt.bg.region
samtools faidx genome.fa -r axt.bg.region --continue |
  perl -p -e "/^>/ and s/:/(+):/" \
    >axt.bg.fasta

cat axt.gl.fasta axt.bg.fasta |
  faops filter -u stdin stdout |
  faops filter -n 250 stdin stdout \
    >axt.all.fasta

# Link paralogs
egaz blastn \
  --parallel 12 \
  axt.all.fasta axt.all.fasta \
  -o axt.all.blast
egaz blastlink \
  --parallel 12 -c 0.95 \
  axt.all.blast \
  -o links.blast.tsv

# Merge paralogs
# Sort links
linkr sort -o links.sort.tsv \
  links.lastz.tsv links.blast.tsv

# Clean links
linkr clean links.sort.tsv -o links.sort.clean.tsv
rgr merge links.sort.clean.tsv -o links.merge.tsv -c 0.95
linkr clean links.sort.clean.tsv -o links.clean.tsv -r links.merge.tsv --bundle 500

# Connect links
linkr connect links.clean.tsv -o links.connect.tsv -r 0.9
linkr filter links.connect.tsv -o links.filter.tsv -r 0.8
```

#### Create links

```shell
fasops create links.filter.tsv -o multi.temp.fas -g genome.fa
fasops refine multi.temp.fas -o multi.refine.fas --msa mafft -p 16 --chop 10
fasops links multi.refine.fas -o stdout |
  linkr sort stdin -o stdout |
  linkr filter stdin -n 2-50 -o links.refine.tsv

fasops links multi.refine.fas -o stdout --best |
  linkr sort stdin -o links.best.tsv
fasops create links.best.tsv -o pair.temp.fas -g genome.fa --name Atha
fasops refine pair.temp.fas -o pair.refine.fas --msa mafft -p 16

perl -nla -F"\t" -e "print for @F" <links.refine.tsv |
  spanr cover stdin -o cover.yml

echo "key,count" >links.count.csv
for n in 2 3 4-50; do
  linkr filter links.refine.tsv -n ${n} -o stdout \
    >links.copy${n}.tsv

  perl -nla -F"\t" -e "print for @F" <links.copy${n}.tsv |
    spanr cover stdin -o copy${n}.temp.yml

  wc -l links.copy${n}.tsv |
    perl -nl -e "
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

mkdir -p ../Results
cp cover.yml ../Results/Atha.cover.yml
cp copy.yml ../Results/Atha.copy.yml
mv cover.yml.csv ../Results/Atha.cover.csv
mv copy.csv ../Results/Atha.copy.csv
cp links.refine.tsv ../Results/Atha.links.tsv
mv multi.refine.fas ../Results/Atha.multi.fas
mv pair.refine.fas ../Results/Atha.pair.fas

find . -type f -name "*genome.fa*" | parallel --no-run-if-empty rm
find . -type f -name "*all.fasta*" | parallel --no-run-if-empty rm
find . -type f -name "*.sep.fasta" | parallel --no-run-if-empty rm
find . -type f -name "axt.*" | parallel --no-run-if-empty rm
find . -type f -name "replace.*.tsv" | parallel --no-run-if-empty rm
find . -type f -name "*.temp.yml" | parallel --no-run-if-empty rm
find . -type f -name "*.temp.fas" | parallel --no-run-if-empty rm

cd ..
```

#### Circos draw

```shell
mkdir Circos && cd Circos
cp ../../data/circos.conf .

# Generate karyotype files
if [ ! -e karyotype.Atha.txt ]; then
  perl -anl -e '$i++; print qq{chr - $F[0] $F[0] 0 $F[1] chr$i}' \
    /home/ivan/fat/MASED/Atha/chr.sizes \
    >karyotype.Atha.txt
fi

# Spaces among chromosomes
if [[ $(perl -n -e '$l++; END{print qq{$l\n}}' /home/ivan/fat/MASED/Atha/chr.sizes) > 1 ]]; then
  perl -nlpi -e 's/    default = 0r/    default = 0.005r/;' circos.conf
  perl -nlpi -e 's/show_label     = no/show_label     = yes/;' circos.conf
fi

# Chromosome units
SIZE=$(perl -an -F'\t' -e '$s += $F[1]; END{print qq{$s\n}}' /home/ivan/fat/MASED/Atha/chr.sizes)
if [ "${SIZE}" -ge 1000000000 ]; then
  echo "    * Set chromosome unit to 1 Mbp"
  perl -nlpi -e 's/chromosomes_units = 1000/chromosomes_units = 100000/;' circos.conf
elif [ "${SIZE}" -ge 100000000 ]; then
  echo "    * Set chromosome unit to 100 kbp"
  perl -nlpi -e 's/chromosomes_units = 1000/chromosomes_units = 100000/;' circos.conf
elif [ "${SIZE}" -ge 10000000 ]; then
  echo "    * Set chromosome unit to 10 kbp"
  perl -nlpi -e 's/chromosomes_units = 1000/chromosomes_units = 10000/;' circos.conf
else
  echo "    * Keep chromosome unit as 1 kbp"
fi

# Gff to highlight
if [ "${SIZE}" -ge 10000000 ]; then
  # avoid errors of too many highlights
  touch highlight.features.Atha.txt
  touch highlight.repeats.Atha.txt
else
  # coding and other features
  perl -anl -e '
        /^#/ and next;
        $F[0] =~ s/\.\d+//;
        $color = q{};
        $F[2] eq q{CDS} and $color = q{chr9};
        $F[2] eq q{ncRNA} and $color = q{dark2-8-qual-1};
        $F[2] eq q{rRNA} and $color = q{dark2-8-qual-2};
        $F[2] eq q{tRNA} and $color = q{dark2-8-qual-3};
        $F[2] eq q{tmRNA} and $color = q{dark2-8-qual-4};
        $color and ($F[4] - $F[3] > 49) and print qq{$F[0] $F[3] $F[4] fill_color=$color};
        ' \
    /home/ivan/fat/MASED/Atha/*.gff \
    >highlight.features.Atha.txt

  # repeats
  perl -anl -e '
        /^#/ and next;
        $F[0] =~ s/\.\d+//;
        $color = q{};
        $F[2] eq q{region} and $F[8] =~ /mobile_element|Transposon/i and $color = q{chr15};
        $F[2] =~ /repeat/ and $F[8] !~ /RNA/ and $color = q{chr15};
        $color and ($F[4] - $F[3] > 49) and print qq{$F[0] $F[3] $F[4] fill_color=$color};
        ' \
    /home/ivan/fat/MASED/Atha/*.gff \
    >highlight.repeats.Atha.txt
fi

# Links of paralog ranges
for n in 2 3 4-50; do
  linkr filter /home/ivan/fat/MASED/ATFSD/Results/Atha.links.tsv -n ${n} -o stdout \
    >links.copy${n}.tsv

  if [ "${n}" = "4-50" ]; then
    linkr circos links.copy${n}.tsv -o Atha.linkN.txt --highlight
  else
    linkr circos links.copy${n}.tsv -o Atha.link${n}.txt
  fi

  rm links.copy${n}.tsv
done

circos -noparanoid -conf circos.conf
#perl /home/linuxbrew/.linuxbrew/Cellar/circos/0.69-9/libexec/bin/circos -noparanoid -conf circos.conf
```

## Timeline identification of duplicated genes

```shell
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_protein.faa.gz -O ./data/Atha_raw_protein.fa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/255/GCF_000004255.2_v.1.0/GCF_000004255.2_v.1.0_protein.faa.gz -O ./data/Alyr_raw_protein.fa.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/375/325/GCF_000375325.1_Caprub1_0/GCF_000375325.1_Caprub1_0_protein.faa.gz -O ./data/Crub_raw_protein.fa.gz

# Brassica rapa GCF_000309985
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/309/985/GCF_000309985.2_CAAS_Brap_v3.01/GCF_000309985.2_CAAS_Brap_v3.01_protein.faa.gz -O ./data/Brap_raw_protein.fa.gz
curl -fsSL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/309/985/GCF_000309985.2_CAAS_Brap_v3.01/GCF_000309985.2_CAAS_Brap_v3.01_genomic.gff.gz |
  pigz -dc |
  awk '$1~"^NC_024"&&$3=="gene"{if($7=="+"){print $1 "\t" $9 "\t" $4 "\t" $5}
    else{print $1 "\t" $9 "\t" $5 "\t" $4}}' |
  perl -lane '$F[1] =~ /Name=(\w+)/;
	    $F[1] = $1;
	    print join("\t", @F)' \
    >./data/Brap.gff.tsv

# Carica papaya
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/535/GCF_000150535.2_Papaya1.0/GCF_000150535.2_Papaya1.0_protein.faa.gz -O ./data/Cpap_raw_protein.fa.gz
curl -fsSL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/535/GCF_000150535.2_Papaya1.0/GCF_000150535.2_Papaya1.0_genomic.gff.gz |
  pigz -dc |
    awk '$3=="gene"{if($7=="+"){print $1 "\t" $9 "\t" $4 "\t" $5}
    else{print $1 "\t" $9 "\t" $5 "\t" $4}}' |
  perl -lane '$F[1] =~ /Name=(\w+)/;
	    $F[1] = $1;
	    print join("\t", @F)' \
    >./data/Cpap.gff.tsv

# Theobroma cacao GCF_000208745
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/208/745/GCF_000208745.1_Criollo_cocoa_genome_V2/GCF_000208745.1_Criollo_cocoa_genome_V2_protein.faa.gz -O ./data/Tcac_raw_protein.fa.gz
curl -fsSL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/208/745/GCF_000208745.1_Criollo_cocoa_genome_V2/GCF_000208745.1_Criollo_cocoa_genome_V2_genomic.gff.gz |
  pigz -dc |
  awk '$1~"^NC_03085[0-9]"&&$3=="gene"{if($7=="+"){print $1 "\t" $9 "\t" $4 "\t" $5}
    else{print $1 "\t" $9 "\t" $5 "\t" $4}}' |
  perl -lane '$F[1] =~ /Name=(\w+)/;
	    $F[1] = $1;
	    print join("\t", @F)' \
    >./data/Tcac.gff.tsv


cd ./data
sed -i 's/Carub\./Carub/g' Crub/Crub.gff3
sed -i 's/Tc\([0-9]\{2\}\)v2_g/Tc\1v2p/g' Tcac/Tcac.gff3
sed -i 's/Carub\./Carub/g' Crub/Crub.pep
sed -i 's/Tc\([0-9]\{2\}\)v2_p/Tc\1v2p/g' Tcac/Tcac.pep

awk '$3 == "gene" {print $1 "\t" $4 "\t" $5 "\t" $9 "\t" $7}' Atha/Atha.gff3 > Atha.gene.gff
awk '$3 == "gene" {print $1 "\t" $4 "\t" $5 "\t" $9 "\t" $7}' Alyr/Alyr.gff3 > Alyr.gene.gff
awk '$3 == "gene" {print $1 "\t" $4 "\t" $5 "\t" $9 "\t" $7}' Crub/Crub.gff3 > Crub.gene.gff
awk '$3 == "gene" {print $1 "\t" $4 "\t" $5 "\t" $9 "\t" $7}' Brap/Brap.gff3 > Brap.gene.gff
awk '$3 == "gene" {print $1 "\t" $4 "\t" $5 "\t" $9 "\t" $7}' Tcac/Tcac.gff3 > Tcac.gene.gff

perl ../gff_pep.pl \
  --abbr "Atha" --pep_tag "Name=" \
  --pep ./Atha/Atha.pep --gff ./Atha.gene.gff \
  --out_gff ./AT.gff --out_pep ./AT.pep

perl ../gff_pep.pl \
  --abbr "Alyr" --pep_tag "Name=" \
  --pep ./Alyr/Alyr.pep --gff ./Alyr.gene.gff \
  --out_gff ./AL.gff --out_pep ./AL.pep
  
perl ../gff_pep.pl \
  --abbr "Crub" --pep_tag "Name=" \
  --pep ./Crub/Crub.pep --gff ./Crub.gene.gff \
  --out_gff ./CR.gff --out_pep ./CR.pep

perl ../gff_pep.pl \
  --abbr "Brap" --pep_tag "ID=" \
  --pep ./Brap/Brap.pep --gff ./Brap.gene.gff \
  --out_gff ./BR.gff --out_pep ./BR.pep

perl ../gff_pep.pl \
  --abbr "Tcac" --pep_tag "Name=" \
  --pep ./Tcac/Tcac.pep --gff ./Tcac.gene.gff \
  --out_gff ./TC.gff --out_pep ./TC.pep

makeblastdb -in AT.pep -dbtype prot -parse_seqids -out ATdb
blastp -query AL.pep -db ATdb -out AT_AL.blast -evalue 1e-10 -num_threads 12 -outfmt 6 -num_alignments 5
blastp -query AL.pep -db ATdb -out AT_AL.blast -evalue 1e-10 -num_threads 12 -outfmt 6 -num_alignments 5
blastp -query CR.pep -db ATdb -out AT_CR.blast -evalue 1e-10 -num_threads 12 -outfmt 6 -num_alignments 5
blastp -query BR.pep -db ATdb -out AT_BR.blast -evalue 1e-10 -num_threads 12 -outfmt 6 -num_alignments 5
blastp -query TC.pep -db ATdb -out AT_TC.blast -evalue 1e-10 -num_threads 12 -outfmt 6 -num_alignments 5

```
