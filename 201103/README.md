```
pigz -dc gencode.v34.annotation.gff3.gz |
    grep "protein_coding" |
    awk '$3=="start_codon"{print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
    perl format.pl > start_codon.tsv

pigz -dc gencode.v34.annotation.gff3.gz |
    grep "protein_coding" |
    awk '$3=="stop_codon"{print $1 "\t" $4 "\t" $5 "\t" $7 "\t" $9}' |
    perl format.pl > stop_codon_format.tsv

perl -ne 's/\n//g;
    ($t,$g)=split"\t";
    ($t_n,undef)=split/\./,$t;
    ($g_n,undef)=split/\./,$g;
    print "$t_n\t$g_n\n"
    ' <~/20200901/data/mmu_main_transcript.txt \
    >mmu_main_transcript.txt

perl -ne 's/\n//g;
    ($t,$g)=split"\t";
    ($t_n,undef)=split/\./,$t;
    ($g_n,undef)=split/\./,$g;
    print "$t_n\t$g_n\n"
    ' <~/20200901/data/mmu_represent_transcript.txt \
    >mmu_represent_transcript.txt

for i in `cut -f1 mmu_main_transcript.txt`; do
    grep "${i}" mmu_start_codon_format.tsv \
        >>mmu_start_codon_main.tsv;
done

for i in `cut -f1 mmu_main_transcript.txt`; do
    grep "${i}" mmu_stop_codon_format.tsv \
        >>mmu_stop_codon_main.tsv;
done

for i in `cut -f1 mmu_represent_transcript.txt`; do
    grep "${i}" mmu_start_codon_format.tsv \
        >>mmu_start_codon_represent.tsv;
done

for i in `cut -f1 mmu_represent_transcript.txt`; do
    grep "${i}" mmu_stop_codon_format.tsv \
        >>mmu_stop_codon_represent.tsv;
done

perl common.pl \
    mmu_start_codon_represent.tsv \
    mmu_start_codon_main.tsv \
    >mmu_start_codon.tsv

perl common.pl \
    mmu_stop_codon_represent.tsv \
    mmu_stop_codon_main.tsv \
    >mmu_stop_codon.tsv

pigz -dc ~/20200901/data/mmu.gff3.gz |
    grep "protein_coding" |
    awk '$3=="transcript"||$3=="exon"' |
    perl pickup_codon.pl \
        mmu_start_codon.tsv \
        >mmu_transcript_wsc.yml

pigz -dc ~/20200901/data/mmu.gff3.gz |
    grep "protein_coding" |
    awk '$3=="transcript"||$3=="exon"' |
    perl pickup_codon.pl \
        mmu_stop_codon.tsv \
        >mmu_transcript_wstopc.yml


```