# 例如，提取基因区域（gene）、外显子（exon）、内含子（intron）等
awk '$3=="gene"' gencode.v*.annotation.gff3 > genes.gff
awk '$3=="exon"' gencode.v*.annotation.gff3 > exons.gff
awk '$3=="CDS"' gencode.v*.annotation.gff3 > cds.gff
awk '$3=="five_prime_UTR"' gencode.v*.annotation.gff3 > utr5.gff
awk '$3=="three_prime_UTR"' gencode.v*.annotation.gff3 > utr3.gff

bedtools intersect -a sites.total.bed -b genes.gff -wa -wb > sites_in_genes.bed
bedtools intersect -a sites.total.bed -b exons.gff -wa -wb > sites_in_exons.bed
bedtools intersect -a sites.total.bed -b cds.gff -wa -wb > sites_in_cds.bed
bedtools intersect -a sites.total.bed -b utr5.gff -wa -wb > sites_in_utr5.bed
bedtools intersect -a sites.total.bed -b utr3.gff -wa -wb > sites_in_utr3.bed


echo "Gene Overlap: $(wc -l < sites_in_genes.bed)"
echo "Exon Overlap: $(wc -l < sites_in_exons.bed)"
echo "CDS Overlap: $(wc -l < sites_in_cds.bed)"


