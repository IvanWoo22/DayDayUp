# Catch and count reads which stopped at target sites.
### It is a short perl script that can normalize and accurately locate given reads at each sites of a given geneome data, and then generate the outcome in a new file.
## Given data.
We have a reads data with high quality:
```bash
head clean.txt
####################################################
# 23	1544686	AGACGTGGCGACCCGCTGAATTT
# 20	546896	CGCGACCTCAGATCAGACGT
# 24	487044	GGAGGTCCGTAGCGGTCCTGACGT
# 38	418798	ACGGCCCTGGCGGAGCGCTGAGAAGACGGTCGAACTTG
# 23	314772	AGACGTGGCGACCCGCTGAATTC
# 31	292285	AAGAACGAAAGTCGGAGGTTCGAAGACGATC
# 35	247034	AAGAACGAAAGTCGGAGGTTCGAAGACGATCAGAT
# 20	238535	CGCGACCTCAGATCAGACGC
# 24	222855	GGAGGTCCGTAGCGGTCCTGACGC
# 26	206660	GGAGGTCCGTAGCGGTCCTGACGTGC
#####################################################
```
The first row is the length of reads, the sencond is count number, the third is the reads sequence.  
And, the gene sequence we want to map to is `18s.txt` :
```bash
TACCTGGTTGATCCTGC......GTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA
```
## Analyse.
```bash
perl stopcounts.pl \
  18s.txt \
  clean.txt \
  output.txt
#Here output.txt is named by yourself.
```

## Output.
```bash
head output.txt
#####################################################
# TACCTGGTTGATCCT	2
# ACCTGGTTGATCCTG	3
# CCTGGTTGATCCTGC	No Match.
# CTGGTTGATCCTGCC	345
# TGGTTGATCCTGCCA	No Match.
# GGTTGATCCTGCCAG	35
# GTTGATCCTGCCAGT	809
# TTGATCCTGCCAGTA	2
# TGATCCTGCCAGTAG	50
# GATCCTGCCAGTAGC	2974
#####################################################
# Note:This is not the result of our given data.
```
Here the first row is the sequence fragment and the second row is the counts of reads ended at this site.
