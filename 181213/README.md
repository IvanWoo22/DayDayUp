# Find which gene the CpG island locate at.
### Using this Perl script, you can get the information which gene each CpG island is nearest to.
## Example.
First, we have a processed `gff3` file `chr1.gff3`:
```
1	.	biological_region	3003629	3003751	0.629	-	.	external_name=rank %3D 1;logic_name=firstef
1	NCBI	biological_region	3009957	3010783	.	+	.	external_name=CH29-168G11 (start);logic_name=mouse_ch29_clones
1	NCBI	biological_region	3009969	3010627	.	+	.	external_name=CH29-168H10 (start);logic_name=mouse_ch29_clones
1	NCBI	biological_region	3028057	3028936	.	+	.	external_name=DN-332C11 (end);logic_name=mouse_dn_clones
1	NCBI	biological_region	3031055	3031498	.	+	.	external_name=WI1-408L1 (end);logic_name=mouse_wi1_clones
1	NCBI	biological_region	3045372	3046053	.	+	.	external_name=CH29-573J4 (start);logic_name=mouse_ch29_clones
1	NCBI	biological_region	3049732	3050154	.	+	.	external_name=RP24-86G11 (start);logic_name=mouse_rp24_clones
1	NCBI	biological_region	3051577	3052653	.	+	.	external_name=B6Ng01-203D22 (start);logic_name=mouse_b6ng01_clones
1	NCBI	biological_region	3054488	3055154	.	+	.	external_name=MHPP-378N8 (end);logic_name=mouse_mhpp_clones
1	NCBI	biological_region	3054488	3055174	.	+	.	external_name=MHPP-378N7 (end);logic_name=mouse_mhpp_clones
```
And then we will get a file `chr1.txt` containing information of each CpG island after running DSS in R:
```
134298875	134299566	692	114	0.145614472815011	0.349156949277884	-0.203542476462872	-867.202994990791
63309370	63309825	456	64	0.196242493673155	0.438249471515221	-0.242006977842066	-467.387266238766
22533946	22534344	399	42	0.074835970813796	0.241678057971851	-0.166842087158055	-401.127016504278
190172744	190173199	456	64	0.274568229706463	0.474962792381054	-0.200394562674591	-382.336657536679
75360328	75360554	227	45	0.124638692994966	0.290299915593009	-0.165661222598043	-316.217688991217
152090379	152090588	210	43	0.047602023748424	0.208987737960803	-0.161385714212379	-304.662431024946
136789097	136789407	311	46	0.549672747263946	0.318564902027382	0.231107845236564	294.953821684057
119255397	119255922	526	39	0.209078923034721	0.520616146335139	-0.311537223300418	-287.030936428374
136740776	136741735	960	44	0.43126136766852	0.70498159706238	-0.273720229393859	-271.509172195656
91405818	91406338	521	41	0.180960449016958	0.406682744927556	-0.225722295910598	-250.334997633041
```
To call the nearest gene for each island, we just run:
```bash
perl CpG2gene.pl chr1.gff3 chr1.txt chr1.loc
```
NOTED: It is up to you to decide how to name the output file. Here I just use `chr1.loc`
## Output explaination.
```bash
head chr1.loc
```
chr1	134296783-134332928	ID=gene:ENSMUSG00000026458;Name=Ppfia4;biotype=protein_coding;description=protein tyrosine phosphatase%2C receptor type%2C f polypeptide (PTPRF)%2C interacting protein (liprin)%2C alpha 4 [Source:MGI Symbol%3BAcc:MGI:1915757];gene_id=ENSMUSG00000026458;logic_name=ensembl_havana_gene;version=13	134298875-134299566	114/692	0.145614472815011	0.349156949277884	-0.203542476462872	-867.202994990791
chr1	63273265-63314576	ID=gene:ENSMUSG00000027520;Name=Zdbf2;biotype=protein_coding;description=zinc finger%2C DBF-type containing 2 [Source:MGI Symbol%3BAcc:MGI:1921134];gene_id=ENSMUSG00000027520;logic_name=ensembl_havana_gene;version=15	63309370-63309825	64/456	0.196242493673155	0.438249471515221	-0.242006977842066	-467.387266238766
chr1	22286251-22805994	ID=gene:ENSMUSG00000041670;Name=Rims1;biotype=protein_coding;description=regulating synaptic membrane exocytosis 1 [Source:MGI Symbol%3BAcc:MGI:2152971];gene_id=ENSMUSG00000041670;logic_name=ensembl_havana_gene;version=16	22533946-22534344	42/399	0.074835970813796	0.241678057971851	-0.166842087158055	-401.127016504278
chr1	190170296-190217276	ID=gene:ENSMUSG00000079045;Name=Prox1os;biotype=antisense;description=prospero homeobox 1%2C opposite strand [Source:MGI Symbol%3BAcc:MGI:4937200];gene_id=ENSMUSG00000079045;logic_name=havana;version=3	190172744-190173199	64/456	0.274568229706463	0.474962792381054	-0.200394562674591	-382.336657536679
chr1	75360329-75368579	ID=gene:ENSMUSG00000026208;Name=Des;biotype=protein_coding;description=desmin [Source:MGI Symbol%3BAcc:MGI:94885];gene_id=ENSMUSG00000026208;logic_name=ensembl_havana_gene;version=9	75360328-75360554	45/227	0.124638692994966	0.290299915593009	-0.165661222598043	-316.217688991217
chr1	152089617-152090619	external_name=oe %3D 0.76;logic_name=cpg	152090379-152090588	43/210	0.047602023748424	0.208987737960803	-0.161385714212379	-304.662431024946
chr1	136789140-136789568	external_name=oe %3D 0.69;logic_name=cpg	136789097-136789407	46/311	0.549672747263946	0.318564902027382	0.231107845236564	294.953821684057
chr1	119255241-119255882	external_name=BMQ-57M24 (end);logic_name=mouse_bmq_clones	119255397-119255922	39/526	0.209078923034721	0.520616146335139	-0.311537223300418	-287.030936428374
chr1	136740484-136741160	external_name=RP24-270D13 (start);logic_name=mouse_rp24_clones	136740776-136741735	44/960	0.43126136766852	0.70498159706238	-0.273720229393859	-271.509172195656
chr1	91405906-91406143	external_name=rank %3D 1;logic_name=firstef	91405818-91406338	41/521	0.180960449016958	0.406682744927556	-0.225722295910598	-250.334997633041
```
`chromosome gene_range  gene_info CpG_range nCG/length meanMethy1 meanMethy2  diff.Methy  areaStat`
