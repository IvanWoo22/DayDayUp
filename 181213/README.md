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
```
perl CpG2gene.pl chr1.gff3 chr1.txt chr1.loc
```
NOTED: It is up to you to decide how to name the output file. Here I just use `chr1.loc`
