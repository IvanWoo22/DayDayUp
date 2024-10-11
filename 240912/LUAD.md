# Steps for LUAD_meth

<!-- toc.levels="1..3" -->

- [Steps for LUAD_meth](#steps-for-luad_meth)
- [Downloads](#downloads)
	- [Manifest](#manifest)
	- [gdc-client](#gdc-client)
	- [Clinical info](#clinical-info)
	- [TIL map](#til-map)
- [Merging samples](#merging-samples)
	- [gene_expression](#gene_expression)
	- [methylation_beta_value](#methylation_beta_value)
	- [mirna_expression](#mirna_expression)
- [Extra validation datasets](#extra-validation-datasets)
	- [GSE39279](#gse39279)
	- [GSE66836](#gse66836)
	- [GSE119144](#gse119144)
	- [XH33](#xh33)
- [Preparations](#preparations)
	- [Clinical features](#clinical-features)
- [To and From HPC](#to-and-from-hpc)
- [1_training](#1_training)
- [1_bootstrap](#1_bootstrap)
- [1_testing](#1_testing)
- [2_training](#2_training)
- [2_bootstrap](#2_bootstrap)
- [2_testing](#2_testing)
- [3_training](#3_training)
- [3_bootstrap](#3_bootstrap)
- [3_testing](#3_testing)
- [3_non_rare](#3_non_rare)
- [Clinical groups](#clinical-groups)
- [3_clinical](#3_clinical)
- [Extra validations](#extra-validations)
	- [GSE39279 LUAD](#gse39279-luad)
	- [GSE119144 PFS](#gse119144-pfs)
	- [GSE119144 benefit](#gse119144-benefit)
	- [XH33](#xh33-1)
- [3_validating](#3_validating)
- [4_non_rare](#4_non_rare)
- [4_combination](#4_combination)
- [4_figure](#4_figure)
- [4_figure_GSE119144](#4_figure_gse119144)
- [Top 300](#top-300)
- [Interesting loci](#interesting-loci)
- [4_loci](#4_loci)
- [Pack up](#pack-up)

# Downloads

[TCGA-LUAD](https://portal.gdc.cancer.gov/projects/TCGA-LUAD)

## gencode

```shell
mkdir -p LUAD
cd LUAD || exit

GENCODE_VERSION=46
aria2c -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_"${GENCODE_VERSION}"/gencode.v"${GENCODE_VERSION}".annotation.gtf.gz

gzip -dcf gencode.v"${GENCODE_VERSION}".annotation.gtf.gz \
	| grep -v "^#" \
	| tsv-summarize --group-by 3 --count
#gene	63086
#transcript	254070
#exon	1668627
#CDS	899146
#start_codon	98995
#stop_codon	92909
#UTR	390193
#Selenocysteine	130

echo -e '#ENSGID\tSYMBOL\tLoc\tGENETYPE' >gencode.tsv
gzip -dcf gencode.v"${GENCODE_VERSION}".annotation.gtf.gz \
	| grep -v "^#" \
	| tsv-filter --str-eq 3:gene \
	| perl -nla -F"\t" -e '
        my $chr = $F[0];
        my $start = $F[3];
        my $end = $F[4];
        my $strand = $F[6];
        my $loc = qq{$chr($strand):$start-$end};
        my $anno = $F[8];
        my %anno_of;
        for my $part ( split(/;\s*/, $anno) ) {
            my ($key, $value) =  split(/ /, $part);
            $value =~ s/"//g;
            $anno_of{$key} = $value;
        }
        print join qq{\t},
            $anno_of{gene_id}, $anno_of{gene_name}, $loc, $anno_of{gene_type};
    ' \
		>>gencode.tsv
keep-header gencode.tsv -- datamash check
#63086 lines, 4 fields
```

## Manifest

```r
library(GenomicDataCommons)
library(readr)
library(stringr)

project_id <- 'TCGA-LUAD'

# counts
resp <- cases() |>
  filter(project.project_id == project_id) |>
  facet('samples.sample_type') |>
  aggregations()
knitr::kable(resp$samples.sample_type[, c("key", "doc_count")], format = "markdown")

resp <- files() |>
  filter(cases.project.project_id == project_id) |>
  filter(access == 'open') |>
  facet('type') |>
  aggregations()
knitr::kable(resp$type[, c("key", "doc_count")], format = "markdown")

resp <- files() |>
  filter(cases.project.project_id == project_id) |>
  filter(access == 'open') |>
  facet('platform') |>
  aggregations()
knitr::kable(resp$platform[, c("key", "doc_count")], format = "markdown")

resp <- files() |>
  filter(cases.project.project_id == project_id) |>
  filter(access == 'open') |>
  facet('analysis.workflow_type') |>
  aggregations()
knitr::kable(resp$analysis.workflow_type[, c("key", "doc_count")], format = "markdown")

# All open access files
manifest <- files() |>
  filter(cases.project.project_id == project_id) |>
  filter(access == 'open') |>
  manifest()
write_tsv(manifest, file = "gdc_manifest.txt")

# All files
records <- files() |>
  filter(cases.project.project_id == project_id) |>
  filter(access == 'open') |>
  GenomicDataCommons::select(
    c(default_fields(files()), 'type', 'cases.submitter_id', 'cases.samples.sample_type')
  ) |>
  response_all() |>
  results()
write_lines(
  str_flatten(c("filename", "data_category", "data_format"), collapse = "\t"),
  file = "biotab.tsv"
)

sink('fileinfo.tsv')
cat("filename\ttype\tcase_id\tsample_type\tdata_type\tdata_format\tdata_category\tplatform\texperimental_strategy\n")

for (i in 1:length(records$file_name)) {
  case_id <- records$cases[i][[1]]$submitter_id
  # clinical info
  if (length(case_id) > 1) {
    fields <- c(
      records$file_name[i],
      records$data_category[i],
      records$data_format[i]
    )
    write_lines(
      str_flatten(fields, collapse = "\t"),
      file = "biotab.tsv",
      append = TRUE
    )
    next
  }
  # flatten multiple sample types
  sample_type <- str_flatten(records$cases[i][[1]]$samples[[1]]$sample_type, collapse = "|")
  fields <- c(
    records$file_name[i],
    records$type[i],
    case_id,
    sample_type,
    records$data_type[i],
    records$data_format[i],
    records$data_category[i],
    records$platform[i],
    records$experimental_strategy[i]
  )
  fields <- str_replace_na(fields, replacement = "NA")
  cat(
    paste(
      str_flatten(fields, collapse = "\t"),
      "\n"
    )
  )
}
sink()
```

| key                  | doc_count |
|:---------------------|----------:|
| primary tumor        |       585 |
| blood derived normal |       424 |
| solid tissue normal  |       274 |
| ffpe scrolls         |         2 |
| recurrent tumor      |         2 |

| key                      | doc_count |
|:-------------------------|----------:|
| copy_number_segment      |      3343 |
| slide_image              |      1608 |
| copy_number_estimate     |      1556 |
| masked_methylation_array |      1420 |
| mirna_expression         |      1134 |
| biospecimen_supplement   |      1123 |
| methylation_beta_value   |       710 |
| clinical_supplement      |       623 |
| masked_somatic_mutation  |       618 |
| gene_expression          |       600 |
| pathology_report         |       523 |
| protein_expression       |       365 |

| key                            | doc_count |
|:-------------------------------|----------:|
| affymetrix snp 6.0             |      4895 |
| illumina                       |      2356 |
| illumina human methylation 450 |      1521 |
| illumina human methylation 27  |       450 |
| rppa                           |       365 |
| illumina methylation epic v2   |       159 |
| _missing                       |      3877 |

| key                                                  | doc_count |
|:-----------------------------------------------------|----------:|
| DNAcopy                                              |      2294 |
| SeSAMe Methylation Beta Estimation                   |      2130 |
| BCGSC miRNA Profiling                                |      1134 |
| ASCAT2                                               |      1088 |
| ASCAT3                                               |      1006 |
| Aliquot Ensemble Somatic Variant Merging and Masking |       618 |
| STAR - Counts                                        |       600 |
| ABSOLUTE LiftOver                                    |       507 |
| AscatNGS                                             |         4 |
| _missing                                             |      4242 |

## File count

```shell
tsv-sort -k2,2 -k1,1 <biotab.tsv >tmp.txt
mv tmp.txt biotab.tsv

datamash check <fileinfo.tsv
# 13602 lines, 9 fields

tsv-summarize -H \
-g data_category,type,data_type --count \
<fileinfo.tsv \
	| keep-header -- tsv-sort -k1,1 -k2,2 \
	| mlr --itsv --omd cat

tsv-filter -H --str-not-in-fld "sample_type:|" \
<fileinfo.tsv  \
	| tsv-summarize -H -g sample_type,type --count \
	| keep-header -- tsv-sort \
	| mlr --itsv --omd cat
```

| data_category               | type                     | data_type                           | count |
|-----------------------------|--------------------------|-------------------------------------|-------|
| Biospecimen                 | biospecimen_supplement   | Biospecimen Supplement              | 1107  |
| Biospecimen                 | slide_image              | Slide Image                         | 1608  |
| Clinical                    | clinical_supplement      | Clinical Supplement                 | 617   |
| Clinical                    | pathology_report         | Pathology Report                    | 523   |
| Copy Number Variation       | copy_number_estimate     | Gene Level Copy Number              | 1556  |
| Copy Number Variation       | copy_number_segment      | Allele-specific Copy Number Segment | 1047  |
| Copy Number Variation       | copy_number_segment      | Copy Number Segment                 | 1149  |
| Copy Number Variation       | copy_number_segment      | Masked Copy Number Segment          | 1147  |
| DNA Methylation             | masked_methylation_array | Masked Intensities                  | 1420  |
| DNA Methylation             | methylation_beta_value   | Methylation Beta Value              | 710   |
| Proteome Profiling          | protein_expression       | Protein Expression Quantification   | 365   |
| Simple Nucleotide Variation | masked_somatic_mutation  | Masked Somatic Mutation             | 618   |
| Transcriptome Profiling     | gene_expression          | Gene Expression Quantification      | 600   |
| Transcriptome Profiling     | mirna_expression         | Isoform Expression Quantification   | 567   |
| Transcriptome Profiling     | mirna_expression         | miRNA Expression Quantification     | 567   |

| sample_type          | type                     | count |
|----------------------|--------------------------|-------|
|                      | biospecimen_supplement   | 1107  |
| Blood Derived Normal | copy_number_segment      | 826   |
|                      | clinical_supplement      | 617   |
| Primary Tumor        | copy_number_estimate     | 505   |
| Primary Tumor        | copy_number_segment      | 1108  |
| Primary Tumor        | gene_expression          | 539   |
| Primary Tumor        | masked_methylation_array | 1304  |
| Primary Tumor        | methylation_beta_value   | 652   |
| Primary Tumor        | mirna_expression         | 1038  |
| Primary Tumor        | pathology_report         | 522   |
| Primary Tumor        | protein_expression       | 365   |
| Primary Tumor        | slide_image              | 1359  |
| Recurrent Tumor      | copy_number_estimate     | 2     |
| Recurrent Tumor      | copy_number_segment      | 4     |
| Recurrent Tumor      | gene_expression          | 2     |
| Recurrent Tumor      | masked_methylation_array | 4     |
| Recurrent Tumor      | methylation_beta_value   | 2     |
| Recurrent Tumor      | mirna_expression         | 4     |
| Recurrent Tumor      | pathology_report         | 1     |
| Recurrent Tumor      | slide_image              | 5     |
| Solid Tissue Normal  | copy_number_segment      | 356   |
| Solid Tissue Normal  | gene_expression          | 59    |
| Solid Tissue Normal  | masked_methylation_array | 112   |
| Solid Tissue Normal  | methylation_beta_value   | 56    |
| Solid Tissue Normal  | mirna_expression         | 92    |
| Solid Tissue Normal  | slide_image              | 244   |

## gdc-client

```shell
for type in \
	copy_number_estimate \
	gene_expression \
	masked_somatic_mutation \
	methylation_beta_value \
	mirna_expression \
	protein_expression \
	slide_image; do
	tsv-join -H \
		--data-fields "file_name" \
		-f <(tsv-filter -H --str-eq "type:${type}" <fileinfo.tsv) \
		--key-fields "filename" \
		<gdc_manifest.txt \
		| keep-header -- sort \
			>"gdc_manifest.${type}.txt"
done

wc -l gdc_manifest.*.txt
#   1553 gdc_manifest.copy_number_estimate.txt
#    600 gdc_manifest.gene_expression.txt
#    618 gdc_manifest.masked_somatic_mutation.txt
#    709 gdc_manifest.methylation_beta_value.txt
#   1133 gdc_manifest.mirna_expression.txt
#    365 gdc_manifest.protein_expression.txt
#   1605 gdc_manifest.slide_image.txt

mkdir -p copy_number_estimate
gdc-client download -n 4 -m gdc_manifest.copy_number_estimate.txt -d copy_number_estimate

mkdir -p gene_expression
gdc-client download -n 4 -m gdc_manifest.gene_expression.txt -d gene_expression

mkdir -p masked_somatic_mutation
gdc-client download -n 4 -m gdc_manifest.masked_somatic_mutation.txt -d masked_somatic_mutation

mkdir -p methylation_beta_value
gdc-client download -n 4 -m gdc_manifest.methylation_beta_value.txt -d methylation_beta_value

mkdir -p mirna_expression
gdc-client download -n 4 -m gdc_manifest.mirna_expression.txt -d mirna_expression

mkdir -p protein_expression
gdc-client download -n 4 -m gdc_manifest.protein_expression.txt -d protein_expression

mkdir -p slide_image
gdc-client download -n 4 -m gdc_manifest.slide_image.txt -d slide_image

find copy_number_estimate -name "*.gene_level_copy_number*tsv" | wc -l
find gene_expression -name "*.augmented_star_gene_counts.tsv" | wc -l
find masked_somatic_mutation -name "*.aliquot_ensemble_masked.maf.gz" | wc -l
find methylation_beta_value -name "*.sesame.level3betas.txt" | wc -l
find mirna_expression -name "*.mirnas.quantification.txt" | wc -l
find mirna_expression -name "*.isoforms.quantification.txt" | wc -l
find slide_image -name "*.svs" | wc -l

```

## Clinical info

* nationwidechildrens.org_clinical_patient_luad.txt: 42bf5eb2-bc49-45be-b18a-290f712b006c
* nationwidechildrens.org_clinical_follow_up_v1.0_luad.txt: c0cb0e6a-b19e-46a5-b750-91ff91385071

* nationwidechildrens.org_clinical_drug_luad.txt: 59015dfd-045b-4a80-9b4a-3cc3a41707f0
* nationwidechildrens.org_clinical_radiation_luad.txt: 27bafc9d-b95c-4b49-b3e3-1c90bb1c15dd

```shell
cd rawdata/LUAD || exit

for category in \
	Biospecimen \
	Clinical; do
	tsv-join -H \
		--data-fields filename \
		-f <(
			tsv-filter -H --str-eq "data_category:${category}" <biotab.tsv
		) \
		--key-fields filename \
		<gdc_manifest.txt \
		| keep-header -- sort \
			>"gdc_manifest.${category}.txt"
done

mkdir -p Biospecimen
gdc-client download -n 4 -m gdc_manifest.Biospecimen.txt -d Biospecimen

mkdir -p Clinical
gdc-client download -n 4 -m gdc_manifest.Clinical.txt -d Clinical

tsv-append -H \
	<(
		find Clinical -name "*_clinical_patient_*.txt" \
			| xargs cat \
			| sed -e '2,3d' \
			| mlr --itsv --otsv cut -o -f bcr_patient_barcode,vital_status,last_contact_days_to,death_days_to
	) \
	<(
		find Clinical -name "*_clinical_follow_up_*.txt" \
			| xargs cat \
			| sed -e '2,3d' \
			| mlr --itsv --otsv cut -o -f bcr_patient_barcode,vital_status,last_contact_days_to,death_days_to
	) \
	| sed '/\[Completed\]/d' \
	| sed '/\[Discrepancy\]/d' \
	| sed 's/\[Not Available\]//g' \
	| sed 's/\[Not Applicable\]//g' \
	| perl -nlp -e '
        s/(Alive|Dead)\t\t(\d+)$/\1\t\2/g and next;
        s/(Alive|Dead)\t(\d+)\t$/\1\t\2/g and next;
        s/(Alive|Dead)\t(\d+)\t(\d+)$/\1\t\3/g and next;
    ' \
	| sed 's/Alive\t/0\t/' \
	| sed 's/Dead\t/1\t/' \
	| perl -nla -e '$F[2] < 0 and next; print' \
	| sed -e '1d' \
	| perl -nla -e '@F == 3 or next; print' \
	| (echo -e '#sample\tstatus\ttime' && cat) \
	| tsv-select -f 1,3,2 \
	| keep-header -- sort -k1,1 -k2,2nr \
	|
	# latest follow up
	tsv-uniq -H -f 1 \
		>basic_clinical.tsv

sed '1d' basic_clinical.tsv | datamash check
#513 lines, 3 fields

find Clinical -name "*_clinical_patient_*.txt" \
	| xargs head -n 1 \
	| tr "\t" "\n" \
	| tsv-uniq -r
#histologic_diagnosis

find Clinical -name "*_clinical_patient_*.txt" \
	| xargs cat \
	| sed '2,3d' \
	| tsv-select -H -f \
		bcr_patient_barcode,gender,ajcc_pathologic_tumor_stage,age_at_initial_pathologic_diagnosis,history_other_malignancy,tobacco_smoking_pack_years_smoked,histologic_diagnosis,location_lung_parenchyma,kras_mutation_found,kras_mutation_identified_type,egfr_mutation_status,egfr_mutation_identified_type \
	| tsv-select -H -e 7 \
	|
	# duplicated headers
	sed 's/\[Not Available\]//g' \
	| sed 's/\[Unknown\]//g' \
		>detailed_clinical.tsv

sed '1d' detailed_clinical.tsv | datamash check
#522 lines, 12 fields

find Clinical -name "*_clinical_drug_*.txt" \
	| xargs cat \
	| sed -e '2,3d' \
	| tsv-select -H -f bcr_patient_barcode,pharmaceutical_therapy_drug_name,pharmaceutical_therapy_type \
	| perl -nla -F'\t' -e '$F[1] = lc $F[1]; print join qq{\t}, @F' \
	| sed 's/platinum/platin/' \
	| sed 's/almita/alimta/' \
		>drug.tsv

tsv-summarize -H --unique-count 1 <drug.tsv
#179

find Clinical -name "*_clinical_radiation_*.txt" \
	| xargs cat \
	| sed -e '2,3d' \
	| tsv-select -H -f bcr_patient_barcode,radiation_therapy_type,radiation_therapy_site \
	| tsv-filter --ff-istr-ne 2:3 \
	|
	# [Not Available]	[Not Available]
	perl -nla -F'\t' -e '$F[1] = lc $F[1]; print join qq{\t}, @F' \
		>radiation.tsv

tsv-summarize -H --unique-count 1 <radiation.tsv
#98

for f in \
	basic_clinical.tsv \
	detailed_clinical.tsv \
	drug.tsv \
	radiation.tsv; do
	cat ${f} \
		| keep-header -- sort >tmp.txt
	mv tmp.txt ${f}
done
```

## TIL map

Tumor-Infiltrating Lymphocytes

```shell
cd rawdata/LUAD || exit

curl -L https://api.gdc.cancer.gov/data/08096f8f-7b56-495a-be45-62d5a56f2ee8 -o TILMap_TableS1.xlsx

plotr xlsx TILMap_TableS1.xlsx --sheet 'TILMap_TableS1.txt' -o stdout \
	| tsv-filter -H --str-eq Study:LUAD \
		>TIL.tsv

```

# Merging samples

## gene_expression

```shell
cd ~/data/cancer/rawdata/LUAD

find gene_expression -type f -name "*.FPKM.txt.gz" |
    head -n 1 |
    xargs gzip -dcf |
    cut -f 1 |
    (echo "#sample" && cat) |
    datamash transpose \
    > gene.tmp

cat fileinfo.tsv |
    tsv-filter -H --str-eq 2:gene_expression --str-in-fld 4:"Primary Tumor" |
    keep-header -- sort -k3,3 -k4,4r |
    tsv-uniq -H -f 3 |
    tsv-select -H -f 1,3 |
    grep -v "^#" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo {2} 1>&2
        filepath=$(find gene_expression -type f -name "{1}")
        gzip -dcf ${filepath} |
            tsv-select -H -f 1,2 |
            (echo -e '\''#sample\t{2}'\'' && cat) |
            datamash transpose |
            grep -v "^#"
    ' |
    perl -nla -F'\t' -e '
        for my $i ( 1 .. $#F) {
            if ($F[$i] ne q{NA}) {
                $F[$i] = sprintf qq{%.6f}, $F[$i];
            }
        }
        print join qq{\t}, @F;
    ' |
    sort \
    >> gene.tmp

sed '1d' gene.tmp | datamash check
#513 lines, 60484 fields

tsv-join \
    basic_clinical.tsv \
    --data-fields 1 \
    -f gene.tmp \
    --key-fields 1 \
    --append-fields 2-60484 \
    > gene.tsv

sed '1d' gene.tsv | datamash check
#504 lines, 60486 fields

cat gene.tsv |
    head -n 5 |
    tsv-select -f 1-5

rm gene.tmp

```

## methylation_beta_value

```shell
cd rawdata/LUAD || exit

find gene_expression -type f -name "*.FPKM.txt.gz" \
	| head -n 1 \
	| xargs gzip -dcf \
	| cut -f 1 \
	| (echo "#sample" && cat) \
	| datamash transpose \
		>gene.tmp

cat fileinfo.tsv \
	| tsv-filter -H --str-eq 2:gene_expression --str-in-fld 4:"Primary Tumor" \
	| keep-header -- sort -k3,3 -k4,4r \
	| tsv-uniq -H -f 3 \
	| tsv-select -H -f 1,3 \
	| grep -v "^#" \
	| parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
		echo {2} 1>&2
		filepath=$(find gene_expression -type f -name "{1}")
		gzip -dcf ${filepath} |
				tsv-select -H -f 1,2 |
				(echo -e '\''#sample\t{2}'\'' && cat) |
				datamash transpose |
				grep -v "^#"
    ' \
	| perl -nla -F'\t' -e '
        for my $i ( 1 .. $#F) {
            if ($F[$i] ne q{NA}) {
                $F[$i] = sprintf qq{%.6f}, $F[$i];
            }
        }
        print join qq{\t}, @F;
    ' \
	| sort \
		>>gene.tmp

sed '1d' gene.tmp | datamash check
#513 lines, 60484 fields

tsv-join \
	basic_clinical.tsv \
	--data-fields 1 \
	-f gene.tmp \
	--key-fields 1 \
	--append-fields 2-60484 \
	>gene.tsv

sed '1d' gene.tsv | datamash check
#504 lines, 60486 fields

cat gene.tsv \
	| head -n 5 \
	| tsv-select -f 1-5

rm gene.tmp

cd rawdata/LUAD

find methylation_beta_value -type f -name "*.gdc_hg38.txt" \
	| head -n 1 \
	| xargs cat \
	| cut -f 1 \
	| sed -e "1d" \
	| (echo "#sample" && cat) \
	| datamash transpose \
		>meth.tmp

tsv-filter -H --str-eq 2:methylation_beta_value \
	--str-in-fld 4:"Primary Tumor" \
	<fileinfo.tsv \
	| keep-header -- sort -k3,3 -k4,4r \
	| tsv-uniq -H -f 3 \
	| tsv-select -H -f 1,3 \
	| grep -v "^#" \
	| parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo {2} 1>&2
        filepath=$(find methylation_beta_value -type f -name "{1}")
        cat ${filepath} |
            tsv-select -H -f 1,2 |
            sed -e "1d" |
            (echo -e '\''#sample\t{2}'\'' && cat) |
            datamash transpose |
            grep -v "^#"
    ' \
	| perl -nla -F'\t' -e '
        for my $i ( 1 .. $#F) {
            if ($F[$i] ne q{NA}) {
                $F[$i] = sprintf qq{%.6f}, $F[$i];
            }
        }
        print join qq{\t}, @F;
    ' \
	| sort \
		>>meth.tmp

sed '1d' meth.tmp | datamash check
#458 lines, 485578 fields

tsv-join \
	basic_clinical.tsv \
	--data-fields 1 \
	-f meth.tmp \
	--key-fields 1 \
	--append-fields 2-485578 \
	>meth.tsv

sed '1d' meth.tsv | datamash check
#449 lines, 485580 fields

head -n 5 meth.tsv \
	| tsv-select -f 1-5

pigz meth.tsv
rm meth.tmp

pigz -dcf meth.tsv.gz \
	| tsv-select -f 1-3 \
	| tsv-join \
		-H \
		--filter-file detailed_clinical.tsv \
		--key-fields 1 \
		--append-fields 2-15 \
		>meth_clinical.tsv

sed '1d' meth_clinical.tsv | datamash check
#449 lines, 14 fields

```

## mirna_expression

```shell
cd ~/data/cancer/rawdata/LUAD

find mirna_expression -type f -name "*.mirnas.quantification.txt" |
    head -n 1 |
    xargs cat |
    cut -f 1 |
    sed -e "1d" |
    perl -nlp -e 's/\-/_/g' | # hyphens in miRNA names confused coxph()
    (echo "#sample" && cat) |
    datamash transpose \
    > mirna.tmp

cat fileinfo.tsv |
    grep ".mirnas.quantification.txt" |
    tsv-filter -H --str-eq 2:mirna_expression --str-in-fld 4:"Primary Tumor" |
    keep-header -- sort -k3,3 -k4,4r |
    tsv-uniq -H -f 3 |
    tsv-select -H -f 1,3 |
    grep -v "^#" |
    parallel --colsep '\t' --no-run-if-empty --linebuffer -k -j 1 '
        echo {2} 1>&2
        filepath=$(find mirna_expression -type f -name "{1}")
        cat ${filepath} |
            tsv-select -H -f 1,3 |
            sed -e "1d" |
            (echo -e '\''#sample\t{2}'\'' && cat) |
            datamash transpose |
            grep -v "^#"
    ' |
    perl -nla -F'\t' -e '
        for my $i ( 1 .. $#F) {
            if ($F[$i] ne q{NA}) {
                $F[$i] = sprintf qq{%.6f}, $F[$i];
            }
        }
        print join qq{\t}, @F;
    ' |
    sort \
    >> mirna.tmp

sed '1d' mirna.tmp | datamash check
#513 lines, 1882 fields

tsv-join \
    basic_clinical.tsv \
    --data-fields 1 \
    -f mirna.tmp \
    --key-fields 1 \
    --append-fields 2-1882 \
    > mirna.tsv

sed '1d' mirna.tsv | datamash check
#504 lines, 1884 fields

cat mirna.tsv |
    head -n 5 |
    tsv-select -f 1-5

rm mirna.tmp

```

# Extra validation datasets

## GSE39279

* <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39279>
* The CURELUNG Project from the Epigenetics Side (Patients dataset)
* LUAD and LUSC

*recurrence* as events

```shell
mkdir -p GEO/GSE39279
cd GEO/GSE39279 || exit
```

```r
library(Biobase)
library(GEOquery)
library(readr)

# load series and platform data from GEO
gset <- getGEO("GSE39279", GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

summary(gset@phenoData@data)

sink('GEO.basic_clinical.tmp')

cat(gset@phenoData@data[["geo_accession"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["title"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["adjuvant chemotherapy:ch1"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["age:ch1"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["gender:ch1"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["nsclc type:ch1"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["recurrence:ch1"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["smoker:ch1"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["Stage:ch1"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["time_rec (years):ch1"]], sep = "\t")
cat("\n")

sink()

assay_df <- data.frame(rownames(gset@assayData[["exprs"]]), gset@assayData[["exprs"]])
colnames(assay_df)[1] <- "#sample"
write_tsv(assay_df, file = "GEO.assay.tmp")

```

```shell
datamash transpose <GEO.basic_clinical.tmp \
	| perl -nla -F"\t" -e '
        $F[9] =~ /NA|NaN/i and next;
        $F[9] = int($F[9] * 365);
        if ($F[6] =~ /yes/i) {
            $F[6] =  q{1};
        } elsif ($F[6] =~ /no/i) {
            $F[6] =  q{0};
        } elsif ($F[13] =~ /NA|NaN/i) {
            next;
        }
        print join qq{\t}, $F[0], $F[9], $F[6], $F[5];
    ' \
	| (echo -e "#sample\ttime\tstatus\tnsclc type" && cat) \
	| tsv-filter -H --istr-eq 4:adenocarcinoma \
	| tsv-select -f 1-3 \
		>GEO.basic_clinical.2.tmp

sed '1d' GEO.basic_clinical.2.tmp | datamash check
#155 lines, 3 fields

sed '1d' GEO.assay.tmp | datamash check
#485577 lines, 445 fields

# 5 didits
perl -nla -F'\t' -e '
        /^#/ and print and next;
        for my $m (@F) {
            $m eq $F[0] and next;
            $m = sprintf "%.5f", $m;
        }
        print join qq{\t}, @F;
    ' <GEO.assay.tmp \
	>GEO.assay.2.tmp

tsv-join \
	GEO.basic_clinical.2.tmp \
	-d 1 \
	-f <(datamash transpose <GEO.assay.2.tmp) \
	-k 1 \
	-a 2-485578 \
	>GEO.tsv

sed '1d' GEO.tsv | datamash check
#155 lines, 485580 fields

pigz GEO.tsv
pigz -dcf GEO.tsv.gz \
	| head -n 5 \
	| tsv-select -f 1-5
rm ./*.tmp

pigz -dcf GEO.tsv.gz \
	| tsv-filter -H --str-eq 3:1 \
	| tsv-summarize --header --median 2 --stdev 2 --count
#time_median	time_stdev	count
#572.5	878.976437386	68

```

## GSE66836

```shell
mkdir -p GEO/GSE66836
cd GEO/GSE66836 || exit
aria2c -x 12 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66836/matrix/GSE66836_series_matrix.txt.gz
```

```r
library(Biobase)
library(GEOquery)
library(readr)

# load series and platform data from GEO
gset <- getGEO("GSE66836", destdir = '.', GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

summary(gset@phenoData@data)

sink('GEO.basic_clinical.tmp')

cat(gset@phenoData@data[["geo_accession"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["title"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["source_name_ch1"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["gender:ch1"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["tissue:ch1"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["Stage:ch1"]], sep = "\t")
cat("\n")

sink()

assay_df <- data.frame(rownames(gset@assayData[["exprs"]]), gset@assayData[["exprs"]])
colnames(assay_df)[1] <- "#sample"
write_tsv(assay_df, file = "GEO.assay.tmp")

```

## GSE119144

<https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-020-00907-4>

*PFS*

```shell
mkdir -p GEO/GSE119144
cd ~/data/cancer/GEO/GSE119144

curl -O https://static-content.springer.com/esm/art%3A10.1186%2Fs13148-020-00907-4/MediaObjects/13148_2020_907_MOESM6_ESM.xlsx

plotr xlsx 13148_2020_907_MOESM6_ESM.xlsx --sheet Table1 -o stdout |
    sed '1,2d' |
    sed '1s/ /_/' \
    > pfs.tsv

cat pfs.tsv |
    tsv-summarize -H --group-by PFS_event,benefit --count
#PFS_event	benefit	count
#1	N	44
#0	Y	9
#1	Y	5
#0	N	2

cat pfs.tsv |
    tsv-summarize -H --mean PFS --median PFS
#PFS_mean	PFS_median
#123.866666667	56

cat pfs.tsv |
    tsv-filter -H --le PFS:30 |
    tsv-summarize -H --count
#count
#12

curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE119nnn/GSE119144/matrix/GSE119144_series_matrix.txt.gz

R

```

```r
library(Biobase)
library(GEOquery)
library(readr)

# load series and platform data from GEO

gset <- getGEO("GSE119144", destdir = '.', GSEMatrix = TRUE, getGPL = FALSE)
if (length(gset) > 1) idx <- grep("GPL23976", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

summary(gset@phenoData@data)

sink('GEO.basic_clinical.tmp')

cat(gset@phenoData@data[["geo_accession"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["title"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["source_name_ch1"]], sep = "\t")
cat("\n")
cat(gset@phenoData@data[["disease state:ch1"]], sep = "\t")
cat("\n")

sink()

assay_df <- data.frame(rownames(gset@assayData[["exprs"]]), gset@assayData[["exprs"]])
colnames(assay_df)[1] <- "#sample"
write_tsv(assay_df, file = "GEO.assay.tmp")

```

```shell
cd ~/data/cancer/GEO/GSE119144

cat GEO.basic_clinical.tmp |
    datamash transpose |
    perl -nla -F"\t" -e '
        $F[2] =~ s/NSCLC //;
        print join qq{\t}, $F[0], $F[2];
    ' |
    (echo -e "#sample\tPatient_ID" && cat) \
    > GEO.basic_clinical.2.tmp

cat GEO.basic_clinical.2.tmp |
    tsv-join -H -f pfs.tsv --key-fields 'Patient_ID'  --append-fields 2,3 |
    tsv-select -f 1,3,4 |
    sed '1d' |
    (echo -e "#sample\ttime\tstatus" && cat) \
    > GEO.basic_clinical.pfs.tmp

# benefit as non-events
cat GEO.basic_clinical.2.tmp |
    tsv-join -H -f pfs.tsv --key-fields 'Patient_ID'  --append-fields 2,4 |
    tsv-select -f 1,3,4 |
    sed '1d' |
    perl -nla -F"\t" -e '
        if ($F[2] =~ /Y/i) {
            $F[2] =  q{0};
        } elsif ($F[2] =~ /N/i) {
            $F[2] =  q{1};
        }
        print join qq{\t}, $F[0], $F[1], $F[2];
    ' |
    (echo -e "#sample\ttime\tstatus" && cat) \
    > GEO.basic_clinical.benefit.tmp

sed '1d' GEO.basic_clinical.pfs.tmp | datamash check
#60 lines, 3 fields

sed '1d' GEO.assay.tmp | datamash check
#866091 lines, 61 fields

# 5 didits
cat GEO.assay.tmp |
    perl -nla -F'\t' -e '
        /^#/ and print and next;
        for my $m (@F) {
            $m eq $F[0] and next;
            $m = sprintf "%.5f", $m;
        }
        print join qq{\t}, @F;
    ' \
    > GEO.assay.2.tmp

tsv-join \
    GEO.basic_clinical.pfs.tmp \
    -d 1 \
    -f <(cat GEO.assay.2.tmp | datamash transpose) \
    -k 1 \
    -a 2-866092 \
    > GEO.pfs.tsv

tsv-join \
    GEO.basic_clinical.benefit.tmp \
    -d 1 \
    -f <(cat GEO.assay.2.tmp | datamash transpose) \
    -k 1 \
    -a 2-866092 \
    > GEO.benefit.tsv

sed '1d' GEO.pfs.tsv | datamash check
#60 lines, 866094 fields

sed '1d' GEO.benefit.tsv | datamash check
#60 lines, 866094 fields

pigz GEO.pfs.tsv
pigz GEO.benefit.tsv

gzip -dcf GEO.pfs.tsv.gz |
    head -n 5 |
    tsv-select -f 1-5

gzip -dcf GEO.benefit.tsv.gz |
    head -n 5 |
    tsv-select -f 1-5

rm *.tmp

pigz -dcf GEO.pfs.tsv.gz |
    tsv-filter -H --str-eq status:1 |
    tsv-summarize --header --median 2 --stdev 2 --count
#time_median	time_stdev	count
#44	126.286515565	49

```

## XH33

33 patients from XH

`Normalized_Beta.csv.gz` from EPIC chip

```shell
mkdir -p ~/data/cancer/ours/XH33
cd ~/data/cancer/ours/XH33

cat <<EOF | mlr --icsv --otsv cat > status.tsv
#sample,time,status
X002,755,1
X003,990,1
X005,1323,0
X007,526,1
X008,278,1
X009,1259,0
X011,296,1
X013,533,1
X014,136,1
X015,1206,0
X016,1186,0
X021,1127,0
X022,1127,0
X023,444,1
X024,1113,0
X025,1113,0
X026,1109,0
X028,1092,0
X030,1087,0
X031,1085,0
X032,1085,0
X033,1078,0
X034,1022,0
X035,1071,0
X036,1043,0
X037,1036,0
X038,1036,0
X040,1017,0
X041,1017,0
X042,196,1
X045,877,0
X046,868,0
X047,171,1
EOF

gzip -dcf Normalized_Beta.csv.gz |
    tr -d '"' |
    sed '1s/^Probe.ID/#sample/' | # add '#' to the first line
    perl -nla -F',' -e '
        if ($. != 1) {
            for my $i ( 1 .. $#F) {
                if ($F[$i] ne q{NA}) {
                    $F[$i] = sprintf qq{%.5f}, $F[$i];
                }
            }
        }

        print join qq{\t}, @F;
    ' |
    datamash transpose \
    > beta.tsv.tmp

cat beta.tsv.tmp | datamash check
#34 lines, 865919 fields

tsv-join \
    status.tsv \
    -d 1 \
    -f beta.tsv.tmp \
    -k 1 \
    -a 2-865919 \
    > ours.tsv

pigz ours.tsv

gzip -dcf ours.tsv.gz |
    head -n 5 |
    tsv-select -f 1-5

rm *.tmp

pigz -dcf ours.tsv.gz |
    tsv-filter -H --str-eq status:1 |
    tsv-summarize --header --median 2 --stdev 2 --count
#time_median	time_stdev	count
#370	276.589647914	10

```

# Preparations

```shell
mkdir -p ~/data/cancer/LUAD
cd ~/data/cancer/LUAD

cp ~/data/cancer/rawdata/LUAD/meth_clinical.tsv .
cp ~/data/cancer/rawdata/LUAD/drug.tsv .
cp ~/data/cancer/rawdata/LUAD/radiation.tsv .
cp ~/data/cancer/rawdata/LUAD/meth.tsv.gz .

cp ~/data/cancer/rawdata/LUAD/gencode.tsv .
cp ~/data/cancer/rawdata/LUAD/gene.tsv .

bmr kb cgc -o Cancer_Gene_Census.tsv
bmr kb epic -o epic.tsv

find ~/data/cancer/rawdata/LUAD/methylation_beta_value -type f -name "*.gdc_hg38.txt" |
    head -n 1 |
    xargs cat |
    cut -f 1,3-11 \
    > site.tsv

# lines of events
# pigz -dcf meth.tsv.gz | perl -nla -F"\t" -e '$F[2] == 1 and print $.'

# 60%
pigz -dcf meth.tsv.gz |
    keep-header -- head -n 269 |
    pigz > training.tsv.gz

# 40%
pigz -dcf meth.tsv.gz |
    keep-header -- tail -n 180 |
    pigz > testing.tsv.gz

# 365 * 5 = 1825
bmr cox \
    --all \
    --threshold 1825 \
    --parallel 20 \
    --sample 135 \
    --count 100

```

```shell
pigz -dcf training.tsv.gz | sed '1d' | datamash check
#269 lines, 485580 fields

pigz -dcf testing.tsv.gz | sed '1d' | datamash check
#180 lines, 485580 fields

pigz -dcf training.tsv.gz |
    tsv-filter -H --str-eq 3:1 |
    tsv-summarize -H --median 2 --stdev 2 --count
#time_median     time_stdev      count
#594     684.725550784   97

pigz -dcf testing.tsv.gz |
    tsv-filter -H --str-eq 3:1 |
    tsv-summarize -H --median 2 --stdev 2 --count
#time_median     time_stdev      count
#589.5   780.609661388   64

pigz -dcf meth.tsv.gz |
    grep -v "^#" |
    tsv-filter --str-eq 3:1 |
    tsv-summarize --median 2 --stdev 2 --count
#593     722.310642028   161

pigz -dcf meth.tsv.gz |
    grep -v "^#" |
    tsv-filter --str-eq 3:0 |
    tsv-summarize --median 2 --stdev 2 --count
#655     1005.22112009   288

pigz -dcf training.tsv.gz |
    tsv-filter -H --str-ne 4:NA |
    tsv-summarize -H --mean 4 --stdev 4
#cg00000029_mean cg00000029_stdev
#0.300949772388  0.0937019049872

```

## Clinical features

```shell
mkdir -p clinical

# gender
cat meth_clinical.tsv |
    tsv-summarize -H --group-by gender --count |
    keep-header -- sort
#gender  count
#FEMALE  238
#MALE    211

cat meth_clinical.tsv |
    tsv-filter -H --str-eq gender:MALE |
    cut -f 1 \
    > clinical/MALE.tsv

cat meth_clinical.tsv |
    tsv-filter -H --str-eq gender:FEMALE |
    cut -f 1 \
    > clinical/FEMALE.tsv

# stage
cat meth_clinical.tsv |
    tsv-summarize -H --group-by ajcc_pathologic_tumor_stage --count |
    keep-header -- sort
#ajcc_pathologic_tumor_stage     count
#        1
#[Discrepancy]   4
#Stage I 5
#Stage IA        119
#Stage IB        120
#Stage II        1
#Stage IIA       49
#Stage IIB       60
#Stage IIIA      61
#Stage IIIB      9
#Stage IV        20

cat meth_clinical.tsv |
    tsv-filter -H --regex 'ajcc_pathologic_tumor_stage:^Stage (I|IA|IB)$' |
    cut -f 1 \
    > clinical/Stage_I.tsv

cat meth_clinical.tsv |
    tsv-filter -H --regex 'ajcc_pathologic_tumor_stage:^Stage (II|IIA|IIB)$'|
    cut -f 1 \
    > clinical/Stage_II.tsv

cat meth_clinical.tsv |
    tsv-filter -H --regex 'ajcc_pathologic_tumor_stage:^Stage (III|IIIA|IIIB|IV)$' |
    cut -f 1 \
    > clinical/Stage_III_IV.tsv

# age at diagnosis
cat meth_clinical.tsv |
    tsv-summarize -H --exclude-missing --quantile age_at_initial_pathologic_diagnosis:0.33,0.66
#age_at_initial_pathologic_diagnosis_pct33       age_at_initial_pathologic_diagnosis_pct66
#60      70

cat meth_clinical.tsv |
    tsv-filter -H --not-blank age_at_initial_pathologic_diagnosis --le age_at_initial_pathologic_diagnosis:59 |
    cut -f 1 \
    > clinical/Age_59.tsv

cat meth_clinical.tsv |
    tsv-filter -H --not-blank age_at_initial_pathologic_diagnosis --ge age_at_initial_pathologic_diagnosis:60 --le age_at_initial_pathologic_diagnosis:69 |
    cut -f 1 \
    > clinical/Age_60_69.tsv

cat meth_clinical.tsv |
    tsv-filter -H --not-blank age_at_initial_pathologic_diagnosis --ge age_at_initial_pathologic_diagnosis:70 |
    cut -f 1 \
    > clinical/Age_70.tsv

# mutation
cat meth_clinical.tsv |
    tsv-summarize -H --group-by kras_mutation_found --count |
    keep-header -- sort
#kras_mutation_found     count
#        387
#NO      39
#YES     23

cat meth_clinical.tsv |
    tsv-summarize -H --group-by egfr_mutation_status --count |
    keep-header -- sort
#egfr_mutation_status    count
#        194
#NO      175
#YES     80

cat meth_clinical.tsv |
    tsv-filter -H --str-eq egfr_mutation_status:YES |
    cut -f 1 \
    > clinical/egfr_yes.tsv

cat meth_clinical.tsv |
    tsv-filter -H --str-eq egfr_mutation_status:NO |
    cut -f 1 \
    > clinical/egfr_no.tsv

# drug
cat drug.tsv |
    tsv-summarize -H --group-by 2 --count |
    keep-header -- sort -t $'\t' -k2nr

cat drug.tsv |
    tsv-summarize -H --group-by 3 --count |
    keep-header -- sort -t $'\t' -k2nr

cat drug.tsv |
    tsv-filter --istr-in-fld 3:Chemotherapy |
    tsv-summarize --group-by 2 --count |
    keep-header -- sort -t $'\t' -k2nr

cat drug.tsv |
    tsv-summarize -H --group-by 1 --unique-values 2 \
    > drug_values.tsv

cat drug_values.tsv |
    tsv-filter \
        -H --or \
        --istr-in-fld 2:cetuximab --istr-in-fld 2:erbitux \
        --istr-in-fld 2:erlotinib --istr-in-fld 2:tarceva --istr-in-fld 2:erlotonib \
        --istr-in-fld 2:gefitinib --istr-in-fld 2:iressa \
    > clinical/drug_egfr.tsv

cat drug_values.tsv |
    tsv-filter \
        -H --or \
        --istr-in-fld 2:cetuximab --istr-in-fld 2:erbitux \
        --istr-in-fld 2:erlotinib --istr-in-fld 2:tarceva --istr-in-fld 2:erlotonib \
        --istr-in-fld 2:gefitinib --istr-in-fld 2:iressa \
        --istr-in-fld 2:bevacizumab --istr-in-fld 2:avastin --istr-in-fld 2:mvasi \
        --istr-in-fld 2:motesanib \
    > clinical/drug_target.tsv

cat drug_values.tsv |
    tsv-filter \
        -H --or \
        --istr-in-fld 2:cisplatin --istr-in-fld 2:platinol \
        --istr-in-fld 2:carboplatin --istr-in-fld 2:paraplatin --istr-in-fld 2:novaplus \
        --istr-in-fld 2:oxaliplatin --istr-in-fld 2:eloxatin |
    tsv-join -H --exclude -k 1 -f clinical/drug_target.tsv \
    > clinical/drug_platin.tsv

cat drug_values.tsv |
    tsv-filter \
        -H --or \
        --istr-in-fld 2:paclitaxel --istr-in-fld 2:taxol --istr-in-fld 2:onxol \
        --istr-in-fld 2:docetaxel --istr-in-fld 2:taxotere --istr-in-fld 2:docefrez \
        --istr-in-fld 2:etoposide --istr-in-fld 2:etopophos --istr-in-fld 2:toposar --istr-in-fld 2:vepeside \
        --istr-in-fld 2:vinorelbine --istr-in-fld 2:navelbine \
        --istr-in-fld 2:vinblastine --istr-in-fld 2:velban \
        --istr-in-fld 2:vincristine --istr-in-fld 2:oncovin --istr-in-fld 2:vincasar |
    tsv-join -H --exclude -k 1 -f clinical/drug_target.tsv \
    > clinical/drug_mitotic.tsv

cat drug_values.tsv |
    tsv-filter \
        -H --or \
        --istr-in-fld 2:fluorouracil --istr-in-fld 2:fluoroplex --istr-in-fld 2:"5 fu" --istr-in-fld 2:"5-fu" --istr-in-fld 2:"5- fu" --istr-in-fld 2:"5fu" \
        --istr-in-fld 2:leucovorin \
        --istr-in-fld 2:"folinic acid" \
        --istr-in-fld 2:floxuridine --istr-in-fld 2:fluorodeoxyuridine \
        --istr-in-fld 2:irinotecan --istr-in-fld 2:campto \
        --istr-in-fld 2:gemcitabine --istr-in-fld 2:gemzar \
        --istr-in-fld 2:pemetrexed --istr-in-fld 2:alimta \
        --istr-in-fld 2:aminopterin \
        --istr-in-fld 2:aminocamptothecin \
        --istr-in-fld 2:capecitabine --istr-in-fld 2:xeloda |
    tsv-join -H --exclude -k 1 -f clinical/drug_target.tsv \
    > clinical/drug_metabolites.tsv

# training and testing
pigz -dcf training.tsv.gz |
    cut -f 1 \
    > clinical/training.tsv

pigz -dcf testing.tsv.gz |
    cut -f 1 \
    > clinical/testing.tsv

for g in \
    MALE FEMALE \
    Stage_I Stage_II Stage_III_IV \
    Age_59 Age_60_69 Age_70 \
    egfr_yes egfr_no \
    drug_target drug_platin \
    ; do
    echo "${g}"

    sed '1d' clinical/training.tsv |
        wc -l

    tsv-join clinical/training.tsv \
        -H -k 1 \
        -f clinical/${g}.tsv |
        sed '1d' |
        wc -l

    sed '1d' clinical/testing.tsv |
        wc -l

    tsv-join clinical/testing.tsv \
        -H -k 1 \
        -f clinical/${g}.tsv |
        sed '1d' |
        wc -l
done |
    paste - - - - - \
    > clinical/cross.tsv

cat clinical/cross.tsv |
    parallel --colsep '\t' -j 1 -k '
        echo "==> {1}"
        Rscript -e "
            x <- matrix(c({2}, {4}, {3}, {5}), nrow=2)
            x
            chisq.test(x)
        "
    '

```

# To and From HPC

Use an LSF cluster for training

```shell
# rsync to hpcc
rsync -avP ~/data/cancer/rawdata/LUAD/ wangq@202.119.37.251:data/cancer/rawdata/LUAD
rsync -avP ~/data/cancer/LUAD/ wangq@202.119.37.251:data/cancer/LUAD

# rsync back to local
rsync -avP wangq@202.119.37.251:data/cancer/rawdata/LUAD/ ~/data/cancer/rawdata/LUAD

rsync -avP \
    --exclude "formula.tsv" \
    --exclude "cox.result.tsv" \
    wangq@202.119.37.251:data/cancer/LUAD/ \
    ~/data/cancer/LUAD

```

# 1_training

```shell
mkdir -p 1_training

bmr split training.tsv.gz -c 1000 --mode column --rc 3 -o split | bash

mkdir -p job
find split -type f -name "*[0-9]" |
    sort |
    split -l 120 -a 3 -d - job/

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo ${f}
    bsub -q mpi -n 24 -J "training-${f}" \
        "
        cat ${f} |
            parallel --no-run-if-empty --line-buffer -k -j 24 '
                echo '\''==> Processing {}'\''
                Rscript univariate.R --valid 0.99 -i {} -o {}.result.tsv
            '
        "
done

#bkill -J "training-*"

#find split -type f -name "*[0-9]" |
#    sort |
#    parallel --no-run-if-empty --line-buffer -k -j 20 '
#        echo "==> Processing {}"
#        Rscript univariate.R --valid 0.99 -i {} -o {}.result.tsv
#    '

tsv-append -H split/*.result.tsv > 1_training/cox.result.tsv

rm -fr split job BS output.*

# 0 for non-checking
bmr cox \
    --coxp 0.05 \
    --hr 0.8,1.25 \
    --kmp 0.5 \
    --rocauc 0.46,0.54 \
    --rocp 0.8

# for compatible with EPIC
cat 1_training/cox.result.tsv |
    keep-header -- bash filter_result.sh |
    keep-header -- tsv-join -f epic.tsv --key-fields 1 \
    > 1_training/result.tsv

bash result_stat.sh 1_training
bash result_dist.sh 1_training/cox.result.tsv

bash select_col.sh -f 1-3 training.tsv.gz 1_training/result.tsv > 1_training/data.tsv

```

| #Item                     | Value   |
|:--------------------------|:--------|
| 1_training/cox.result.tsv | 396066  |
| 1_training/result.tsv     | 17862   |
| count                     | 17861   |
| hr_max                    | 2.78937 |
| hr_median                 | 1.48049 |
| kmp_min                   | 0.00000 |
| kmp_median                | 0.0604  |
| rocauc_max                | 0.75559 |
| rocauc_median             | 0.59949 |
| rocp_min                  | 0.00000 |
| rocp_median               | 0.1052  |

| #Item    | Value  |
|:---------|:-------|
| #Cox P   | Count  |
| 0.05     | 32094  |
| 0.01     | 10526  |
| #HR      | Count  |
| 2        | 415    |
| 1        | 298934 |
| 0        | 96716  |
| #KM P    | Count  |
| 0.5      | 200816 |
| 0.05     | 20243  |
| 0.01     | 3960   |
| #ROC AUC | Count  |
| 0.7      | 398    |
| 0.6      | 49564  |
| 0.5      | 229008 |
| 0.4      | 111730 |
| 0.3      | 5356   |
| 0.2      | 9      |
| #ROC P   | Count  |
| 0.5      | 214897 |
| 0.05     | 28064  |
| 0.01     | 6909   |

# 1_bootstrap

```shell
mkdir -p 1_bootstrap

bash BS.sh 1_training/data.tsv

bash bootstrap.sh BS 1_training/result.tsv 1_bootstrap

mv 1_bootstrap/*.bootstrap.tsv 1_bootstrap/result.tsv
rm 1_bootstrap/*.count.tsv

rm -fr split job BS output.*

bash result_stat.sh -b 1_bootstrap

```

| #Item                  | Value   |
|:-----------------------|:--------|
| 1_bootstrap/result.tsv | 17862   |
| count                  | 17861   |
| hr_max                 | 2.78937 |
| hr_median              | 1.48049 |
| kmp_min                | 0.00000 |
| kmp_median             | 0.0604  |
| rocauc_max             | 0.75559 |
| rocauc_median          | 0.59949 |
| rocp_min               | 0.00000 |
| rocp_median            | 0.1052  |
| bin(BS)                | count   |
| 100                    | 40      |
| 80                     | 4984    |
| 60                     | 8103    |
| 40                     | 4523    |
| 20                     | 211     |

# 1_testing

```shell
BS_PASS=70

mkdir -p 1_testing

bash select_col.sh -f 1-3 testing.tsv.gz 1_training/result.tsv > 1_testing/data.tsv

cat 1_training/result.tsv |
    sed '1d' |
    parallel --pipe --block 100k -j 20 'Rscript validate.R -i 1_testing/data.tsv -f stdin --no' |
    keep-header -- grep -v '^#' \
    > 1_testing/validate.result.tsv

# 0 for non-checking
bmr cox \
    --coxp 0.05 \
    --hr 0.83,1.2 \
    --kmp 0.5 \
    --rocauc 0.48,0.52 \
    --rocp 0.8

cat 1_testing/validate.result.tsv |
    keep-header -- bash filter_result.sh \
    > 1_testing/result.tsv

cat 1_bootstrap/result.tsv |
    tsv-filter -H --ge 9:${BS_PASS} |
    tsv-join \
        -H -f 1_testing/result.tsv \
        -k 1 -a 5,6,7,8 --prefix testing- |
    keep-header -- perl -nla -F"\t" -e '
        if ($F[4] >= 1) {
            print if $F[9] >= 1;
        } else {
            print if $F[9] < 1;
        }
    ' |
    keep-header -- perl -nla -F"\t" -e '
        if ($F[6] >= 0.5) {
            print if $F[11] >= 0.5;
        } else {
            print if $F[11] < 0.5;
        }
    ' \
    > 1_testing/join.result.tsv

bash result_stat.sh 1_testing
bash result_dist.sh 1_testing/validate.result.tsv

bash select_col.sh -f 1-3 training.tsv.gz 1_testing/join.result.tsv > 1_testing/data.training.tsv
bash select_col.sh -f 1-3 testing.tsv.gz 1_testing/join.result.tsv > 1_testing/data.testing.tsv

```

| #Item                         | Value    |
|:------------------------------|:---------|
| 1_testing/data.tsv            | 181      |
| 1_testing/join.result.tsv     | 2540     |
| 1_testing/result.tsv          | 8159     |
| 1_testing/validate.result.tsv | 17862    |
| count                         | 8158     |
| hr_max                        | 3.59807  |
| hr_median                     | 1.319995 |
| kmp_min                       | 0.00000  |
| kmp_median                    | 0.19704  |
| rocauc_max                    | 0.78171  |
| rocauc_median                 | 0.54234  |
| rocp_min                      | 0.00000  |
| rocp_median                   | 0.342665 |

| #Item    | Value |
|:---------|:------|
| #Cox P   | Count |
| 0.05     | 17861 |
| 0.01     | 6449  |
| #HR      | Count |
| 3        | 2     |
| 2        | 311   |
| 1        | 11045 |
| 0        | 6503  |
| #KM P    | Count |
| 0.5      | 10092 |
| 0.05     | 1537  |
| 0.01     | 417   |
| #ROC AUC | Count |
| 0.7      | 64    |
| 0.6      | 1960  |
| 0.5      | 7525  |
| 0.4      | 6955  |
| 0.3      | 1334  |
| 0.2      | 23    |
| #ROC P   | Count |
| 0.5      | 9045  |
| 0.05     | 857   |
| 0.01     | 202   |

# 2_training

```shell
mkdir -p 2_training

bmr replace 1_testing/data.training.tsv -o 2_training/data.tsv

bmr nextstep 2_training/data.replace.tsv |
    tsv-uniq > 2_training/formula.tsv

bmr split 2_training/formula.tsv -c 5000 --mode row --rr 1 -o split | bash

mkdir -p job
find split -type f -name "*[0-9]" |
    sort |
    split -l 120 -a 3 -d - job/

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo ${f}
    bsub -q mpi -n 24 -J "training-${f}" \
        "
        cat ${f} |
            parallel --no-run-if-empty --line-buffer -k -j 24 '
                echo '\''==> Processing {}'\''
                Rscript multivariate.R -i 2_training/data.tsv -f {} -o {}.result.tsv
            '
        "
done

#bkill -J "training-*"

#find split -type f -name "*[0-9]" |
#    sort |
#    parallel --no-run-if-empty --line-buffer -k -j 20 '
#        echo "==> Processing {}"
#        Rscript multivariate.R -i 2_training/data.tsv -f {} -o {}.result.tsv
#    '

tsv-append -H split/*.result.tsv > 2_training/cox.result.tsv

rm -fr split job BS output.*

bmr cox \
    --coxp 0.05 \
    --hr 0.5,2 \
    --kmp 0.05 \
    --rocauc 0.42,0.58 \
    --rocp 0.05

cat 2_training/cox.result.tsv |
    keep-header -- parallel --pipe 'bash filter_result.sh' \
    > 2_training/result.tsv

bash result_stat.sh 2_training
bash result_dist.sh 2_training/cox.result.tsv

```

| #Item                       | Value   |
|:----------------------------|:--------|
| 2_training/cox.result.tsv   | 3221992 |
| 2_training/data.replace.tsv | 2539    |
| 2_training/data.tsv         | 270     |
| 2_training/formula.tsv      | 3221992 |
| 2_training/result.tsv       | 508082  |
| count                       | 2539    |
| hr_max                      | 4.18277 |
| hr_median                   | 2.17381 |
| kmp_min                     | 0.00000 |
| kmp_median                  | 0.00023 |
| rocauc_max                  | 0.83619 |
| rocauc_median               | 0.67168 |
| rocp_min                    | 0.00000 |
| rocp_median                 | 0.00283 |

| #Item    | Value   |
|:---------|:--------|
| #Cox P   | Count   |
| 0.05     | 1555491 |
| 0.01     | 415267  |
| #HR      | Count   |
| 4        | 2       |
| 3        | 2340    |
| 2        | 876389  |
| 1        | 2343246 |
| 0        | 14      |
| #KM P    | Count   |
| 0.5      | 3220953 |
| 0.05     | 2948363 |
| 0.01     | 2152182 |
| #ROC AUC | Count   |
| 0.8      | 124     |
| 0.7      | 261105  |
| 0.6      | 2796358 |
| 0.5      | 164398  |
| 0.4      | 6       |
| #ROC P   | Count   |
| 0.5      | 3221757 |
| 0.05     | 2705661 |
| 0.01     | 1544422 |

# 2_bootstrap

```shell
mkdir -p 2_bootstrap

bash BS.sh 2_training/data.tsv

bmr split 2_training/result.tsv -c 10000 --mode row --rr 1 -o split | bash

for f in $(find split -maxdepth 1 -type f -name "result*[0-9]" | sort); do
    echo ${f}
    bsub -q mpi -n 24 -J "bs-${f}" \
        bash bootstrap.sh BS ${f} split
done

# bkill -J "bs-*"

tsv-append -H split/*.bootstrap.tsv > 2_bootstrap/result.tsv

rm -fr split job BS output.*

bash result_stat.sh -b 2_bootstrap

```

| #Item                  | Value   |
|:-----------------------|:--------|
| 2_bootstrap/result.tsv | 508082  |
| count                  | 2539    |
| hr_max                 | 4.18277 |
| hr_median              | 2.17381 |
| kmp_min                | 0.00000 |
| kmp_median             | 0.00023 |
| rocauc_max             | 0.83619 |
| rocauc_median          | 0.67168 |
| rocp_min               | 0.00000 |
| rocp_median            | 0.00283 |
| bin(BS)                | count   |
| 100                    | 6       |
| 80                     | 13478   |
| 60                     | 93863   |
| 40                     | 238960  |
| 20                     | 158650  |
| 0                      | 3124    |

# 2_testing

```shell
BS_PASS=60

mkdir -p 2_testing

bmr replace 1_testing/data.testing.tsv -o 2_testing/data.tsv

#bmr split 2_training/result.tsv -c 5000 --mode row --rr 1 -o split | bash
#
#mkdir -p job
#find split -type f -name "*[0-9]" |
#    sort |
#    split -l 120 -a 3 -d - job/
#
#for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
#    echo ${f}
#    bsub -q mpi -n 24 -J "validate-${f}" \
#        "
#        cat ${f} |
#            parallel --no-run-if-empty --line-buffer -k -j 24 '
#                echo '\''==> Processing {}'\''
#                Rscript validate.R -i 2_testing/data.tsv -f {} -o {}.result.tsv
#            '
#        "
#done
#
#tsv-append -H split/*.result.tsv > 3_testing/validate.result.tsv
#
#rm -fr split job BS output.*

cat 2_training/result.tsv |
    sed '1d' |
    parallel --pipe --block 200k -j 20 'Rscript validate.R -i 2_testing/data.tsv -f stdin --no' |
    keep-header -- grep -v '^#' \
    > 2_testing/validate.result.tsv

# 0 for non-checking
bmr cox \
    --coxp 0.05 \
    --hr 0.55,1.80 \
    --kmp 0.05 \
    --rocauc 0.44,0.56 \
    --rocp 0.5

cat 2_testing/validate.result.tsv |
    keep-header -- parallel --pipe 'bash filter_result.sh' \
    > 2_testing/result.tsv

tsv-join \
    2_training/result.tsv \
    -d 1 \
    -f 2_testing/result.tsv \
    -k 1 \
    > 2_testing/training.result.tsv

cat 2_bootstrap/result.tsv |
    tsv-filter -H --ge 9:${BS_PASS} |
    tsv-join \
        -H -f 2_testing/result.tsv \
        -k 1 -a 5,6,7,8 --prefix testing- |
    keep-header -- perl -nla -F"\t" -e '
        if ($F[4] >= 1) {
            print if $F[9] >= 1;
        } else {
            print if $F[9] < 1;
        }
    ' |
    keep-header -- perl -nla -F"\t" -e '
        if ($F[6] >= 0.5) {
            print if $F[11] >= 0.5;
        } else {
            print if $F[11] < 0.5;
        }
    ' \
    > 2_testing/join.result.tsv

bash result_stat.sh 2_testing
bash result_dist.sh 2_testing/validate.result.tsv

```

Check `2_training/data.replace.tsv` and `2_testing/data.replace.tsv`

```shell
diff 2_training/data.replace.tsv 2_testing/data.replace.tsv

```

| #Item                         | Value   |
|:------------------------------|:--------|
| 2_testing/data.replace.tsv    | 2539    |
| 2_testing/data.tsv            | 181     |
| 2_testing/join.result.tsv     | 37047   |
| 2_testing/result.tsv          | 172434  |
| 2_testing/training.result.tsv | 172434  |
| 2_testing/validate.result.tsv | 508082  |
| count                         | 2539    |
| hr_max                        | 4.57812 |
| hr_median                     | 2.0321  |
| kmp_min                       | 0.00000 |
| kmp_median                    | 0.00602 |
| rocauc_max                    | 0.84587 |
| rocauc_median                 | 0.63716 |
| rocp_min                      | 0.00000 |
| rocp_median                   | 0.06622 |

| #Item    | Value  |
|:---------|:-------|
| #Cox P   | Count  |
| 0.05     | 508081 |
| 0.01     | 188539 |
| #HR      | Count  |
| 4        | 14     |
| 3        | 1748   |
| 2        | 95988  |
| 1        | 406565 |
| 0        | 3766   |
| #KM P    | Count  |
| 0.5      | 479241 |
| 0.05     | 251694 |
| 0.01     | 116390 |
| #ROC AUC | Count  |
| 0.8      | 177    |
| 0.7      | 23411  |
| 0.6      | 292059 |
| 0.5      | 191539 |
| 0.4      | 895    |
| #ROC P   | Count  |
| 0.5      | 471613 |
| 0.05     | 130601 |
| 0.01     | 42072  |

# 3_training

```shell
mkdir -p 3_training

bmr nextstep 2_testing/join.result.tsv 2_training/data.replace.tsv |
    tsv-uniq \
    > 3_training/formula.tsv

bmr split 3_training/formula.tsv -c 10000 --mode row --rr 1 -o split | bash

mkdir -p job
find split -type f -name "*[0-9]" |
    sort |
    split -l 120 -a 3 -d - job/

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo ${f}
    bsub -q mpi -n 24 -J "cox-${f}" \
        "
        cat ${f} |
            parallel --no-run-if-empty --line-buffer -k -j 24 '
                echo '\''==> Processing {}'\''
                Rscript multivariate.R -i 2_training/data.tsv -f {} -o {}.result.tsv
            '
        "
done

tsv-append -H split/*.result.tsv > 3_training/cox.result.tsv

rm -fr split job BS output.*

bmr cox \
    --coxp 0.01 \
    --hr 0.4,2.5 \
    --kmp 0.05 \
    --rocauc 0.30,0.70 \
    --rocp 0.05

cat 3_training/cox.result.tsv |
    keep-header -- parallel --pipe 'bash filter_result.sh' \
    > 3_training/result.tsv

bash result_stat.sh 3_training
bash result_dist.sh 3_training/cox.result.tsv

```

| #Item                     | Value    |
|:--------------------------|:---------|
| 3_training/cox.result.tsv | 89789132 |
| 3_training/formula.tsv    | 89789132 |
| 3_training/result.tsv     | 1835910  |
| count                     | 2536     |
| hr_max                    | 5.85728  |
| hr_median                 | 2.76504  |
| kmp_min                   | 0.00000  |
| kmp_median                | 0        |
| rocauc_max                | 0.88639  |
| rocauc_median             | 0.74004  |
| rocp_min                  | 0.00000  |
| rocp_median               | 0        |

| #Item    | Value     |
|:---------|:----------|
| #Cox P   | Count     |
| 0.05     | 63132457  |
| 0.01     | 9052990   |
| #HR      | Count     |
| 6        | 1         |
| 5        | 127       |
| 4        | 22132     |
| 3        | 6510852   |
| 2        | 205643384 |
| 1        | 32902784  |
| #KM P    | Count     |
| 0.5      | 245079257 |
| 0.05     | 244931894 |
| 0.01     | 242207554 |
| #ROC AUC | Count     |
| 0.8      | 294479    |
| 0.7      | 105345276 |
| 0.6      | 139414804 |
| 0.5      | 24721     |
| #ROC P   | Count     |
| 0.5      | 245079280 |
| 0.05     | 244594073 |
| 0.01     | 227360956 |

# 3_bootstrap

```shell
mkdir -p 3_bootstrap

bash BS.sh 2_training/data.tsv

bmr split 3_training/result.tsv -c 20000 --mode row --rr 1 -o split | bash

for f in $(find split -maxdepth 1 -type f -name "result*[0-9]" | sort); do
    echo ${f}
    bsub -q mpi -n 24 -J "bs-${f}" \
        bash bootstrap.sh BS ${f} split
done

# bkill -J "bs-*"

tsv-append -H split/*.bootstrap.tsv > 3_bootstrap/result.tsv

rm -fr split job BS output.*

bash result_stat.sh -b 3_bootstrap

```

| #Item                  | Value   |
|:-----------------------|:--------|
| 3_bootstrap/result.tsv | 1835910 |
| count                  | 2536    |
| hr_max                 | 5.85728 |
| hr_median              | 2.76504 |
| kmp_min                | 0.00000 |
| kmp_median             | 0       |
| rocauc_max             | 0.88639 |
| rocauc_median          | 0.74004 |
| rocp_min               | 0.00000 |
| rocp_median            | 0       |
| bin(BS)                | count   |
| 100                    | 225     |
| 80                     | 140869  |
| 60                     | 654224  |
| 40                     | 941574  |
| 20                     | 99016   |
| 0                      | 1       |

# 3_testing

```shell
BS_PASS=60

mkdir -p 3_testing

bmr split 3_training/result.tsv -c 5000 --mode row --rr 1 -o split | bash

mkdir -p job
find split -type f -name "*[0-9]" |
    sort |
    split -l 120 -a 3 -d - job/

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo ${f}
    bsub -q mpi -n 24 -J "validate-${f}" \
        "
        cat ${f} |
            parallel --no-run-if-empty --line-buffer -k -j 24 '
                echo '\''==> Processing {}'\''
                Rscript validate.R -i 2_testing/data.tsv -f {} -o {}.result.tsv
            '
        "
done

tsv-append -H split/*.result.tsv > 3_testing/validate.result.tsv

rm -fr split job BS output.*

# 0 for non-checking
bmr cox \
    --coxp 0.01 \
    --hr 0.5,2 \
    --kmp 0.05 \
    --rocauc 0.35,0.65 \
    --rocp 0.05

cat 3_testing/validate.result.tsv |
    keep-header -- parallel --pipe 'bash filter_result.sh' \
    > 3_testing/result.tsv

tsv-join \
    3_training/result.tsv \
    -d 1 \
    -f 3_testing/result.tsv \
    -k 1 \
    > 3_testing/training.result.tsv

cat 3_bootstrap/result.tsv |
    tsv-filter -H --ge 9:${BS_PASS} |
    tsv-join \
        -H -f 3_testing/result.tsv \
        -k 1 -a 5,6,7,8 --prefix testing- |
    keep-header -- perl -nla -F"\t" -e '
        if ($F[4] >= 1) {
            print if $F[9] >= 1;
        } else {
            print if $F[9] < 1;
        }
    ' |
    keep-header -- perl -nla -F"\t" -e '
        if ($F[6] >= 0.5) {
            print if $F[11] >= 0.5;
        } else {
            print if $F[11] < 0.5;
        }
    ' \
    > 3_testing/join.result.tsv

bash result_stat.sh 3_testing
bash result_dist.sh 3_testing/validate.result.tsv

```

| #Item                         | Value   |
|:------------------------------|:--------|
| 3_testing/join.result.tsv     | 275218  |
| 3_testing/result.tsv          | 635576  |
| 3_testing/training.result.tsv | 635576  |
| 3_testing/validate.result.tsv | 1835910 |
| count                         | 2512    |
| hr_max                        | 5.93357 |
| hr_median                     | 2.35423 |
| kmp_min                       | 0.00000 |
| kmp_median                    | 0.00095 |
| rocauc_max                    | 0.87255 |
| rocauc_median                 | 0.68628 |
| rocp_min                      | 0.00000 |
| rocp_median                   | 0.00902 |

| #Item    | Value   |
|:---------|:--------|
| #Cox P   | Count   |
| 0.05     | 1835909 |
| 0.01     | 1835909 |
| #HR      | Count   |
| 5        | 26      |
| 4        | 1195    |
| 3        | 53998   |
| 2        | 959259  |
| 1        | 821301  |
| 0        | 130     |
| #KM P    | Count   |
| 0.5      | 1832417 |
| 0.05     | 1590222 |
| 0.01     | 1105768 |
| #ROC AUC | Count   |
| 0.8      | 1770    |
| 0.7      | 278387  |
| 0.6      | 1358119 |
| 0.5      | 197568  |
| 0.4      | 65      |
| #ROC P   | Count   |
| 0.5      | 1826150 |
| 0.05     | 1026750 |
| 0.01     | 449392  |

# 3_non_rare

```shell
RARE_PASS=5

mkdir -p 3_non_rare

# Count markers
cat 3_testing/join.result.tsv |
    cut -f 1 |
    perl -nl -e 'print for split /\+/' |
    tsv-uniq --number |
    tsv-summarize --group-by 1 --max 2 |
    sort -k2,2nr \
    > 3_non_rare/count.tsv

# remove rares
head -n 1 3_testing/join.result.tsv > 3_non_rare/result.tsv
cat 3_testing/join.result.tsv |
    grep -F -f <(tsv-filter --le 2:${RARE_PASS} 3_non_rare/count.tsv | cut -f 1 ) -w -v \
    >> 3_non_rare/result.tsv

cat 3_non_rare/result.tsv |
    cut -f 1 |
    perl -nl -e 'print for split /\+/' |
    tsv-uniq --number |
    tsv-summarize --group-by 1 --max 2 |
    sort -k2,2nr \
    > 3_non_rare/non-rare.count.tsv

bash result_stat.sh 3_non_rare

```

| #Item                         | Value   |
|:------------------------------|:--------|
| 3_non_rare/count.tsv          | 2484    |
| 3_non_rare/non-rare.count.tsv | 2263    |
| 3_non_rare/result.tsv         | 274609  |
| count                         | 2262    |
| hr_max                        | 5.85728 |
| hr_median                     | 2.97758 |
| kmp_min                       | 0.00000 |
| kmp_median                    | 0       |
| rocauc_max                    | 0.87934 |
| rocauc_median                 | 0.754   |
| rocp_min                      | 0.00000 |
| rocp_median                   | 0       |

# Clinical groups

```shell
for g in \
    MALE FEMALE \
    Stage_I Stage_II Stage_III_IV \
    Age_59 Age_60_69 Age_70 \
    egfr_yes egfr_no \
    drug_egfr drug_target drug_platin drug_mitotic drug_metabolites \
    ; do
    echo "==> ${g}"
    mkdir -p clinical/${g}
    tsv-append -H 2_training/data.tsv 2_testing/data.tsv |
        tsv-join \
            -H -k 1 \
            -f clinical/${g}.tsv \
        > clinical/${g}/data.tsv
    datamash check < clinical/${g}/data.tsv
done

# chemo
mkdir -p clinical/chemo_yes
tsv-append -H 2_training/data.tsv 2_testing/data.tsv |
    tsv-join \
        -H -k 1 \
        -f drug_values.tsv \
    > clinical/chemo_yes/data.tsv

mkdir -p clinical/chemo_no
tsv-append -H 2_training/data.tsv 2_testing/data.tsv |
    tsv-join \
        -H -k 1 \
        -f drug_values.tsv \
        --exclude \
    > clinical/chemo_no/data.tsv

# radiation
cat radiation.tsv |
    tsv-summarize -H --group-by 1 --unique-values 3 \
    > radiation_values.tsv

mkdir -p clinical/radiation_yes
tsv-append -H 2_training/data.tsv 2_testing/data.tsv |
    tsv-join \
        -H -k 1 \
        -f radiation_values.tsv \
    > clinical/radiation_yes/data.tsv

mkdir -p clinical/radiation_no
tsv-append -H 2_training/data.tsv 2_testing/data.tsv |
    tsv-join \
        -H -k 1 \
        -f radiation_values.tsv \
        --exclude \
    > clinical/radiation_no/data.tsv

```

# 3_clinical

```shell
mkdir -p 3_clinical

bmr cox \
    --coxp 0.05 \
    --hr 0.01,1.5 \
    --kmp 0.5 \
    --rocauc 0.01,0.60 \
    --rocp 0.5

for c in \
    MALE FEMALE \
    Stage_I Stage_II Stage_III_IV \
    Age_59 Age_60_69 Age_70 \
    egfr_yes egfr_no \
    chemo_yes chemo_no \
    radiation_yes radiation_no \
    drug_target drug_platin drug_mitotic drug_metabolites \
    ; do
    echo "==> ${c}"

    if [ -f 3_clinical/${c}.result.tsv ]; then
        continue
    fi

    bsub -q mpi -n 24 -J "clinical-${c}" \
        "
        cat 3_non_rare/result.tsv |
            sed '1d' |
            parallel --pipe --block 200k -j 24 'Rscript validate.R -i clinical/${c}/data.tsv -f stdin --no' |
            keep-header -- grep -v '^#' |
            keep-header -- parallel --pipe 'bash filter_result.sh' \
            > 3_clinical/${c}.result.tsv
        "
done

# bkill -J "clinical-*"

cp 3_non_rare/result.tsv 3_clinical/tmp.pass.result.tsv
for c in \
    MALE FEMALE \
    Stage_I Stage_II Stage_III_IV \
    Age_59 Age_60_69 Age_70 \
    egfr_yes egfr_no \
    chemo_yes chemo_no \
    radiation_yes radiation_no \
    drug_target drug_platin drug_mitotic drug_metabolites \
    ; do
    echo "==> ${c}"
    wc -l 3_clinical/tmp.pass.result.tsv

    if [ ! -f 3_clinical/${c}.result.tsv ]; then
        continue
    fi

    tsv-join \
        3_clinical/tmp.pass.result.tsv \
        -d 1 \
        -f 3_clinical/${c}.result.tsv \
        -k 1 \
        > 3_clinical/tmp.result.tsv

    mv 3_clinical/tmp.result.tsv 3_clinical/tmp.pass.result.tsv
    wc -l 3_clinical/tmp.pass.result.tsv
done | tee 3_clinical/clinical.log

mv 3_clinical/tmp.pass.result.tsv 3_clinical/result.tsv
rm 3_clinical/tmp.*.tsv

bash result_stat.sh 3_clinical

```

| #Item                                  | Value   |
|:---------------------------------------|:--------|
| 3_clinical/Age_59.result.tsv           | 263574  |
| 3_clinical/Age_60_69.result.tsv        | 270746  |
| 3_clinical/Age_70.result.tsv           | 261952  |
| 3_clinical/chemo_no.result.tsv         | 274609  |
| 3_clinical/chemo_yes.result.tsv        | 249488  |
| 3_clinical/drug_metabolites.result.tsv | 151316  |
| 3_clinical/drug_mitotic.result.tsv     | 193192  |
| 3_clinical/drug_platin.result.tsv      | 163969  |
| 3_clinical/egfr_no.result.tsv          | 273871  |
| 3_clinical/egfr_yes.result.tsv         | 235675  |
| 3_clinical/FEMALE.result.tsv           | 270823  |
| 3_clinical/MALE.result.tsv             | 274515  |
| 3_clinical/radiation_no.result.tsv     | 274609  |
| 3_clinical/radiation_yes.result.tsv    | 269635  |
| 3_clinical/result.tsv                  | 58150   |
| 3_clinical/Stage_III_IV.result.tsv     | 247286  |
| 3_clinical/Stage_II.result.tsv         | 216794  |
| 3_clinical/Stage_I.result.tsv          | 274176  |
| count                                  | 2177    |
| hr_max                                 | 5.65191 |
| hr_median                              | 3.00791 |
| kmp_min                                | 0.00000 |
| kmp_median                             | 0       |
| rocauc_max                             | 0.87724 |
| rocauc_median                          | 0.76136 |
| rocp_min                               | 0.00000 |
| rocp_median                            | 0       |

# Extra validations

## GSE39279 LUAD

```shell
mkdir -p validating/GSE39279

cat 2_training/data.replace.tsv |
    tsv-select -f 2,1 \
    > validating/reverse.replace.tsv

bash select_col.sh -f 1-3 \
    ~/data/cancer/GEO/GSE39279/GEO.tsv.gz \
    validating/reverse.replace.tsv \
    > validating/tmp.data.tsv

bmr replace validating/tmp.data.tsv validating/reverse.replace.tsv \
    > validating/GSE39279/data.tsv

cat validating/GSE39279/data.tsv |
    datamash check
#156 lines, 2753 fields

rm validating/tmp.data.tsv

```

## GSE119144 PFS

```shell
mkdir -p validating/GSE119144_PFS

bash select_col.sh -f 1-3 \
    ~/data/cancer/GEO/GSE119144/GEO.pfs.tsv.gz \
    validating/reverse.replace.tsv \
    > validating/tmp.data.tsv

bmr replace validating/tmp.data.tsv validating/reverse.replace.tsv \
    > validating/GSE119144_PFS/data.tsv

cat validating/GSE119144_PFS/data.tsv |
    datamash check
#61 lines, 2752 fields

rm validating/tmp.data.tsv

```

## GSE119144 benefit

```shell
mkdir -p validating/GSE119144_benefit

bash select_col.sh -f 1-3 \
    ~/data/cancer/GEO/GSE119144/GEO.benefit.tsv.gz \
    validating/reverse.replace.tsv \
    > validating/tmp.data.tsv

bmr replace validating/tmp.data.tsv validating/reverse.replace.tsv \
    > validating/GSE119144_benefit/data.tsv

cat validating/GSE119144_benefit/data.tsv |
    datamash check
#61 lines, 2752 fields

rm validating/tmp.data.tsv

```

## XH33

```shell
mkdir -p validating/XH33

bash select_col.sh -f 1-3 \
    ~/data/cancer/ours/XH33/ours.tsv.gz \
    validating/reverse.replace.tsv \
    > validating/tmp.data.tsv

bmr replace validating/tmp.data.tsv validating/reverse.replace.tsv \
    > validating/XH33/data.tsv

cat validating/XH33/data.tsv |
    datamash check
#34 lines, 2753 fields

rm validating/tmp.data.tsv

```

# 3_validating

GSE39279 and GSE119144 are PFS.

XH33 is OS. Data collected from 2016 to 2018.

GSE119144 has a short follow-up period.

> If the patient experienced a partial response (PR) or stable disease (SD) of more than 6 months,
> he or she was categorized as responder receiving durable clinical benefit (DCB). Conversely, if a
> patient experienced progressive disease (PD) or SD of less than 6 months, he or she was
> categorized as non-responder with non-durable benefit (NDB).

So set -t to 180

```shell
mkdir -p 3_validating

bmr cox \
    --coxp 0.05 \
    --hr 0.67,1.5 \
    --kmp 0.5 \
    --rocauc 0.01,0.55 \
    --rocp 0

# some markers don't present
for c in \
    GSE39279 GSE119144_PFS XH33 \
    ; do
    head -n 1 validating/${c}/data.tsv |
        datamash transpose
done |
    tsv-uniq --at-least 3 |
    grep '^m' \
    > 3_validating/valid.txt

cat validating/reverse.replace.tsv |
    tsv-select -f 2 |
    grep -v -F -f 3_validating/valid.txt \
    > 3_validating/invalid.txt

# continues with 3_clinical
cat 3_clinical/result.tsv |
    grep -v -F -f 3_validating/invalid.txt \
    > 3_validating/tmp.pass.result.tsv

#
rm 3_validating/validation.log

for c in \
    GSE39279 \
    ; do
    echo "==> ${c}"
    wc -l 3_validating/tmp.pass.result.tsv

    cat 3_validating/tmp.pass.result.tsv |
        sed '1d' |
        parallel --pipe --block 200k -j 20 "Rscript validate.R -i validating/${c}/data.tsv -f stdin --no" |
        keep-header -- grep -v '^#' |
        keep-header -- parallel --pipe 'bash filter_result.sh' \
        > 3_validating/tmp.filter.result.tsv

    tsv-join \
        3_validating/tmp.pass.result.tsv \
        -d 1 \
        -f 3_validating/tmp.filter.result.tsv \
        -k 1 \
        > 3_validating/tmp.result.tsv

    mv 3_validating/tmp.result.tsv 3_validating/tmp.pass.result.tsv
    wc -l 3_validating/tmp.pass.result.tsv
done | tee -a 3_validating/validation.log

for c in \
    XH33 \
    ; do
    echo "==> ${c}"
    wc -l 3_validating/tmp.pass.result.tsv

    cat 3_validating/tmp.pass.result.tsv |
        sed '1d' |
        parallel --pipe --block 200k -j 20 "Rscript validate.R -t 1095 -i validating/${c}/data.tsv -f stdin --no" |
        keep-header -- grep -v '^#' |
        keep-header -- parallel --pipe 'bash filter_result.sh' \
        > 3_validating/tmp.filter.result.tsv

    tsv-join \
        3_validating/tmp.pass.result.tsv \
        -d 1 \
        -f 3_validating/tmp.filter.result.tsv \
        -k 1 \
        > 3_validating/tmp.result.tsv

    mv 3_validating/tmp.result.tsv 3_validating/tmp.pass.result.tsv
    wc -l 3_validating/tmp.pass.result.tsv
done | tee -a 3_validating/validation.log

mv 3_validating/tmp.pass.result.tsv 3_validating/result.tsv
rm 3_validating/tmp.*.tsv

bash result_stat.sh 3_validating

```

| #Item                   | Value   |
|:------------------------|:--------|
| 3_validating/result.tsv | 14186   |
| count                   | 1797    |
| hr_max                  | 4.83823 |
| hr_median               | 3.01711 |
| kmp_min                 | 0.00000 |
| kmp_median              | 0       |
| rocauc_max              | 0.86246 |
| rocauc_median           | 0.75881 |
| rocp_min                | 0.00000 |
| rocp_median             | 0       |

# 4_non_rare

```shell
RARE_PASS=2

mkdir -p 4_non_rare

# Count markers
cat 3_validating/result.tsv |
    cut -f 1 |
    perl -nl -e 'print for split /\+/' |
    tsv-uniq --number |
    tsv-summarize --group-by 1 --max 2 |
    sort -k2,2nr \
    > 4_non_rare/count.tsv

# remove rares
head -n 1 3_validating/result.tsv > 4_non_rare/result.tsv
cat 3_validating/result.tsv |
    grep -F -f <(tsv-filter --le 2:${RARE_PASS} 4_non_rare/count.tsv | cut -f 1 ) -w -v \
    >> 4_non_rare/result.tsv

cat 4_non_rare/result.tsv |
    cut -f 1 |
    perl -nl -e 'print for split /\+/' |
    tsv-uniq --number |
    tsv-summarize --group-by 1 --max 2 |
    sort -k2,2nr \
    > 4_non_rare/non-rare.count.tsv

bash result_stat.sh 4_non_rare

```

| #Item                         | Value    |
|:------------------------------|:---------|
| 4_non_rare/count.tsv          | 1798     |
| 4_non_rare/non-rare.count.tsv | 1309     |
| 4_non_rare/result.tsv         | 13499    |
| count                         | 1308     |
| hr_max                        | 4.83823  |
| hr_median                     | 3.018555 |
| kmp_min                       | 0.00000  |
| kmp_median                    | 0        |
| rocauc_max                    | 0.86246  |
| rocauc_median                 | 0.75899  |
| rocp_min                      | 0.00000  |
| rocp_median                   | 0        |

# 4_combination

```shell
mkdir -p 4_combination

# Count markers
cat 4_non_rare/result.tsv |
    cut -f 1 |
    perl -nl -e 'print for split /\+/' |
    tsv-uniq --number |
    tsv-summarize --group-by 1 --max 2 \
    > 4_combination/count.tsv

tsv-join \
    4_combination/count.tsv \
    -d 1 \
    -f 2_training/data.replace.tsv \
    -k 1 \
    -a 2 |
    tsv-filter --ge 2:5 |
    tsv-join \
        -d 3 \
        -f site.tsv \
        -k 1 \
        -a 2-10 |
    sort -k2,2nr \
    > 4_combination/sites.tsv

# count occurrences of genes
cat 4_combination/sites.tsv |
    parallel -j 1 -k --pipe '
        cut -f 7 |
            perl -nla -F"\;" -MList::Util -e '\''
                @symbols = List::Util::uniq (sort @F);
                print for @symbols;
            '\'' |
            grep -v "\."
    ' |
    tsv-uniq --at-least 2 --number |
    tsv-summarize --group-by 1 --max 2 |
    tsv-filter --ge 2:3 \
    > 4_combination/interest_occur.tsv
cat 4_combination/interest_occur.tsv

# count occurrences of gene families
cat 4_combination/sites.tsv |
    parallel -j 1 -k --pipe '
        cut -f 7 |
            perl -nla -F"\;" -MList::Util -e '\''
                @symbols = List::Util::uniq (sort @F);
                s/\d+\w?$// for @symbols;
                print for @symbols;
            '\'' |
            grep -v "\."
    ' |
    tsv-uniq --at-least 2 --number |
    tsv-summarize --group-by 1 --max 2 |
    tsv-filter --ge 2:4 \
    > 4_combination/interest_family.tsv
cat 4_combination/interest_family.tsv

# Cancer Genes
cat 4_combination/sites.tsv |
    parallel -j 1 -k --pipe '
        cut -f 7 |
            perl -nla -F"\;" -MList::Util -e '\''
                @symbols = List::Util::uniq (sort @F);
                print for @symbols;
            '\'' |
            grep -v "\."
    ' |
    tsv-uniq |
    (echo "Gene Symbol" && cat) |
    tsv-join \
        -H -d 1 \
        -f Cancer_Gene_Census.tsv \
        -k 1 \
        -a 2,5,6,8-11,15 \
    > 4_combination/interest_cgc.tsv
cat 4_combination/interest_cgc.tsv

# Full combinations
tsv-filter --ge 2:100 4_combination/sites.tsv |
    tsv-filter --str-ne 7:. |
    cut -f 1 |
    sort > 4_combination/step1.tsv

cat 4_combination/sites.tsv |
    grep -Fw -f <(cut -f 1 4_combination/interest_cgc.tsv) |
    tsv-filter --lt 2:100 --ge 2:10 \
    >> 4_combination/step1.tsv

tsv-join \
    4_combination/sites.tsv -d 1 \
    -f 4_combination/step1.tsv -k 1 |
    cut -f 3 |
    sort \
    > 4_combination/sites.step1.tsv

wc -l 4_combination/step1.tsv
#87 4_combination/step1.tsv

bmr nextstep 4_combination/step1.tsv > 4_combination/step2.tsv
bmr nextstep 4_combination/step2.tsv 4_combination/step1.tsv |
    tsv-uniq > 4_combination/step3.tsv
bmr nextstep 4_combination/step3.tsv 4_combination/step1.tsv |
    tsv-uniq > 4_combination/step4.tsv
bmr nextstep 4_combination/step4.tsv 4_combination/step1.tsv |
    tsv-uniq > 4_combination/step5.tsv

rm 4_combination/step1.tsv
cat 4_combination/step*.tsv |
    grep -v '#' |
    (echo '#marker' && cat) |
    tsv-uniq > 4_combination/formula.tsv
rm 4_combination/step*.tsv

# training
bmr split 4_combination/formula.tsv -c 10000 --mode row --rr 1 -o split | bash

mkdir -p job
find split -type f -name "*[0-9]" |
    sort |
    split -l 120 -a 3 -d - job/

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo ${f}
    bsub -q mpi -n 24 -J "training-${f}" \
        "
        cat ${f} |
            parallel --no-run-if-empty --line-buffer -k -j 24 '
                echo '\''==> Processing {}'\''
                Rscript multivariate.R -i 2_training/data.tsv -f {} -o {}.result.tsv
            '
        "
done

#find split -type f -name "*[0-9]" |
#    sort |
#    parallel --no-run-if-empty --line-buffer -k -j 20 '
#        echo "==> Processing {}"
#        Rscript multivariate.R -i 2_training/data.tsv -f {} -o {}.result.tsv
#    '

tsv-append -H split/*.result.tsv > 4_combination/cox.result.tsv

rm -fr split job BS output.*

bmr cox \
    --coxp 0.01 \
    --hr 0.01,3 \
    --kmp 0.01 \
    --rocauc 0.1,0.80 \
    --rocp 0.01

cat 4_combination/cox.result.tsv |
    keep-header -- parallel --pipe 'bash filter_result.sh' \
    > 4_combination/training.result.tsv

# testing
bmr cox \
    --coxp 0.05 \
    --hr 0.01,2.5 \
    --kmp 0.01 \
    --rocauc 0.01,0.75 \
    --rocp 0.01

cat 4_combination/training.result.tsv |
    sed '1d' |
    parallel --pipe  --block 200k -j 20 'Rscript validate.R -i 2_testing/data.tsv -f stdin --no' |
    keep-header -- grep -v '^#' \
    > 4_combination/validate.result.tsv

cat 4_combination/validate.result.tsv |
    keep-header -- parallel --pipe 'bash filter_result.sh' \
    > 4_combination/testing.result.tsv

tsv-join \
    4_combination/training.result.tsv \
    -d 1 \
    -f 4_combination/testing.result.tsv \
    -k 1 \
    > 4_combination/pass.result.tsv

# clinical
# --rocp set to 0 as drug_target may fall back to pROC
bmr cox \
    --coxp 0.05 \
    --hr 0.01,2 \
    --kmp 0.05 \
    --rocauc 0.01,0.65 \
    --rocp 0

cp 4_combination/pass.result.tsv 4_combination/tmp.pass.result.tsv
for c in \
    MALE FEMALE \
    Stage_I Stage_II Stage_III_IV \
    Age_59 Age_60_69 Age_70 \
    egfr_yes egfr_no \
    drug_target drug_platin \
    chemo_yes chemo_no \
    radiation_yes radiation_no \
    ; do
    echo "==> ${c}"
    wc -l 4_combination/tmp.pass.result.tsv

    cat 4_combination/tmp.pass.result.tsv |
        sed '1d' |
        parallel --pipe --block 200k -j 20 "Rscript validate.R -i clinical/${c}/data.tsv -f stdin --no" |
        keep-header -- grep -v '^#' |
        keep-header -- parallel --pipe 'bash filter_result.sh' \
        > 4_combination/tmp.filter.result.tsv

    tsv-join \
        4_combination/tmp.pass.result.tsv \
        -d 1 \
        -f 4_combination/tmp.filter.result.tsv \
        -k 1 \
        > 4_combination/tmp.result.tsv

    mv 4_combination/tmp.result.tsv 4_combination/tmp.pass.result.tsv
    wc -l 4_combination/tmp.pass.result.tsv
done | tee 4_combination/clinical.log

mv 4_combination/tmp.pass.result.tsv 4_combination/clinical.result.tsv
rm 4_combination/tmp.*.tsv

# extra validation
# --rocp set to 0 as XH33 may fall back to pROC
bmr cox \
    --coxp 0.05 \
    --hr 0.01,2 \
    --kmp 0.05 \
    --rocauc 0.01,0.65 \
    --rocp 0

rm 4_combination/validation.log

cp 4_combination/clinical.result.tsv 4_combination/tmp.pass.result.tsv
for c in \
    GSE39279 \
    ; do
    echo "==> ${c}"
    wc -l 4_combination/tmp.pass.result.tsv

    cat 4_combination/tmp.pass.result.tsv |
        sed '1d' |
        parallel --pipe --block 200k -j 20 "Rscript validate.R -i validating/${c}/data.tsv -f stdin --no" |
        keep-header -- grep -v '^#' |
        keep-header -- parallel --pipe 'bash filter_result.sh' \
        > 4_combination/tmp.filter.result.tsv

    tsv-join \
        4_combination/tmp.pass.result.tsv \
        -d 1 \
        -f 4_combination/tmp.filter.result.tsv \
        -k 1 \
        > 4_combination/tmp.result.tsv

    mv 4_combination/tmp.result.tsv 4_combination/tmp.pass.result.tsv
    wc -l 4_combination/tmp.pass.result.tsv
done | tee -a 4_combination/validation.log

for c in \
    XH33 \
    ; do
    echo "==> ${c}"
    wc -l 4_combination/tmp.pass.result.tsv

    cat 4_combination/tmp.pass.result.tsv |
        sed '1d' |
        parallel --pipe --block 200k -j 20 "Rscript validate.R -t 1095 -i validating/${c}/data.tsv -f stdin --no" |
        keep-header -- grep -v '^#' |
        keep-header -- parallel --pipe 'bash filter_result.sh' \
        > 4_combination/tmp.filter.result.tsv

    tsv-join \
        4_combination/tmp.pass.result.tsv \
        -d 1 \
        -f 4_combination/tmp.filter.result.tsv \
        -k 1 \
        > 4_combination/tmp.result.tsv

    mv 4_combination/tmp.result.tsv 4_combination/tmp.pass.result.tsv
    wc -l 4_combination/tmp.pass.result.tsv
done | tee -a 4_combination/validation.log

#for c in \
#    GSE119144_PFS GSE119144_benefit \
#    ; do
#    echo "==> ${c}"
#    wc -l 4_combination/tmp.pass.result.tsv
#
#    cat 4_combination/tmp.pass.result.tsv |
#        sed '1d' |
#        parallel --pipe --block 200k -j 20 "Rscript validate.R -t 180 -i validating/${c}/data.tsv -f stdin --no" |
#        keep-header -- grep -v '^#' |
#        keep-header -- parallel --pipe 'bash filter_result.sh' \
#        > 4_combination/tmp.filter.result.tsv
#
#    tsv-join \
#        4_combination/tmp.pass.result.tsv \
#        -d 1 \
#        -f 4_combination/tmp.filter.result.tsv \
#        -k 1 \
#        > 4_combination/tmp.result.tsv
#
#    mv 4_combination/tmp.result.tsv 4_combination/tmp.pass.result.tsv
#    wc -l 4_combination/tmp.pass.result.tsv
#done | tee -a 4_combination/validation.log

mv 4_combination/tmp.pass.result.tsv 4_combination/validating.result.tsv
rm 4_combination/tmp.*.tsv

# bootstrap
BS_PASS=80

bmr cox \
    --coxp 0.01 \
    --hr 0.33,3 \
    --kmp 0.01 \
    --rocauc 0.25,0.75 \
    --rocp 0.01

bash BS.sh 2_training/data.tsv

bash bootstrap.sh BS 4_combination/validating.result.tsv 4_combination

rm -fr split job BS output.*

mv 4_combination/validating.result.tsv.bootstrap.tsv 4_combination/bootstrap.result.tsv

cat 4_combination/bootstrap.result.tsv |
    tsv-filter -H --ge 9:${BS_PASS} |
    tsv-join \
        -H -f 4_combination/testing.result.tsv \
        -k 1 -a 5,6,7,8 --prefix testing- |
    keep-header -- perl -nla -F"\t" -e '
        if ($F[4] >=1) {
            print if $F[9] >= 1;
        } else {
            print if $F[9] < 1;
        }
    ' |
    keep-header -- perl -nla -F"\t" -e '
        if ($F[6] >=0.5) {
            print if $F[11] >= 0.5;
        } else {
            print if $F[11] < 0.5;
        }
    ' \
    > 4_combination/result.tsv

bash result_stat.sh -b 4_combination
#bash result_dist.sh 4_combination/cox.result.tsv
#bash result_dist.sh 4_combination/validate.result.tsv

# attach univariate results
cat 4_combination/result.tsv |
    cut -f 1 |
    perl -nl -e 'print for split /\+/' |
    sort |
    uniq -c |
    perl -nl -e 's/^(\s*)(\d*)(\s*)(.*)/$4\t$2/g; print $_' |
    sort -k2,2nr |
    tsv-join \
        -f 4_combination/sites.tsv \
        -k 1 -a 3,11,12,7-8 |
    perl -nla -F"\t" -MList::Util -e '
        @symbols = List::Util::uniq split /;/, $F[5];
        $F[5] = join qq{;}, @symbols;
        @functions = List::Util::uniq split /;/, $F[6];
        $F[6] = join qq{;}, @functions;
        print join qq{\t}, @F;
    ' \
    > 4_combination/annotation.tsv

cat 1_testing/join.result.tsv |
    sed '1d' |
    tsv-join -d 1 -k 3 -f 4_combination/annotation.tsv -a 1 |
    tsv-select -f 14,2-13 \
    >> 4_combination/result.tsv

cat gencode.tsv |
    tsv-join -H \
        -d 2 \
        -f <(
            tsv-select 4_combination/annotation.tsv -f 6 |
                tr ';' '\n' |
                tsv-uniq |
                (echo SYMBOL && cat)
        ) \
        -k 1 \
    > 4_combination/gencode.tsv

# XH33 individuals, single blinded
#Rscript individual.R -i validating/XH33/data.tsv -f 4_combination/result.tsv \
#    > 4_combination/XH33.individuals.tsv

```

| #Item                                         | Value    |
|:----------------------------------------------|:---------|
| 4_combination/annotation.tsv                  | 82       |
| 4_combination/bootstrap.result.tsv            | 765      |
| 4_combination/clinical.result.tsv             | 7287     |
| 4_combination/count.tsv                       | 1309     |
| 4_combination/cox.result.tsv                  | 39285489 |
| 4_combination/formula.tsv                     | 39285489 |
| 4_combination/gencode.tsv                     | 81       |
| 4_combination/interest_cgc.tsv                | 43       |
| 4_combination/interest_family.tsv             | 33       |
| 4_combination/interest_occur.tsv              | 39       |
| 4_combination/pass.result.tsv                 | 22640    |
| 4_combination/result.tsv                      | 438      |
| 4_combination/sites.step1.tsv                 | 87       |
| 4_combination/sites.tsv                       | 1042     |
| 4_combination/testing.result.tsv              | 22640    |
| 4_combination/training.result.tsv             | 74291    |
| 4_combination/validate.result.tsv             | 74291    |
| 4_combination/validating.result.tsv           | 765      |
| 4_combination/validating.result.tsv.count.tsv | 765      |
| count                                         | 77       |
| hr_max                                        | 5.58141  |
| hr_median                                     | 3.96253  |
| kmp_min                                       | 0.00000  |
| kmp_median                                    | 0        |
| rocauc_max                                    | 0.88173  |
| rocauc_median                                 | 0.82264  |
| rocp_min                                      | 0.00000  |
| rocp_median                                   | 0        |
| bin(BS)                                       | count    |
| 100                                           | 12       |
| 80                                            | 425      |
| 60                                            | 274      |
| 40                                            | 53       |

# 4_figure

It should be run locally.

```shell
rm -fr 4_figure

parallel --no-run-if-empty --linebuffer -k -j 4 '
    SEQ=$(printf "%02d" {#})
    mkdir -p 4_figure/${SEQ}
    Rscript kmroc.R \
        -i {} \
        -f 4_combination/result.tsv \
        -o 4_figure/${SEQ}
    ' ::: \
    2_training/data.tsv 2_testing/data.tsv \
    clinical/MALE/data.tsv clinical/FEMALE/data.tsv \
    clinical/Stage_I/data.tsv clinical/Stage_II/data.tsv clinical/Stage_III_IV/data.tsv \
    clinical/Age_59/data.tsv clinical/Age_60_69/data.tsv clinical/Age_70/data.tsv \
    clinical/egfr_yes/data.tsv clinical/egfr_no/data.tsv \
    clinical/chemo_yes/data.tsv clinical/chemo_no/data.tsv \
    clinical/radiation_yes/data.tsv clinical/radiation_no/data.tsv \
    clinical/drug_target/data.tsv clinical/drug_platin/data.tsv \
    clinical/drug_mitotic/data.tsv clinical/drug_metabolites/data.tsv \
    validating/GSE39279/data.tsv

find plot/01 -maxdepth 1 -type f -name "*.pdf" |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {/}
        gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=plot/{/} \
            $(find plot -maxdepth 2 -mindepth 2 -type f -name {/} | sort)
        ../../.TinyTeX/bin/x86_64-linux/pdfjam plot/{/} --nup 3x11 --suffix nup -o plot
        mv plot/{/.}-nup.pdf plot/{/}
    '

find plot/ -maxdepth 1 -mindepth 1 -type d | xargs rm -fr

```

# 4_figure_GSE119144

```shell
rm -fr 4_figure_GSE119144

parallel --no-run-if-empty --linebuffer -k -j 4 '
    SEQ=$(printf "%02d" {#})
    mkdir -p 4_figure_GSE119144/${SEQ}
    Rscript kmroc.R \
        -i {} \
        -f 4_combination/result.tsv \
        -t 180 \
        -o 4_figure_GSE119144/${SEQ}
    ' ::: \
    validating/GSE119144_PFS/data.tsv validating/GSE119144_benefit/data.tsv

find 4_figure_GSE119144/01 -maxdepth 1 -type f -name "*.pdf" |
    parallel --no-run-if-empty --linebuffer -k -j 4 '
        echo {/}
        gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=4_figure_GSE119144/{/} \
            $(find 4_figure_GSE119144 -maxdepth 2 -mindepth 2 -type f -name {/} | sort)
        pdfjam 4_figure_GSE119144/{/} --nup 2x1 --suffix nup -o 4_figure_GSE119144
        mv 4_figure_GSE119144/{/.}-nup.pdf 4_figure_GSE119144/{/}
    '

find 4_figure_GSE119144/ -maxdepth 1 -mindepth 1 -type d | xargs rm -fr

```

# Top 300

Top 300 loci

```shell
rm -fr top_300
mkdir -p top_300

cat 4_combination/sites.tsv |
    tsv-filter --str-ne 7:. |
    head -n 300 |
    cut -f 3-7 |
    sort -k2,2 -k3,3n \
    > top_300/list.tsv

bash select_col.sh -f 1-3 meth.tsv.gz top_300/list.tsv > top_300/data.tsv

sed '1d' top_300/data.tsv | datamash check
# 449 lines, 303 fields

Rscript univariate.R --valid 0.99 -i top_300/data.tsv -o top_300/result.tsv

rm -fr top_300/figure
mkdir -p top_300/figure
Rscript kmroc.R \
    -s \
    -i top_300/data.tsv \
    -f top_300/result.tsv \
    -o top_300/figure

mkdir -p top_300/split
sed '1d' top_300/result.tsv |
    cut -f 1 |
    split -l 10 -a 3 - top_300/split/

mkdir -p top_300/merge
find top_300/split/ -type f |
    sort |
    parallel --no-run-if-empty --linebuffer -k -j 1 '
        echo {/}
        gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=top_300/merge/{/}.pdf \
            $(
            cat {} |
                sed '\''s/.*/top_300\/figure\/&.pdf/'\''
            )
        pdfjam top_300/merge/{/}.pdf --nup 2x5 --suffix nup -o top_300/merge
        mv top_300/merge/{/}-nup.pdf top_300/merge/{/}.pdf
    '

rm -fr top_300/figure top_300/split

# selected one
cat <<EOF > top_300/selected.tsv
#marker
cg03046247+cg03356689+cg07109046+cg08813325+cg11770080
EOF

Rscript multivariate.R -i top_300/data.tsv \
    -f top_300/selected.tsv -o top_300/selected.result.tsv

mkdir -p top_300/selected
Rscript kmroc.R \
    -s \
    -i top_300/data.tsv \
    -f top_300/selected.result.tsv \
    -o top_300/selected

TAB=$'\t'
cat <<EOF > top_300/selected.1.tsv
#marker${TAB}coef
cg03046247+cg03356689+cg07109046+cg08813325+cg11770080${TAB}1,1,1,-1,1
EOF

Rscript quantile.R \
    -i top_300/data.tsv \
    -f top_300/selected.1.tsv \
    -q 0,50,100 \
    -o top_300/
# median 1.07453

cat <<EOF > top_300/selected.1.tsv
#marker${TAB}coef${TAB}p${TAB}median${TAB}hr
cg03046247+cg03356689+cg07109046+cg08813325+cg11770080${TAB}1,1,1,-1,1${TAB}${TAB}1.07453${TAB}
EOF

mkdir -p top_300/selected-1
Rscript kmroc.R \
    -s \
    -i top_300/data.tsv \
    -f top_300/selected.1.tsv \
    -o top_300/selected-1

```

# Interesting loci

Looking for signatures of ecDNA

```shell
rm -fr loci
mkdir -p loci

cat 4_combination/annotation.tsv |
    tsv-filter --or \
        --istr-in-fld 6:FAM53B \
        --istr-in-fld 6:HMGA2 \
        --istr-in-fld 6:KATNAL1 \
        --istr-in-fld 6:HMBOX1 \
        --istr-in-fld 6:BHMT \
        --istr-in-fld 6:DMGDH \
        --istr-in-fld 6:SORCS2 \
        --istr-in-fld 6:EGFR \
        --istr-in-fld 6:ETS1 \
        --istr-in-fld 6:FLI1 |
    tsv-select -f 1,3,4,6 |
    perl -nla -F'\t' -e '
        $F[2] =~ s/CGI://g;
        $F[2] =~ s/[:\-]/\t/g;
        ($F[3]) = split /;/, $F[3];
        print join qq{\t}, @F;
    ' |
    sort -k3,3 -k6,6 \
    > loci/loci.tsv
cat loci/loci.tsv
#m0721	cg07081759	chr10	124615225	124615644	FAM53B
#m1598	cg16404170	chr10	124615225	124615644	FAM53B
#m2428	cg26508444	chr10	124615225	124615644	FAM53B
#m2294	cg24662823	chr11	128451766	128452000	ETS1
#m0345	cg03356689	chr11	128692776	128695170	FLI1
#m0888	cg08825225	chr11	128692776	128695170	FLI1
#m1205	cg12141052	chr12	65823737	65826401	HMGA2
#m2146	cg23146197	chr12	65823737	65826401	HMGA2
#m0860	cg08538509	chr13	30114858	30115198	KATNAL1
#m0884	cg08813325	chr4	7664292	7664620	SORCS2
#m2512	cg27488680	chr4	7646028	7646233	SORCS2
#m0251	cg02286091	chr5	79069475	79069888	BHMT
#m0606	cg05890484	chr5	79069475	79069888	BHMT
#m1066	cg10660256	chr5	79069475	79069888	BHMT
#m1159	cg11770080	chr5	79069475	79069888	BHMT
#m1347	cg13695646	chr5	79069475	79069888	BHMT
#m2220	cg23987322	chr5	79069475	79069888	BHMT
#m0316	cg03046247	chr7	55018448	55020412	EGFR
#m1172	cg11893955	chr8	29071554	29072201	HMBOX1

# sites around loci
cat loci/loci.tsv |
    parallel --colsep '\t' -j 1 -k '
        lower=$(perl -MPOSIX -e "print 1_000_000 * (POSIX::floor({4} / 1_000_000) - 2)")
        upper=$(perl -MPOSIX -e "print 1_000_000 * (POSIX::ceil({5} / 1_000_000) + 2)")

        >&2 echo {6} {4} {5} ${lower} ${upper}

        cat site.tsv |
            tsv-filter --istr-eq 2:{3} --ge 3:${lower} --le 3:${upper}  |
            sort -k3,3n
    ' |
    tsv-uniq \
    > loci/sites.tsv

bash select_col.sh -f 1-3 meth.tsv.gz loci/sites.tsv > loci/meth.tsv

Rscript univariate.R -p 8 -i loci/meth.tsv -o loci/meth.result.tsv

# clinical
for g in \
    Stage_I Stage_II Stage_III_IV \
    egfr_yes egfr_no \
    drug_target drug_platin \
    ; do \
    echo "==> ${g}"
    cat loci/meth.tsv |
        tsv-join \
            -H -k 1 \
            -f clinical/${g}/data.tsv\
        > loci/${g}.data.tsv
done

for g in \
    Stage_I Stage_II Stage_III_IV \
    egfr_yes egfr_no \
    drug_target drug_platin \
    ; do \
    echo "==> ${g}"
    Rscript validate.R -p 8 -i loci/${g}.data.tsv \
        -f loci/meth.result.tsv -o loci/${g}.result.tsv
done

# join
cp loci/meth.result.tsv loci/tmp.join.result.tsv
for g in \
    Stage_I Stage_II Stage_III_IV \
    egfr_yes egfr_no \
    drug_target drug_platin \
    ; do
    echo "==> ${g}"

    tsv-join \
        loci/tmp.join.result.tsv \
        -H -k 1 \
        -f loci/${g}.result.tsv \
        -a 6,7 --prefix "${g}-" \
        > loci/tmp.result.tsv

    mv loci/tmp.result.tsv loci/tmp.join.result.tsv
done

cat loci/tmp.join.result.tsv |
    tsv-join \
        --write-all HEAD -k 1 \
        -f loci/sites.tsv \
        -a 2,3,5 |
    keep-header -- sort -k23,23 -k24,24n \
    > loci/join.result.tsv

rm loci/tmp.*

# excel
plotr tsv loci/join.result.tsv \
    --header --font_size 11 \
    $(
    cat loci/loci.tsv |
        parallel --colsep '\t' -j 1 -k '
            echo "--contain 25:{6} --contain 1:{2}"
        '
    ) \
    --le 3:0.05  --le 6:0.05 --ge 7:0.65  \
    --le 9:0.05  --le 11:0.05 --le 13:0.05 --le 15:0.05 --le 17:0.05 --le 19:0.05 --le 21:0.05 \
    --ge 10:0.65 --ge 12:0.65 --ge 14:0.65 --ge 16:0.65 --ge 18:0.65 --ge 20:0.65 --ge 22:0.65

```

# 4_loci

```shell
mkdir -p 4_loci

cut -d' ' -f 1 loci/loci.tsv > 4_loci/step1.tsv

wc -l 4_loci/step1.tsv
#19 4_loci/step1.tsv

bmr nextstep 4_loci/step1.tsv > 4_loci/step2.tsv
bmr nextstep 4_loci/step2.tsv 4_loci/step1.tsv |
    tsv-uniq > 4_loci/step3.tsv
bmr nextstep 4_loci/step3.tsv 4_loci/step1.tsv |
    tsv-uniq > 4_loci/step4.tsv
bmr nextstep 4_loci/step4.tsv 4_loci/step1.tsv |
    tsv-uniq > 4_loci/step5.tsv
bmr nextstep 4_loci/step5.tsv 4_loci/step1.tsv |
    tsv-uniq > 4_loci/step6.tsv
bmr nextstep 4_loci/step6.tsv 4_loci/step1.tsv |
    tsv-uniq > 4_loci/step7.tsv
bmr nextstep 4_loci/step7.tsv 4_loci/step1.tsv |
    tsv-uniq > 4_loci/step8.tsv
bmr nextstep 4_loci/step8.tsv 4_loci/step1.tsv |
    tsv-uniq > 4_loci/step9.tsv
bmr nextstep 4_loci/step9.tsv 4_loci/step1.tsv |
    tsv-uniq > 4_loci/step10.tsv

rm 4_loci/step1.tsv
cat 4_loci/step*.tsv |
    grep -v '#' |
    (echo '#marker' && cat) |
    tsv-uniq > 4_loci/formula.tsv
rm 4_loci/step*.tsv

# training
bmr split 4_loci/formula.tsv -c 10000 --mode row --rr 1 -o split | bash

find split -type f -name "*[0-9]" |
    sort |
    parallel --no-run-if-empty --line-buffer -k -j 20 '
        echo "==> Processing {}"
        Rscript multivariate.R -i 2_training/data.tsv -f {} -o {}.result.tsv
    '

tsv-append -H split/*.result.tsv \
    > 4_loci/cox.result.tsv

rm -fr split job BS output.*

bmr cox \
    --coxp 0.05 \
    --hr 0.33,3 \
    --kmp 0.05 \
    --rocauc 0.20,0.80 \
    --rocp 0.5

cat 4_loci/cox.result.tsv |
    keep-header -- parallel --pipe 'bash filter_result.sh' \
    > 4_loci/training.result.tsv

# testing
bmr cox \
    --coxp 0.05 \
    --hr 0.4,2.5 \
    --kmp 0.05 \
    --rocauc 0.25,0.75 \
    --rocp 0.5

cat 4_loci/training.result.tsv |
    sed '1d' |
    parallel --pipe -j 20 'Rscript validate.R -i 2_testing/data.tsv -f stdin --no' |
    keep-header -- grep -v '^#' \
    > 4_loci/validate.result.tsv

cat 4_loci/validate.result.tsv |
    keep-header -- parallel --pipe 'bash filter_result.sh' \
    > 4_loci/testing.result.tsv

tsv-join \
    4_loci/training.result.tsv \
    -d 1 \
    -f 4_loci/testing.result.tsv \
    -k 1 \
    > 4_loci/result.tsv

bash result_stat.sh 4_loci

```

| #Item                      | Value   |
|:---------------------------|:--------|
| 4_loci/cox.result.tsv      | 354503  |
| 4_loci/formula.tsv         | 354503  |
| 4_loci/result.tsv          | 200     |
| 4_loci/testing.result.tsv  | 200     |
| 4_loci/training.result.tsv | 704     |
| 4_loci/validate.result.tsv | 704     |
| count                      | 19      |
| hr_max                     | 5.49804 |
| hr_median                  | 3.71991 |
| kmp_min                    | 0.00000 |
| kmp_median                 | 0       |
| rocauc_max                 | 0.87556 |
| rocauc_median              | 0.82759 |
| rocp_min                   | 0.02000 |
| rocp_median                | 0.16    |

# Pack up

```shell
cd ~/data/cancer/

tar cvfz LUAD.$(date +"%Y%m%d").tar.gz \
    LUAD/4_figure/ \
    LUAD/loci/join.result.tsv.xlsx \
    LUAD/4_combination/result.tsv \
    LUAD/4_combination/annotation.tsv \
    LUAD/4_combination/sites.tsv \
    LUAD/4_combination/gencode.tsv \
    LUAD/meth_clinical.tsv \
    LUAD/drug_values.tsv \
    LUAD/radiation_values.tsv

```

