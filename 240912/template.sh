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
