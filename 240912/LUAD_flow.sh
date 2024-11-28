sort ours.probe.list others.probe.list | sort | uniq >allsite.list
bash select_col.sh -f 1-3 training.tsv.gz allsite.list >training/data.tsv
bash select_col.sh -f 1-3 testing.tsv.gz allsite.list >testing/data.tsv

for g in \
	MALE FEMALE \
	Stage_I Stage_II Stage_III_IV \
	Age_59 Age_60_69 Age_70 \
	egfr_yes egfr_no \
	drug_egfr drug_target drug_platin drug_mitotic drug_metabolites; do
	echo "==> ${g}"
	mkdir -p clinical/${g}
	tsv-append -H training/data.tsv testing/data.tsv \
		| tsv-join \
			-H -k 1 \
			-f clinical/${g}.tsv \
			>clinical/${g}/data.tsv
	datamash check <clinical/${g}/data.tsv
done
#==> MALE
#212 lines, 164 fields
#==> FEMALE
#239 lines, 164 fields
#==> Stage_I
#245 lines, 164 fields
#==> Stage_II
#111 lines, 164 fields
#==> Stage_III_IV
#91 lines, 164 fields
#==> Age_59
#127 lines, 164 fields
#==> Age_60_69
#144 lines, 164 fields
#==> Age_70
#171 lines, 164 fields
#==> egfr_yes
#81 lines, 164 fields
#==> egfr_no
#176 lines, 164 fields
#==> drug_egfr
#23 lines, 164 fields
#==> drug_target
#32 lines, 164 fields
#==> drug_platin
#113 lines, 164 fields
#==> drug_mitotic
#70 lines, 164 fields
#==> drug_metabolites
#54 lines, 164 fields

# chemo
mkdir -p clinical/chemo_yes clinical/chemo_no
tsv-append -H training/data.tsv testing/data.tsv \
	| tsv-join \
		-H -k 1 \
		-f drug_values.tsv \
		>clinical/chemo_yes/data.tsv
tsv-append -H training/data.tsv testing/data.tsv \
	| tsv-join \
		-H -k 1 \
		-f drug_values.tsv \
		--exclude \
		>clinical/chemo_no/data.tsv
datamash check <clinical/chemo_yes/data.tsv
datamash check <clinical/chemo_no/data.tsv
#160 lines, 164 fields
#291 lines, 164 fields

# radiation
tsv-summarize -H --group-by 1 --unique-values 3 \
	<radiation.tsv \
	>radiation_values.tsv

mkdir -p clinical/radiation_yes clinical/radiation_no
tsv-append -H training/data.tsv testing/data.tsv \
	| tsv-join \
		-H -k 1 \
		-f radiation_values.tsv \
		>clinical/radiation_yes/data.tsv
tsv-append -H training/data.tsv testing/data.tsv \
	| tsv-join \
		-H -k 1 \
		-f radiation_values.tsv \
		--exclude \
		>clinical/radiation_no/data.tsv
datamash check <clinical/radiation_yes/data.tsv
datamash check <clinical/radiation_no/data.tsv
#89 lines, 164 fields
#362 lines, 164 fields
