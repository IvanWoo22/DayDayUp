for i in NJU-BSP000-{1..4}; do
	mkdir "${i}"
	ln -sf /home/ivan/cater_data/NJU-BSP000/"${i}"*_R1_*.fastq.gz "${i}"/R1.fq.gz
	ln -sf /home/ivan/cater_data/NJU-BSP000/"${i}"*_R2_*.fastq.gz "${i}"/R2.fq.gz
	cutadapt -a AGATCGGAAGAGCACA -A AGATCGGAAGAGCGTC \
		-O 5 -m 10 -e 0.18 -j 8 \
		-o "${i}"/clean_R1.fq.gz -p "${i}"/clean_R2.fq.gz \
		"${i}"/R1.fq.gz "${i}"/R2.fq.gz >"${i}"/cutadapt.log \
		2>&1
done

parallel --keep-order --xapply -j 4 \
	'pear -j 4 -n 150 -f {}/clean_R1.fq.gz -r {}/clean_R2.fq.gz -o {}/merge' \
	::: FXH139 FXH149 FXH22 FXH28 FXH36 FXH37 FXH4 FXH42 FXH45 FXH5 FXH52 FXH63 FXH67 FXH68 FXH69 FXH81 FXH83 FXH86 FXH87 FXH88 FXH92 FXH95 FXH96 FXHN139 FXHN149 FXHN2 FXHN22 FXHN28 FXHN36 FXHN37 FXHN4 FXHN42 FXHN45 FXHN5 FXHN52 FXHN68 FXHN83 FXHN86 FXHN88 FXHN96

rm ./*/merge.discarded.fastq ./*/merge.unassembled.forward.fastq ./*/merge.unassembled.reverse.fastq

for i in FXH139 FXH149 FXH22 FXH28 FXH36 FXH37 FXH4 FXH42 FXH45 FXH5 FXH52 FXH63 FXH67 FXH68 FXH69 FXH81 FXH83 FXH86 FXH87 FXH88 FXH92 FXH95 FXH96 FXHN139 FXHN149 FXHN2 FXHN22 FXHN28 FXHN36 FXHN37 FXHN4 FXHN42 FXHN45 FXHN5 FXHN52 FXHN68 FXHN83 FXHN86 FXHN88 FXHN96; do
	sed '/^@/s/ /_/g' "${i}"/merge.assembled.fastq \
		| umi_tools extract \
			--umi-separator=: --extract-method=string --3prime \
			--bc-pattern=NNNNNN \
			-L "${i}"/umi_extract.log \
			-S "${i}"/merge.umi.fastq

	rm "${i}"/merge.assembled.fastq

	/home/ivan/bismark/bismark \
		--genome /home/ivan/fat/index/hg38_bs_index/ \
		--non_directional -q -p 4 \
		"${i}"/merge.umi.fastq \
		--gzip --sam -B "${i}"/bsk

	/home/ivan/bismark/bismark_methylation_extractor \
		--parallel 6 --bedGraph --counts --report --comprehensive \
		--genome_folder /home/ivan/fat/index/hg38_bs_index/ \
		--output_dir "${i}"/. \
		"${i}"/bsk.sam
done

rm ./*/bsk.sam ./*/bsk.deduplicated.sam ./*/bsk.deduplicated.CpG_report.txt ./*/*_context_bsk.deduplicated.txt
pigz ./*/merge.umi.fastq

/home/ivan/bismark/deduplicate_bismark \
	--sam --barcode "${i}"/bsk.sam \
	--output_dir "${i}"/.

/home/ivan/bismark/bismark_methylation_extractor \
	--parallel 6 --bedGraph --counts --report --cytosine_report --comprehensive \
	--genome_folder /home/ivan/fat/index/hg38_bs_index/ \
	--output_dir "${i}"/. \
	"${i}"/bsk.deduplicated.sam

for i in FXH12 FXH15 FXH16 FXH24 FXH26 FXH27 FXH30 FXH32 FXH55 FXH64 FXH71 FXH72 FXH73 FXH74 FXH80 FXH89 FXHN12 FXHN15 FXHN16 FXHN17 FXHN18 FXHN24 FXHN26 FXHN27 FXHN30 FXHN55 FXHN64 FXHN71 FXHN72 FXHN73 FXHN74 FXHN80 FXHN89; do
	sed '/^@/s/ /_/g' "${i}"/merge.assembled.fastq \
		| umi_tools extract \
			--umi-separator=: --extract-method=string --3prime \
			--bc-pattern=NNNNNN \
			-L "${i}"/umi_extract.log \
			-S "${i}"/merge.umi.fastq

	rm "${i}"/merge.assembled.fastq

	/home/ivan/bismark/bismark \
		--genome /home/ivan/fat/index/hg38_bs_index/ \
		--non_directional -q -p 4 \
		"${i}"/merge.umi.fastq \
		--gzip --sam -B "${i}"/bsk

	/home/ivan/bismark/bismark_methylation_extractor \
		--parallel 6 --bedGraph --counts --report --cytosine_report --comprehensive \
		--genome_folder /home/ivan/fat/index/hg38_bs_index/ \
		--output_dir "${i}"/. \
		"${i}"/bsk.sam
done
