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
	::: NJU-BSP000-{1..4}

rm ./*/merge.discarded.fastq ./*/merge.unassembled.forward.fastq ./*/merge.unassembled.reverse.fastq

for i in NJU-BSP000-{1..4}; do
	sed '/^@/s/ /_/g' "${i}"/merge.assembled.fastq \
		| umi_tools extract \
			--umi-separator=: --extract-method=string --3prime \
			--bc-pattern=NNNNNN \
			-L "${i}"/umi_extract.log \
			-S "${i}"/merge.umi.fastq

	rm "${i}"/merge.assembled.fastq

	/home/ivan/bismark/bismark \
		--genome ../index/hg38_bs/ \
		--non_directional -q -p 4 \
		"${i}"/merge.umi.fastq \
		--gzip --sam -B "${i}"/bsk

	/home/ivan/bismark/deduplicate_bismark \
		--sam --barcode "${i}"/bsk.sam \
		-o "${i}"/bsk

	/home/ivan/bismark/bismark_methylation_extractor \
		--parallel 6 --bedGraph --counts --report --cytosine_report --comprehensive \
		--genome_folder /home/ivan/fat/index/hg38_bs/ \
		-o "${i}"/bsk \
		"${i}"/bsk.deduplicated.sam
done
