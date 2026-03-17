for acc in Cvi-0 Cant-1 Had-6b Blh-1 Ishikawa Yo-0 Ri-0 Are-1 Col-0; do
	mkdir -p methylseq/${acc}_rep{1,2,3}
done

parallel -j 10 --keep-order '
	ln -sf /share/home/wangq/wyf/nfmethylseq/data/{1}.1.fq.gz \
		/share/home/wangq/wyf/methylseq/{2}/R1.fq.gz
	ln -sf /share/home/wangq/wyf/nfmethylseq/data/{1}.2.fq.gz \
		/share/home/wangq/wyf/methylseq/{2}/R2.fq.gz
' ::: NJU_EM_00{{01..05},{07..27}} \
	::: Ishikawa_rep{1..3} Ri-0_rep{1,2} Are-1_rep{1..3} Yo-0_rep{1..3} Cvi-0_rep{1..3} Cant-1_rep{1..3} Had-6b_rep{1..3} Blh-1_rep{1..3} Col-0_rep{1..3}

for acc in Cvi-0 Cant-1 Blh-1 Yo-0 Are-1 Col-0 Ishikawa Had-6b; do
	for rep in rep{1,2,3}; do
		bsub -n 24 -J "${acc}_${rep}" "
			/share/home/wangq/miniconda3/bin/bismark \
				-p 7 -X 1200 -I 50 --gzip -q \
				--genome index/Ath_${acc}_index/ --temp_dir methylseq/${acc}_${rep}/ \
				-B methylseq/${acc}_${rep}/bismark \
				-1 methylseq/${acc}_${rep}/R1.fq.gz \
				-2 methylseq/${acc}_${rep}/R2.fq.gz
			"
	done
done

for acc in Ri-0; do
	for rep in rep{1,2}; do
		bsub -n 24 -J "${acc}_${rep}" "
			/share/home/wangq/miniconda3/bin/bismark \
				-p 7 -X 1200 -I 50 --gzip -q \
				--genome index/Ath_${acc}_index/ --temp_dir methylseq/${acc}_${rep}/ \
				-B methylseq/${acc}_${rep}/bismark \
				-1 methylseq/${acc}_${rep}/R1.fq.gz \
				-2 methylseq/${acc}_${rep}/R2.fq.gz
			"
	done
done