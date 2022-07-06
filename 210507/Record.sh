time parallel --keep-order --xapply -j 2 '
perl NJU_seq/vRNA_analysis/score.pl output/{1}/mhv.tsv output/{2}/mhv.tsv output/{3}/mhv.tsv output/{4}/mhv.tsv >output/{5}_scored_v3.tsv
' ::: NJU{6390..6401..4} ::: NJU{6391..6401..4} ::: NJU{6392..6401..4} ::: NJU{6393..6401..4} ::: MHV_particle MHV_infect MHV_uninfect

time parallel --keep-order --xapply -j 2 '
perl score4virus_beta3.pl output/{1}/mhv_rev.tsv output/{2}/mhv_rev.tsv output/{3}/mhv_rev.tsv output/{4}/mhv_rev.tsv >output/{5}_rev_scored_v3.tsv
' ::: NJU{6390..6401..4} ::: NJU{6391..6401..4} ::: NJU{6392..6401..4} ::: NJU{6393..6401..4} ::: MHV_particle MHV_infect MHV_uninfect

time parallel --keep-order --xapply -j 2 '
perl NJU_seq/vRNA_analysis/score.pl output/{1}/MCMV.tsv output/{2}/MCMV.tsv output/{3}/MCMV.tsv output/{4}/MCMV.tsv >output/{5}_scored_v3.tsv
' ::: NJU{6356..6363..4} ::: NJU{6357..6363..4} ::: NJU{6358..6363..4} ::: NJU{6359..6363..4} ::: MCMV_inf MCMV_uni

for PREFIX in FF180702001 FF180702002 FF180702007 FF180702008 FF180702009 FF180702014; do
  mkdir ${PREFIX}
  bsub -n 10 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
trim_galore --illumina --path_to_cutadapt ~/.linuxbrew/bin/cutadapt -j 8 \
--rrbs --stringency 3 --length 35 -e 0.1 --gzip -o ${PREFIX} \
--paired puh/PUH-${PREFIX}*_R1_*.gz puh/PUH-${PREFIX}*_R2_*.gz
"
done

bsub -n 24 -J "bismark_genome_preparation" "
bismark/bismark_genome_preparation --bowtie2 --genomic_composition --parallel 8 hg19
"

for PREFIX in FF180702001 FF180702002 FF180702007 FF180702008 FF180702009 FF180702014; do
  mkdir output/${PREFIX}
done

for PREFIX in FF180702001 FF180702002 FF180702007 FF180702008 FF180702009 FF180702014; do
  mkdir temp/${PREFIX}
done

parallel --keep-order --xapply -j 6 '
bsub -n 24 -o log/{}_hg19meth_alignment.log -J "{}" "
bismark/bismark --genome hg19/ -1 data/{}/R1.fq.gz -2 data/{}/R2.fq.gz \
--parallel 6 -o output/{}/ --temp_dir temp/{}/ --rg_tag --rg_id {} --rg_sample {} \
--unmapped --nucleotide_coverage
"
' ::: FF180702001 FF180702002 FF180702007 FF180702008 FF180702009 FF180702014

parallel --keep-order --xapply -j 6 '
bsub -n 24 -o log/{}_hg19meth_extract.log -J "{}" "
Bismark-0.23.0/bismark_methylation_extractor --paired-end --report --comprehensive --output /Users/yumh/Cleandata/Cleandata/${PREFIX}/ --gzip --bedGraph --ucsc /Users/yumh/Cleandata/Cleandata/${PREFIX}/${PREFIX}_R1_bismark_bt2_pe.bam
"
' ::: FF180702001 FF180702002 FF180702007 FF180702008 FF180702009 FF180702014

parallel --keep-order --xapply -j 5 '
pigz -dc output/{}/R1_bismark_bt2_pe.bismark.cov.gz | awk '\''$5+$6>10{print $1 "\t" $2}'\'' > output/{}/CpG_sites.tsv
' ::: FF180702002 FF180702007 FF180702008 FF180702009 FF180702014

for PREFIX in FF18070200{3..6}; do
  mkdir ${PREFIX}
  bsub -n 10 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
trim_galore --illumina --path_to_cutadapt ~/.linuxbrew/bin/cutadapt -j 8 \
--rrbs --stringency 3 --length 35 -e 0.1 --gzip -o ${PREFIX} \
--paired puh/PUH-${PREFIX}*_R1_*.gz puh/PUH-${PREFIX}*_R2_*.gz
"
done

for PREFIX in FF18070200{3..6}; do
  mkdir output/${PREFIX} temp/${PREFIX}
done

for PREFIX in FF18070200{3..6}; do
  mv data/${PREFIX}/*_val_1.fq.gz data/${PREFIX}/R1.fq.gz
  mv data/${PREFIX}/*_val_2.fq.gz data/${PREFIX}/R2.fq.gz
done

parallel --keep-order --xapply -j 4 '
bsub -n 24 -o log/{}_hg19meth_alignment.log -J "{}" "
bismark/bismark --genome hg19/ -1 data/{}/R1.fq.gz -2 data/{}/R2.fq.gz \
--parallel 6 -o output/{}/ --temp_dir temp/{}/ --rg_tag --rg_id {} --rg_sample {} \
--unmapped --nucleotide_coverage
"
' ::: FF18070200{3..6}

parallel --keep-order --xapply -j 4 '
bsub -n 24 -o log/{}_hg19meth_extract.log -J "{}" "
bismark/bismark_methylation_extractor --paired-end --parallel 8 --report --comprehensive --output output/{}/ --gzip --bedGraph --ucsc output/{}/R1_bismark_bt2_pe.bam
"
' ::: FF18070200{3..6}

parallel --keep-order --xapply -j 4 '
pigz -dc output/{}/R1_bismark_bt2_pe.bismark.cov.gz | awk '\''$5+$6>10{print $1 "\t" $2}'\'' > output/{}/CpG_sites.tsv
' ::: FF18070200{3..6}

for PREFIX in FF180710001 FF180710006 FF180729012 FF180729013 FF180729018; do
  mkdir ${PREFIX}
  bsub -n 10 -o ../log/${PREFIX}_cutadapt.log -J "${PREFIX}" "
trim_galore --illumina --path_to_cutadapt ~/.linuxbrew/bin/cutadapt -j 8 \
--rrbs --stringency 3 --length 35 -e 0.1 --gzip -o ${PREFIX} \
--paired puh/PUH-${PREFIX}*_R1_*.gz puh/PUH-${PREFIX}*_R2_*.gz
"
done

for PREFIX in FF180710001 FF180710006 FF180729012 FF180729013 FF180729018; do
  mkdir output/${PREFIX} temp/${PREFIX}
  mv data/${PREFIX}/*_val_1.fq.gz data/${PREFIX}/R1.fq.gz
  mv data/${PREFIX}/*_val_2.fq.gz data/${PREFIX}/R2.fq.gz
done

parallel --keep-order --xapply -j 5 '
bsub -n 24 -o log/{}_hg19meth_alignment.log -J "{}" "
bismark/bismark --genome hg19/ -1 data/{}/R1.fq.gz -2 data/{}/R2.fq.gz \
--parallel 6 -o output/{}/ --temp_dir temp/{}/ --rg_tag --rg_id {} --rg_sample {} \
--unmapped --nucleotide_coverage
"
' ::: FF180710001 FF180710006 FF180729012 FF180729013 FF180729018

parallel --keep-order --xapply -j 5 '
bsub -n 24 -o log/{}_hg19meth_extract.log -J "{}" "
bismark/bismark_methylation_extractor --paired-end --parallel 8 --report --comprehensive --output output/{}/ --gzip --bedGraph --ucsc output/{}/R1_bismark_bt2_pe.bam
"
' ::: FF180710001 FF180710006 FF180729012 FF180729013 FF180729018

parallel --keep-order --xapply -j 5 '
pigz -dc output/{}/R1_bismark_bt2_pe.bismark.cov.gz | awk '\''$5+$6>10{print $1 "\t" $2}'\'' > output/{}/CpG_sites.tsv
' ::: FF180710001 FF180710006 FF180729012 FF180729013 FF180729018

for PREFIX in X2 X3 X4 X7 X8 X9 X11 X12 X13 X14 X15 X16 X19 X22 X24 X26 X27 X28 X29 X30 X32 X33 X35 X40 X41 X42 X46 X47; do
  echo "${PREFIX}"
  pigz -dc /Users/yumh/Cleandata/Cleandata/${PREFIX}/${PREFIX}_R1_bismark_bt2_pe.bismark.cov.gz | perl test.pl test.tsv
  echo
done

for PREFIX in X2 X3 X4 X7 X8 X9 X11 X12 X13 X14 X15 X16 X19 X22 X24 X26 X27 X28 X29 X30 X32 X33 X35 X40 X41 X42 X46 X47; do
  echo "${PREFIX}"
  parallel --xapply --keep-order -j 4 '
pigz -dc Cleandata/Cleandata/{3}/CHG_context_{3}_R1_bismark_bt2_pe.txt.gz Cleandata/Cleandata/{3}/CHH_context_{3}_R1_bismark_bt2_pe.txt.gz |
awk -v a={1} -v b={2} '\''
BEGIN{sum1=0; sum2=0};
$3==a&&$4==b&&($5=="x"||$5=="h"){sum1++};
$3==a&&$4==b&&($5=="X"||$5=="H"){sum2++};
END{print a "\t" b "\t" sum1 "\t" sum2}
'\''
' ::: cg03046247 cg03046247 cg03046247 cg03046247 cg06457011 cg06457011 cg06457011 cg06457011 cg13695646 cg13695646 cg13695646 cg13695646 cg26508444 cg26508444 cg26508444 cg26508444 cg26508444 cg26508444 cg03356689 cg03356689 cg03356689 cg03356689 cg01899620 cg01899620 cg01899620 cg01899620 cg07109046 cg07109046 cg07109046 cg07109046 cg01798157 cg01798157 cg01798157 cg01798157 cg00161124 cg00161124 cg00161124 cg00161124 cg10535478 cg10535478 cg10535478 cg10535478 cg06580065 cg06580065 cg06580065 cg06580065 cg11893955 cg11893955 cg11893955 cg11893955 cg12141052 cg12141052 cg12141052 cg12141052 cg11770080 cg11770080 cg11770080 cg11770080 cg08813325 cg08813325 cg08813325 cg08813325 cg24662823 cg24662823 cg24662823 cg24662823 cg23146197 cg23146197 cg23146197 cg23146197 cg24127345 cg24127345 cg24127345 cg24127345 cg04481181 cg04481181 cg04481181 cg04481181 cg17450815 cg17450815 cg17450815 cg17450815 ::: 405 415 425 503 294 296 308 360 269 294 306 317 299 323 334 353 366 408 237 295 310 326 274 297 307 379 237 291 299 313 259 300 305 320 255 293 300 309 253 294 306 316 219 253 307 309 272 297 306 379 263 299 309 381 241 297 309 347 246 294 304 348 250 289 310 372 311 315 359 410 299 305 319 405 283 288 309 391 243 300 307 348 ::: ${PREFIX}
  echo
done
