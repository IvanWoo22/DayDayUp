for PREFIX in NJU2029 NJU2030 NJU2031; do
  Bismark-0.23.0/bismark --genome BSref/ \
    -1 Cleandata/${PREFIX}/${PREFIX}_R1.fq.gz \
    -2 Cleandata/${PREFIX}/${PREFIX}_R2.fq.gz \
    -o Cleandata/${PREFIX}/ \
    --temp_dir temp/ --rg_tag --rg_id ${PREFIX} \
    --rg_sample ${PREFIX} --unmapped \
    --nucleotide_coverage
done

for PREFIX in NJU2029 NJU2030 NJU2031; do
  Bismark-0.23.0/bismark_methylation_extractor \
    --paired-end --report --comprehensive \
    --output Cleandata/${PREFIX}/ --gzip \
    --bedGraph Cleandata/${PREFIX}/${PREFIX}_R1_bismark_bt2_pe.bam
done

for PREFIX in NJU2029 NJU2030 NJU2031; do
  parallel --xapply --keep-order -j 4 '
pigz -dc Cleandata/Cleandata/{3}/CHG_context_{3}_R1_bismark_bt2_pe.txt.gz \
Cleandata/Cleandata/{3}/CHH_context_{3}_R1_bismark_bt2_pe.txt.gz |
awk -v a={1} -v b={2} '\''
BEGIN{sum1=0; sum2=0};
$3==a&&$4==b&&($5=="x"||$5=="h"){sum1++};
$3==a&&$4==b&&($5=="X"||$5=="H"){sum2++};
END{print sum1 "\t" sum2}
'\''
' ::: cg03046247 cg03046247 cg03046247 cg03046247 cg06457011 cg06457011 cg06457011 cg06457011 cg13695646 cg13695646 cg13695646 cg13695646 cg26508444 cg26508444 cg26508444 cg26508444 cg26508444 cg26508444 cg03356689 cg03356689 cg03356689 cg03356689 cg01899620 cg01899620 cg01899620 cg01899620 cg07109046 cg07109046 cg07109046 cg07109046 cg01798157 cg01798157 cg01798157 cg01798157 cg00161124 cg00161124 cg00161124 cg00161124 cg10535478 cg10535478 cg10535478 cg10535478 cg06580065 cg06580065 cg06580065 cg06580065 cg11893955 cg11893955 cg11893955 cg11893955 cg12141052 cg12141052 cg12141052 cg12141052 cg11770080 cg11770080 cg11770080 cg11770080 cg08813325 cg08813325 cg08813325 cg08813325 cg24662823 cg24662823 cg24662823 cg24662823 cg23146197 cg23146197 cg23146197 cg23146197 cg24127345 cg24127345 cg24127345 cg24127345 cg04481181 cg04481181 cg04481181 cg04481181 cg17450815 cg17450815 cg17450815 cg17450815 ::: 405 415 425 503 294 296 308 360 269 294 306 317 299 323 334 353 366 408 237 295 310 326 274 297 307 379 237 291 299 313 259 300 305 320 255 293 300 309 253 294 306 316 219 253 307 309 272 297 306 379 263 299 309 381 241 297 309 347 246 294 304 348 250 289 310 372 311 315 359 410 299 305 319 405 283 288 309 391 243 300 307 348 ::: ${PREFIX}
  echo
done

for PREFIX in NJU2029 NJU2030 NJU2031; do
  echo "${PREFIX}"
  pigz -dc Cleandata/${PREFIX}/${PREFIX}_R1_bismark_bt2_pe.bismark.cov.gz |
    perl test.pl test.tsv
  echo
done
