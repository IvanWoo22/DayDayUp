mkdir -p 1_training
bmr split training.tsv.gz -c 2000 --mode column --rc 2 -o split | bash
bmr split testing.tsv.gz -c 2000 --mode column --rc 2 -o split | bash

mkdir -p job
find split -type f -name "*[0-9]" |
  sed "s:split/training.::g" |
  sed "s:split/testing.::g" |
  sort | uniq |
  split -l 400 -a 1 -d - job/

mkdir -p 1_training
for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "train-${f}" "
    parallel --no-run-if-empty --line-buffer -k -j 128 '
    echo '\''==> Processing {}'\''
    Rscript univariate_Twist.R Cancer split/training.{} split/testing.{} 1_training/{}.tsv
    ' <${f} "
done

tsv-append -H 1_training/*.tsv >1.training.result.tsv

awk '{h[sprintf("%.1f",$8) "\t" sprintf("%.1f",$9)]++};END{for (i in h){print i "\t" h[i]}}' \
  1.training.result.tsv |
  sort -nrk 1,2
#Train_AUC	Test_AUC	Count
#1.0	0.6	6
#1.0	0.5	2
#1.0	0.4	1
#1.0	0.3	1
#0.9	0.8	11
#0.9	0.7	49
#0.9	0.6	120
#0.9	0.5	133
#0.9	0.4	78
#0.9	0.3	24
#0.9	0.2	12
#0.8	0.9	12
#0.8	0.8	156
#0.8	0.7	934
#0.8	0.6	2728
#0.8	0.5	3302
#0.8	0.4	1769
#0.8	0.3	474
#0.8	0.2	71
#0.8	0.1	4
#0.7	1.0	1
#0.7	0.9	63
#0.7	0.8	1064
#0.7	0.7	12677
#0.7	0.6	59574
#0.7	0.5	69124
#0.7	0.4	25308
#0.7	0.3	4440
#0.7	0.2	542
#0.7	0.1	40
#0.7	0.0	1
#0.6	1.0	3
#0.6	0.9	120
#0.6	0.8	3453
#0.6	0.7	66136
#0.6	0.6	472662
#0.6	0.5	735322
#0.6	0.4	239108
#0.6	0.3	25574
#0.6	0.2	1871
#0.6	0.1	88
#0.6	0.0	7
#0.5	1.0	3
#0.5	0.9	87
#0.5	0.8	2353
#0.5	0.7	50698
#0.5	0.6	467788
#0.5	0.5	1110517
#0.5	0.4	333259
#0.5	0.3	27919
#0.5	0.2	1586
#0.5	0.1	70
#0.4	0.9	9
#0.4	0.8	145
#0.4	0.7	1674
#0.4	0.6	8811
#0.4	0.5	13801
#0.4	0.4	5861
#0.4	0.3	852
#0.4	0.2	74
#0.4	0.1	4
#0.3	0.8	7
#0.3	0.7	29
#0.3	0.6	76
#0.3	0.5	87
#0.3	0.4	70
#0.3	0.3	13
#0.3	0.2	1
#0.2	0.6	2
#0.2	0.5	3
#0.2	0.4	1
#0.2	0.3	1
#0.2	0.2	1
#0.0	0.0	1

train_upper=$(cut -f 8 1.training.result.tsv |
  grep -v "rocauc" | sort -n |
  awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
test_upper=$(cut -f 9 1.training.result.tsv |
  grep -v "testauc" | sort -n |
  awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
echo "Train_AUC" "Test_AUC"
echo "$train_upper" "$test_upper"
#Train_AUC Test_AUC
#0.64887 0.63871

# shellcheck disable=SC2016
grep -v "NA" 1.training.result.tsv |
  keep-header -- awk \
    '$2>1&&$3>1&&$5<0.05' |
  cut -f 1,4-9 |
  keep-header -- awk \
    -va1=0.6 -va2=0.4 -va3=0.6 -va4=0.4 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    >1.training.result.filter.tsv
#rm output.*

bash select_col.sh \
  -f 1-3 training.tsv.gz 1.training.result.filter.tsv \
  >1.training.data.tsv
bash select_col.sh \
  -f 1-3 testing.tsv.gz 1.training.result.filter.tsv \
  >1.testing.data.tsv

cat 1.training.data.tsv >1.data.tsv
grep -v "Sample" 1.testing.data.tsv >>1.data.tsv
bash BS.sh 1.data.tsv 147 1_bootstrap

bsub -n 128 -q amd_milan -J "BS" "
  bash unibootstrap.sh \
    1_bootstrap 1.training.result.filter.tsv \
    Cancer 0.6 0.4 1_bootstrap
"
bash result_stat.sh -b 1_bootstrap/result.tsv
#| #Item | Value |
#| --- | --- |
#| 1_bootstrap/result.tsv | 24737 |
#| count | 24736 |
#| reg_p_max | 0.0436114014 |
#| reg_p_min | 3.017e-06 |
#| rocauc_min | 0.60015 |
#| rocauc_max | 0.96364 |
#| testauc_min | 0.60008 |
#| testauc_max | 0.91818 |
#| bin(100) | count |
#| 100 | 15902 |
#| 95 | 3383 |
#| 90 | 1551 |
#| 85 | 1105 |
#| 80 | 858 |
#| 75 | 685 |
#| 70 | 553 |
#| 65 | 364 |
#| 60 | 208 |
#| 55 | 78 |
#| 50 | 39 |
#| 45 | 9 |
#| 40 | 1 |

BS_PASS=95
tsv-filter -H --ge 8:${BS_PASS} \
  <1_bootstrap/result.tsv \
  >1_training/bs.tsv

rm -fr ./*split ./*job
rm ./output.* 1_training/[0-9]*.tsv 1_bootstrap/*[0-9]

bash select_col.sh -f 1-3 \
  training.tsv.gz 1_training/bs.tsv \
  >1.training.tsv
bash select_col.sh -f 1-3 \
  testing.tsv.gz 1_training/bs.tsv \
  >1.testing.tsv
sed -i "s:Sample:#sample:g" 1.training.tsv
sed -i "s:Sample:#sample:g" 1.testing.tsv

###

mkdir -p 2_training
bmr nextstep 1_training/bs.tsv |
  tsv-uniq >2_training/formula.tsv

bmr split 2_training/formula.tsv \
  -c 20000 --mode row --rr 1 -o split | bash

mkdir -p job
find split -type f -name "*[0-9]" |
  sort |
  split -l 600 -a 3 -d - job/

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "train-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript multivariate_Twist.R Cancer 1.training.tsv 1.testing.tsv {} {}.result.tsv
    ' <${f}
    "
done

tsv-append -H split/*.result.tsv >2.training.result.tsv
rm -fr ./*split ./*job
rm ./output.*

for target in ad mci hc; do
  train_upper=$(cut -f 6 2_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 2_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
done

#0.62727 0.67051
#0.66520 0.67680
#0.61759 0.62731

## reset filter for mci
for target in ad mci hc; do
  keep-header -- awk \
    -va1=0.61759 -va2=0.4 -va3=0.62731 -va4=0.4 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <2_training/${target}.result.tsv \
    >2_training/${target}.result.filter.tsv
done

for target in ad mci hc; do
  bash result_stat.sh 2_training/${target}.result.filter.tsv
done
