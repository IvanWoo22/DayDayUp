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

awk '{h[sprintf("%.1f",$6) "\t" sprintf("%.1f",$7)]++};END{for (i in h){print i "\t" h[i]}}' \
  2.training.result.tsv |
  sort -nrk 1,2
#Train_AUC	Test_AUC	Count
#1.0	1.0	11021
#1.0	0.9	50311
#1.0	0.8	191718
#1.0	0.7	317466
#1.0	0.6	215471
#1.0	0.5	54276
#1.0	0.4	12287
#1.0	0.3	3360
#1.0	0.2	1301
#1.0	0.1	273
#1.0	0.0	596
#0.9	1.0	26542
#0.9	0.9	291100
#0.9	0.8	1932560
#0.9	0.7	4426707
#0.9	0.6	2400128
#0.9	0.5	296022
#0.9	0.4	35112
#0.9	0.3	5966
#0.9	0.2	1670
#0.9	0.1	290
#0.9	0.0	588
#0.8	1.0	32412
#0.8	0.9	623546
#0.8	0.8	7995512
#0.8	0.7	37301345
#0.8	0.6	16385699
#0.8	0.5	703829
#0.8	0.4	49556
#0.8	0.3	6668
#0.8	0.2	1648
#0.8	0.1	267
#0.8	0.0	445
#0.7	1.0	8384
#0.7	0.9	198045
#0.7	0.8	6986846
#0.7	0.7	71987893
#0.7	0.6	29917342
#0.7	0.5	313252
#0.7	0.4	16031
#0.7	0.3	2407
#0.7	0.2	608
#0.7	0.1	126
#0.7	0.0	161
#0.6	1.0	857
#0.6	0.9	7940
#0.6	0.8	185928
#0.6	0.7	1976299
#0.6	0.6	941721
#0.6	0.5	16074
#0.6	0.4	2053
#0.6	0.3	614
#0.6	0.2	228
#0.6	0.1	54
#0.6	0.0	45
#0.5	1.0	64
#0.5	0.9	291
#0.5	0.8	1063
#0.5	0.7	1944
#0.5	0.6	1618
#0.5	0.5	574
#0.5	0.4	268
#0.5	0.3	144
#0.5	0.2	73
#0.5	0.1	17
#0.5	0.0	14
#0.4	1.0	4
#0.4	0.9	18
#0.4	0.8	43
#0.4	0.7	60
#0.4	0.6	66
#0.4	0.5	34
#0.4	0.4	13
#0.4	0.3	11
#0.4	0.2	9
#0.4	0.1	4
#0.3	0.9	1
#0.3	0.8	5
#0.3	0.7	5
#0.3	0.6	5
#0.3	0.5	5
#0.3	0.4	1
#0.3	0.3	4
#0.3	0.2	2
#0.3	0.1	3
#0.2	1.0	1
#0.2	0.6	1
#0.2	0.4	1
#0.2	0.0	1
#0.0	1.0	18
#0.0	0.9	18
#0.0	0.8	44
#0.0	0.7	52
#0.0	0.6	59
#0.0	0.5	102
#0.0	0.4	82
#0.0	0.3	52
#0.0	0.2	47
#0.0	0.1	24
#0.0	0.0	506

train_upper=$(cut -f 6 2.training.result.tsv |
  grep -v "rocauc" | sort -n |
  awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
test_upper=$(cut -f 7 2.training.result.tsv |
  grep -v "testauc" | sort -n |
  awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
echo "Train_AUC" "Test_AUC"
echo "$train_upper" "$test_upper"
#Train_AUC Test_AUC
#0.85526 0.77515

rm -fr ./*split ./*job
rm ./output.*

# shellcheck disable=SC2016
grep -v "NA" 2.training.result.tsv |
  keep-header -- awk \
    -va1=0.8 -va2=0.2 -va3=0.8 -va4=0.2 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    >2.training.result.filter.tsv

bash result_stat.sh 2.training.result.filter.tsv
#| #Item | Value |
#| --- | --- |
#| 2.training.result.filter.tsv | 5467398 |
#| count | 19285 |
#| reg_p_max | 1.0000000000 |
#| reg_p_min | 0.0000000000 |
#| rocauc_min | 0.80007 |
#| rocauc_max | 1.00000 |
#| testauc_min | 0.75019 |
#| testauc_max | 1.00000 |

mkdir -p 2_bootstrap
bash BS.sh 1.data.tsv 147 2_bootstrap
bmr split 2.training.result.filter.tsv \
  -c 1000000 --mode row --rr 1 \
  -o split | bash

for f in $(find split -maxdepth 1 -type f -name "*[0-9]" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "bs-${f}" \
    bash multibootstrap.sh 2_bootstrap "${f}" Cancer 0.78 0.22 2_bootstrap
done

rm output.*
rm -fr ./*split

tsv-append -H 2_bootstrap/*.count.tsv \
  >2_bootstrap/result.count
tsv-append -H 2_bootstrap/*.bootstrap.tsv \
  >2.bootstrap.result.tsv

bash result_stat.sh -b 2.bootstrap.result.tsv

BS_PASS=95
tsv-filter -H --ge 8:${BS_PASS} \
  <2.bootstrap.result.tsv >2_training/bs.tsv

###

mkdir -p 3_training
bmr nextstep 2_training/bs.tsv 1_training/bs.tsv |
  tsv-uniq >3_training/formula.tsv

bmr split 3_training/formula.tsv \
  -c 300000 --mode row --rr 1 -o split | bash

mkdir -p job
find split -type f -name "*[0-9]" |
  sort |
  split -l 1450 -a 1 -d - job/

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

## hc might change 15000 to 30000
for target in ad mci hc; do
  bmr split 3_training/"${target}".formula.tsv \
    -c 15000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad mci hc; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 200 -a 3 -d - "${target}"_job/
done

for target in ad mci hc; do
  for f in $(find "${target}"_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo "${f}"
    bsub -n 24 -J "training-${f}" \
      "
        parallel --no-run-if-empty --line-buffer -k -j 24 '
        echo '\''==> Processing {}'\''
        Rscript multivariate.R ${target} 1_${target}.training.tsv 1_${target}.testing.tsv {} {}.result.tsv
    ' < ${f}
    "
  done
done
