cp training.tsv training.tsv.bak
cp testing.tsv testing.tsv.bak
printf 'loci\tno\n' >site.list
head -1 training.tsv | awk '
    {for(i=3;i<=NF;i++){print $i}}
  ' | awk '
    BEGIN{sum==0000000};{print $1 "\ts" sprintf("%07d", sum); sum++}' \
  >>site.list

head -1 training.tsv >training.head.tmp
awk '
  FNR==NR{dict[$1]=$2; next};{for(i=1;i<=NF;i++){$i=($i in dict) ? dict[$i] : $i}print}
  ' site.list training.head.tmp | sed 's/ /\t/g' >training.tsv
cp training.tsv testing.tsv
awk 'NR>1' training.tsv.bak >>training.tsv
awk 'NR>1' testing.tsv.bak >>testing.tsv
pigz training.tsv testing.tsv
rm training.head.tmp

mkdir -p 1_training
bmr split training.tsv.gz -c 2000 --mode column --rc 2 -o split | bash
bmr split testing.tsv.gz -c 2000 --mode column --rc 2 -o split | bash

mkdir -p job
find split -type f -name "*[0-9]" |
  sed "s:split/training.::g" |
  sed "s:split/testing.::g" |
  sort | uniq |
  split -l 300 -a 1 -d - job/

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
#0.9	0.5	1
#0.8	0.8	5
#0.8	0.7	69
#0.8	0.6	421
#0.8	0.5	420
#0.8	0.4	121
#0.8	0.3	5
#0.8	0.2	1
#0.7	0.9	1
#0.7	0.8	80
#0.7	0.7	5297
#0.7	0.6	40012
#0.7	0.5	43382
#0.7	0.4	10551
#0.7	0.3	649
#0.7	0.2	11
#0.6	0.8	460
#0.6	0.7	38177
#0.6	0.6	378745
#0.6	0.5	629347
#0.6	0.4	168713
#0.6	0.3	10371
#0.6	0.2	129
#0.5	0.9	1
#0.5	0.8	307
#0.5	0.7	28411
#0.5	0.6	359997
#0.5	0.5	937620
#0.5	0.4	248553
#0.5	0.3	14519
#0.5	0.2	162
#0.4	0.8	7
#0.4	0.7	480
#0.4	0.6	4894
#0.4	0.5	9365
#0.4	0.4	3482
#0.4	0.3	240
#0.4	0.2	2
#0.3	0.6	19
#0.3	0.5	45
#0.3	0.4	14
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
#0.63857 0.62782

# shellcheck disable=SC2016
grep -v "NA" 1.training.result.tsv |
  keep-header -- awk \
    '$2>2&&$3>2&&$5<0.05' |
  cut -f 1,4-9 |
  keep-header -- awk \
    -va1=0.55 -va2=0.4 -va3=0.55 -va4=0.4 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    >1.training.result.filter.tsv
#rm output.*

bash select_col.sh \
  -f 1-2 training.tsv.gz 1.training.result.filter.tsv \
  >1.training.data.tsv
bash select_col.sh \
  -f 1-2 testing.tsv.gz 1.training.result.filter.tsv \
  >1.testing.data.tsv

cat 1.training.data.tsv >1.data.tsv
grep -v "Sample" 1.testing.data.tsv >>1.data.tsv
bash BS.sh 1.data.tsv 147 1_bootstrap

bsub -n 128 -q amd_milan -J "BS" "
  bash unibootstrap.sh \
    1_bootstrap 1.training.result.filter.tsv \
    Cancer 0.55 0.4 1_bootstrap
"
bash result_stat.sh -b 1_bootstrap/result.tsv
#| #Item | Value |
#| --- | --- |
#| 1_bootstrap/result.tsv | 21497 |
#| count | 21496 |
#| reg_p_max | 0.0439892474 |
#| reg_p_min | 1.8402e-06 |
#| rocauc_min | 0.57822 |
#| rocauc_max | 0.83796 |
#| testauc_min | 0.55006 |
#| testauc_max | 0.83684 |
#| bin(100) | count |
#| 100 | 14151 |
#| 95 | 4860 |
#| 90 | 1498 |
#| 85 | 616 |
#| 80 | 255 |
#| 75 | 83 |
#| 70 | 28 |
#| 65 | 5 |

BS_PASS=100
tsv-filter -H --ge 8:${BS_PASS} \
  <1_bootstrap/result.tsv \
  >1_training/bs.tsv

rm -fr ./*split ./*job
rm ./output.* 1_training/[0-9]*.tsv 1_bootstrap/*[0-9]

bash select_col.sh -f 1-2 \
  training.tsv.gz 1_training/bs.tsv \
  >1.training.tsv
bash select_col.sh -f 1-2 \
  testing.tsv.gz 1_training/bs.tsv \
  >1.testing.tsv
sed -i "s:Sample:#sample:g" 1.training.tsv
sed -i "s:Sample:#sample:g" 1.testing.tsv

###

mkdir -p 2_training
bmr nextstep 1_training/bs.tsv |
  tsv-uniq >2_training/formula.tsv

bmr split 2_training/formula.tsv \
  -c 12000 --mode row --rr 1 -o split | bash

mkdir -p job
find split -type f -name "*[0-9]" |
  sort |
  split -l 256 -a 1 -d - job/

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
#1.0	0.8	10
#1.0	0.7	42
#1.0	0.6	66
#1.0	0.5	16
#0.9	0.9	84
#0.9	0.8	13182
#0.9	0.7	126275
#0.9	0.6	185067
#0.9	0.5	16571
#0.9	0.4	91
#0.8	0.9	1763
#0.8	0.8	503150
#0.8	0.7	10650119
#0.8	0.6	18315374
#0.8	0.5	966971
#0.8	0.4	2081
#0.8	0.3	1
#0.7	0.9	2801
#0.7	0.8	979053
#0.7	0.7	24280822
#0.7	0.6	41675741
#0.7	0.5	1685374
#0.7	0.4	1938
#0.7	0.3	1
#0.6	0.9	37
#0.6	0.8	9180
#0.6	0.7	229250
#0.6	0.6	451854
#0.6	0.5	21377
#0.6	0.4	25
#0.5	0.7	4
#0.5	0.6	5
#0.0	0.0	1

train_upper=$(cut -f 6 2.training.result.tsv |
  grep -v "rocauc" | sort -n |
  awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
test_upper=$(cut -f 7 2.training.result.tsv |
  grep -v "testauc" | sort -n |
  awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
echo "Train_AUC" "Test_AUC"
echo "$train_upper" "$test_upper"
#Train_AUC Test_AUC
#0.79901 0.71983

rm -fr ./*split ./*job
rm ./output.*

# shellcheck disable=SC2016
grep -v "NA" 2.training.result.tsv |
  keep-header -- awk \
    -F ',|\t' '$5<0.05&&$6<0.05&&$7<0.05' |
  keep-header -- awk \
    -va1=0.7 -va2=0.3 -va3=0.7 -va4=0.3 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    >2.training.result.filter.tsv

bash result_stat.sh 2.training.result.filter.tsv
#| #Item | Value |
#| --- | --- |
#| 2.training.result.filter.tsv | 1371160 |
#| count | 14146 |
#| reg_p_max | 0.0079691848 |
#| reg_p_min | 0.0000000000 |
#| rocauc_min | 0.70009 |
#| rocauc_max | 0.97732 |
#| testauc_min | 0.70008 |
#| testauc_max | 0.92292 |

mkdir -p 2_bootstrap
bash BS.sh 1.data.tsv 147 2_bootstrap
bmr split 2.training.result.filter.tsv \
  -c 275000 --mode row --rr 1 \
  -o split | bash

for f in $(find split -maxdepth 1 -type f -name "*[0-9]" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "bs-${f}" \
    bash multibootstrap.sh 2_bootstrap "${f}" Cancer 0.7 0.3 2_bootstrap
done

rm output.*
rm -fr ./*split

tsv-append -H 2_bootstrap/*.count.tsv \
  >2_bootstrap/result.count
tsv-append -H 2_bootstrap/*.bootstrap.tsv \
  >2.bootstrap.result.tsv

bash result_stat.sh -b 2.bootstrap.result.tsv
#| #Item | Value |
#| --- | --- |
#| 2.bootstrap.result.tsv | 1371160 |
#| count | 14146 |
#| reg_p_max | 0.0079691848 |
#| reg_p_min | 0.0000000000 |
#| rocauc_min | 0.70009 |
#| rocauc_max | 0.97732 |
#| testauc_min | 0.70008 |
#| testauc_max | 0.92292 |
#| bin(BS) | count |
#| 100 | 14845 |
#| 95 | 155689 |
#| 90 | 229600 |
#| 85 | 249501 |
#| 80 | 229522 |
#| 75 | 186664 |
#| 70 | 136518 |
#| 65 | 86946 |
#| 60 | 49067 |
#| 55 | 22527 |
#| 50 | 7817 |
#| 45 | 2053 |
#| 40 | 366 |
#| 35 | 41 |
#| 30 | 3 |

BS_PASS=95
tsv-filter -H --ge 8:${BS_PASS} \
  <2.bootstrap.result.tsv >2_training/bs.tsv

###

mkdir -p 3_training
bmr nextstep 2_training/bs.tsv 1_training/bs.tsv |
  tsv-uniq >3_training/formula.tsv

bmr split 3_training/formula.tsv \
  -c 30000 --mode row --rr 1 -o split | bash

mkdir -p job
find split -type f -name "*[0-9]" |
  sort |
  split -l 400 -a 1 -d - job/

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

tsv-append -H split/*.result.tsv >3.training.result.tsv

awk '{h[sprintf("%.1f",$6) "\t" sprintf("%.1f",$7)]++};END{for (i in h){print i "\t" h[i]}}' \
  3.training.result.tsv |
  sort -nrk 1,2
#Train_AUC	Test_AUC	Count
#1.0	1.0	1062132
#1.0	0.9	7042181
#1.0	0.8	7408598
#1.0	0.7	1291797
#1.0	0.6	299108
#1.0	0.5	86141
#1.0	0.4	25523
#1.0	0.3	9724
#1.0	0.2	5088
#1.0	0.1	1018
#1.0	0.0	7044
#0.9	1.0	2176958
#0.9	0.9	20751362
#0.9	0.8	13450143
#0.9	0.7	954690
#0.9	0.6	125006
#0.9	0.5	23063
#0.9	0.4	5241
#0.9	0.3	1458
#0.9	0.2	483
#0.9	0.1	103
#0.9	0.0	734
#0.8	1.0	300695
#0.8	0.9	1637239
#0.8	0.8	562180
#0.8	0.7	64985
#0.8	0.6	11914
#0.8	0.5	2930
#0.8	0.4	878
#0.8	0.3	346
#0.8	0.2	144
#0.8	0.1	35
#0.8	0.0	166
#0.7	1.0	4132
#0.7	0.9	10679
#0.7	0.8	6392
#0.7	0.7	2014
#0.7	0.6	652
#0.7	0.5	245
#0.7	0.4	107
#0.7	0.3	49
#0.7	0.2	16
#0.7	0.1	13
#0.7	0.0	7
#0.6	1.0	189
#0.6	0.9	432
#0.6	0.8	359
#0.6	0.7	189
#0.6	0.6	123
#0.6	0.5	55
#0.6	0.4	26
#0.6	0.3	20
#0.6	0.2	13
#0.6	0.1	8
#0.6	0.0	3
#0.5	1.0	21
#0.5	0.9	52
#0.5	0.8	58
#0.5	0.7	44
#0.5	0.6	36
#0.5	0.5	25
#0.5	0.4	12
#0.5	0.3	9
#0.5	0.2	5
#0.5	0.1	2
#0.4	1.0	4
#0.4	0.9	4
#0.4	0.8	5
#0.4	0.7	17
#0.4	0.6	16
#0.4	0.5	9
#0.4	0.4	8
#0.4	0.3	3
#0.4	0.2	2
#0.4	0.1	2
#0.4	0.0	2
#0.3	0.9	1
#0.3	0.8	5
#0.3	0.7	11
#0.3	0.6	4
#0.3	0.5	2
#0.3	0.4	3
#0.3	0.3	2
#0.3	0.2	2
#0.2	0.9	1
#0.2	0.7	3
#0.2	0.6	1
#0.2	0.3	1
#0.2	0.2	2
#0.2	0.1	1
#0.1	0.4	1
#0.1	0.3	2
#0.1	0.2	2
#0.1	0.1	2
#0.0	1.0	214
#0.0	0.9	125
#0.0	0.8	253
#0.0	0.7	193
#0.0	0.6	147
#0.0	0.5	415
#0.0	0.4	133
#0.0	0.3	159
#0.0	0.2	212
#0.0	0.1	136
#0.0	0.0	1623

train_upper=$(cut -f 6 3.training.result.tsv |
  grep -v "rocauc" | sort -n |
  awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
test_upper=$(cut -f 7 3.training.result.tsv |
  grep -v "testauc" | sort -n |
  awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
echo "Train_AUC" "Test_AUC"
echo "$train_upper" "$test_upper"
#Train_AUC Test_AUC
#1.00000 0.95833

rm -fr ./*split ./*job
rm ./output.*

# shellcheck disable=SC2016
grep -v "NA" 3.training.result.tsv |
  keep-header -- awk \
    -va1=0.9 -va2=0.1 -va3=0.9 -va4=0.1 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    >3.training.result.filter.tsv

bash result_stat.sh 3.training.result.filter.tsv
#| #Item | Value |
#| --- | --- |
#| 3.training.result.filter.tsv | 8159499 |
#| count | 4996 |
#| reg_p_max | 1.0000000000 |
#| reg_p_min | 0.0000000000 |
#| rocauc_min | 0.06667 |
#| rocauc_max | 1.00000 |
#| testauc_min | 0.06250 |
#| testauc_max | 1.00000 |

mkdir -p 3_bootstrap
bash BS.sh 1.data.tsv 147 3_bootstrap
bmr split 3.training.result.filter.tsv \
  -c 1700000 --mode row --rr 1 \
  -o split | bash

for f in $(find split -maxdepth 1 -type f -name "*[0-9]" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "bs-${f}" \
    bash multibootstrap.sh 3_bootstrap "${f}" Cancer 0.9 0.1 3_bootstrap
done

rm output.*
rm -fr ./*split

tsv-append -H 3_bootstrap/*.count.tsv \
  >3_bootstrap/result.count
tsv-append -H 3_bootstrap/*.bootstrap.tsv \
  >3.bootstrap.result.tsv

bash result_stat.sh -b 3.bootstrap.result.tsv
#| #Item | Value |
#| --- | --- |
#| 3.bootstrap.result.tsv | 8159499 |
#| count | 4996 |
#| reg_p_max | 1.0000000000 |
#| reg_p_min | 0.0000000000 |
#| rocauc_min | 0.06667 |
#| rocauc_max | 1.00000 |
#| testauc_min | 0.06250 |
#| testauc_max | 1.00000 |
#| bin(BS) | count |
#| 100 | 928284 |
#| 95 | 1885235 |
#| 90 | 1434093 |
#| 85 | 1181592 |
#| 80 | 941926 |
#| 75 | 696913 |
#| 70 | 478992 |
#| 65 | 299072 |
#| 60 | 167211 |
#| 55 | 82090 |
#| 50 | 38598 |
#| 45 | 16911 |
#| 40 | 6243 |
#| 35 | 1735 |
#| 30 | 438 |
#| 25 | 124 |
#| 20 | 36 |
#| 15 | 5 |

BS_PASS=100
tsv-filter -H --ge 8:${BS_PASS} \
  <3.bootstrap.result.tsv >3_training/bs.tsv

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "train-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      bash 2.sh  {}
    ' <${f}
    "
done
