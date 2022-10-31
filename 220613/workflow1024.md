mkdir -p 1_training

bmr split training.tsv.gz -c 2000 --mode column --rc 3 -o split | bash
bmr split testing.tsv.gz -c 2000 --mode column --rc 3 -o split | bash

mkdir -p job

find split -type f -name "*[0-9]" |
  sed "s:split/training.::g" |
  sed "s:split/testing.::g" |
  sort | uniq |
  split -l 120 -a 3 -d - job/

mkdir -p 1_ad_training

for target in ad; do
  for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo "${target}-${f}"
    bsub -n 24 -J "training-${f}" "
    parallel --no-run-if-empty --line-buffer -k -j 24 '
      echo '\''==> Processing {}'\''
      Rscript univariate.R ${target} split/training.{} split/testing.{} 1_${target}_training/{}.tsv
    ' <${f} "
  done
done

tsv-append -H 1_ad_training/*.tsv >1_training/ad.result.tsv

for target in ad; do
  cut -f 1,4-9 1_training/${target}.result.tsv \
  >1_training/${target}.result
  bash result_stat.sh 1_training/${target}.result
  rm 1_training/${target}.result
done

| #Item | Value |
| --- | --- |
| 1_training/ad.result | 837611 |
| count | 837610 |
| rocauc_min | 0.36238 |
| rocauc_median | 0.52229 |
| rocauc_max | 0.67772 |
| testauc_min | 0.27941 |
| testauc_median | 0.50619 |
| testauc_max | 0.73188 |

for target in ad; do
  train_upper=$(cut -f 8 1_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 9 1_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    '$2>0.01&&$3>0.01&&$5<0.05' \
    <1_training/${target}.result.tsv |
    cut -f 1,4-9 \
      >1_training/${target}.result.filter.tmp
  keep-header -- awk \
    -va1=${train_upper} -va2=0.45 -va3=${train_upper} -va4=0.45 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <1_training/${target}.result.filter.tmp \
    >1_training/${target}.result.filter.tsv
done

for target in ad; do
  bash result_stat.sh 1_training/${target}.result.filter.tmp
done

| #Item | Value |
| --- | --- |
| 1_training/ad.result.filter.tmp | 23954 |
| count | 23953 |
| rocauc_min | 0.52111 |
| rocauc_median | 0.57788 |
| rocauc_max | 0.67336 |
| testauc_min | 0.32628 |
| testauc_median | 0.59512 |
| testauc_max | 0.73176 |

for target in ad; do
  bash result_stat.sh 1_training/${target}.result.filter.tsv
done

| #Item | Value |
| --- | --- |
| 1_training/ad.result.filter.tsv | 10889 |
| count | 10888 |
| rocauc_min | 0.57051 |
| rocauc_median | 0.58463 |
| rocauc_max | 0.67336 |
| testauc_min | 0.57060 |
| testauc_median | 0.61832 |
| testauc_max | 0.73176 |

for target in ad; do
  bash select_col.sh \
    -f 1-3 training.tsv.gz 1_training/${target}.result.filter.tsv \
    >1_training/${target}.data.tsv
done

for target in ad; do
  bash BS.sh 1_training/${target}.data.tsv 422 1_${target}_bootstrap
done

for target in ad; do
  bsub -n 24 -q mpi -J "bs_${target}" \
    bash unibootstrap.sh \
    1_${target}_bootstrap \
    1_training/${target}.result.filter.tsv \
    ${target} 0.55 0.45 \
    1_${target}_bootstrap
done

for target in ad; do
  mv 1_${target}_bootstrap/${target}.result.filter.tsv.bootstrap.tsv \
    1_${target}_bootstrap/result.tsv
done

for target in ad; do
  bash result_stat.sh -b 1_${target}_bootstrap/result.tsv
done

BS_PASS=90
for target in ad; do
    tsv-filter -H --ge 8:${BS_PASS} \
      <1_"${target}"_bootstrap/result.tsv \
      >1_training/"${target}".bs.tsv
done

rm -fr ./*split ./*job
rm ./output.*

for target in ad; do
  bash select_col.sh -f 1-3 \
    training.tsv.gz 1_training/"${target}".bs.tsv \
    >1_"${target}".training.tsv
  bash select_col.sh -f 1-3 \
    testing.tsv.gz 1_training/"${target}".bs.tsv \
    >1_"${target}".testing.tsv
  sed -i "s:Row.names:#sample:g" 1_"${target}".training.tsv
  sed -i "s:Row.names:#sample:g" 1_"${target}".testing.tsv
done

======

mkdir -p 2_training

for target in ad; do
  sed -i "s:ProbeID:#marker:g" 1_training/"${target}".bs.tsv
  bmr nextstep 1_training/"${target}".bs.tsv |
    tsv-uniq >2_training/"${target}".formula.tsv
done

for target in ad; do
  bmr split 2_training/"${target}".formula.tsv \
    -c 8000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 120 -a 3 -d - "${target}"_job/
done

for target in ad; do
  for f in $(find "${target}"_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo "${f}"
    bsub -n 24 -q mpi -J "training-${f}" \
      "
        parallel --no-run-if-empty --line-buffer -k -j 24 '
        echo '\''==> Processing {}'\''
        Rscript multivariate.R ${target} 1_${target}.training.tsv 1_${target}.testing.tsv {} {}.result.tsv
    ' < ${f}
    "
  done
done


tsv-append -H ad_split/*.result.tsv >2_training/ad.result.tsv

rm -fr ./*_split ./*job
rm ./output.*

for target in ad; do
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
for target in ad; do
  keep-header -- awk \
    -va1=0.63258 -va2=0.4 -va3=0.63258 -va4=0.4 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <2_training/${target}.result.tsv \
    >2_training/${target}.result.filter.tsv
done

for target in ad; do
  bash result_stat.sh 2_training/${target}.result.filter.tsv
done

for target in ad; do
  mkdir -p 2_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 422 2_${target}_bootstrap
  bmr split 2_training/${target}.result.filter.tsv \
    -c 10000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo "${f}"
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh 2_${target}_bootstrap "${f}" ${target} 0.6 0.4 2_${target}_bootstrap
  done
done
## ......
## ......
## ......

rm output.*
rm -fr ./*_split

for target in ad; do
  tsv-append -H 2_${target}_bootstrap/*.count.tsv \
    >2_${target}_bootstrap/${target}.result.count
  tsv-append -H 2_${target}_bootstrap/*.bootstrap.tsv \
    >2_${target}_bootstrap/result.tsv
done

for target in ad; do
  bash result_stat.sh -b 2_${target}_bootstrap/result.tsv
done

BS_PASS=95
for target in ad; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <2_"${target}"_bootstrap/result.tsv >2_training/"${target}".bs.tsv
done

----

mkdir -p 3_training

for target in ad; do
  sed -i "s:ProbeID:#marker:g" 2_training/"${target}".bs.tsv
  bmr nextstep 2_training/"${target}".bs.tsv 1_training/"${target}".bs.tsv |
    tsv-uniq >3_training/"${target}".formula.tsv
done

## hc might change 15000 to 30000
for target in ad; do
  bmr split 3_training/"${target}".formula.tsv \
    -c 15000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 200 -a 3 -d - "${target}"_job/
done

for target in ad; do
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

tsv-append -H ad_split/*.result.tsv >3_training/ad.result.tsv

rm -fr ./*_split ./*_job

for target in ad; do
  train_upper=$(cut -f 6 3_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 3_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.35 -va3="${train_upper}" -va4=0.35 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <3_training/${target}.result.tsv \
    >3_training/${target}.result.filter.tsv
done

#0.70247 0.70004

for target in ad; do
  bash result_stat.sh 3_training/${target}.result.filter.tsv
done

for target in ad; do
  mkdir -p 3_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 422 3_${target}_bootstrap
  bmr split 3_training/${target}.result.filter.tsv \
    -c 6000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh 3_${target}_bootstrap ${f} ${target} 0.7 0.35 3_${target}_bootstrap
  done
done

rm ./output.*
rm -fr ./*_split

for target in ad; do
  tsv-append -H 3_${target}_bootstrap/*.count.tsv \
    >3_${target}_bootstrap/${target}.result.count
  tsv-append -H 3_${target}_bootstrap/*.bootstrap.tsv \
    >3_${target}_bootstrap/result.tsv
done

for target in ad; do
  bash result_stat.sh -b 3_${target}_bootstrap/result.tsv
done

BS_PASS=85
for target in ad; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <3_"${target}"_bootstrap/result.tsv >3_training/"${target}".bs.tsv
done

---

mkdir -p 4_training

for target in ad; do
  sed -i "s:ProbeID:#marker:g" 3_training/"${target}".bs.tsv
  bmr nextstep 3_training/"${target}".bs.tsv 1_training/"${target}".bs.tsv |
    tsv-uniq >4_training/"${target}".formula.tsv
done

for target in ad; do
  bmr split 4_training/"${target}".formula.tsv \
    -c 20000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 200 -a 3 -d - "${target}"_job/
done

for target in ad; do
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

tsv-append -H ad_split/*.result.tsv >4_training/ad.result.tsv

rm -fr ./*_split ./*_job

for target in ad; do
  train_upper=$(cut -f 6 4_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 4_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.3 -va3="${train_upper}" -va4=0.3 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <4_training/${target}.result.tsv \
    >4_training/${target}.result.filter.tsv
done

#0.74351 0.73200

for target in ad; do
  bash result_stat.sh 4_training/${target}.result.filter.tsv
done

for target in ad; do
  mkdir -p 4_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 422 4_${target}_bootstrap
  bmr split 4_training/${target}.result.filter.tsv \
    -c 6000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh 4_${target}_bootstrap ${f} ${target} 0.74 0.3 4_${target}_bootstrap
  done
done

rm ./output.*
rm -fr ./*_split

for target in ad; do
  tsv-append -H 4_${target}_bootstrap/*.count.tsv \
    >4_${target}_bootstrap/${target}.result.count
  tsv-append -H 4_${target}_bootstrap/*.bootstrap.tsv \
    >4_${target}_bootstrap/result.tsv
done

for target in ad; do
  bash result_stat.sh -b 4_${target}_bootstrap/result.tsv
done

BS_PASS=75
for target in ad; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <4_"${target}"_bootstrap/result.tsv >4_training/"${target}".bs.tsv
done

---

mkdir -p 5_training

for target in ad; do
  sed -i "s:ProbeID:#marker:g" 4_training/"${target}".bs.tsv
  bmr nextstep 4_training/"${target}".bs.tsv 1_training/"${target}".bs.tsv |
    tsv-uniq >5_training/"${target}".formula.tsv
done

for target in ad; do
  bmr split 5_training/"${target}".formula.tsv \
    -c 30000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 200 -a 3 -d - "${target}"_job/
done

for target in ad; do
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

tsv-append -H ad_split/*.result.tsv >5_training/ad.result.tsv

rm -fr ./*_split ./*_job

for target in ad; do
  train_upper=$(cut -f 6 5_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 5_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.25 -va3="${train_upper}" -va4=0.25 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <5_training/${target}.result.tsv \
    >5_training/${target}.result.filter.tsv
done

#0.76622 0.76073

for target in ad; do
  bash result_stat.sh 5_training/${target}.result.filter.tsv
done

| #Item | Value |
| --- | --- |
| 5_training/ad.result.filter.tsv | 356077 |
| count | 10514 |
| rocauc_min | 0.76625 |
| rocauc_median | 0.76917 |
| rocauc_max | 0.81337 |
| testauc_min | 0.76625 |
| testauc_median | 0.76986 |
| testauc_max | 0.82382 |

for target in ad; do
  mkdir -p 5_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 422 5_${target}_bootstrap
  bmr split 5_training/${target}.result.filter.tsv \
    -c 8000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh 5_${target}_bootstrap ${f} ${target} 0.76 0.25 5_${target}_bootstrap
  done
done

rm ./output.*
rm -fr ./*_split

for target in ad; do
  tsv-append -H 5_${target}_bootstrap/*.count.tsv \
    >5_${target}_bootstrap/${target}.result.count
  tsv-append -H 5_${target}_bootstrap/*.bootstrap.tsv \
    >5_${target}_bootstrap/result.tsv
done

for target in ad; do
  bash result_stat.sh -b 5_${target}_bootstrap/result.tsv
done

BS_PASS=85
for target in ad; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <5_"${target}"_bootstrap/result.tsv >5_training/"${target}".bs.tsv
done

---

mkdir -p 6_training

for target in ad; do
  sed -i "s:ProbeID:#marker:g" 5_training/"${target}".bs.tsv
  bmr nextstep 5_training/"${target}".bs.tsv 1_training/"${target}".bs.tsv |
    tsv-uniq >6_training/"${target}".formula.tsv
done

for target in ad; do
  bmr split 6_training/"${target}".formula.tsv \
    -c 20000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 120 -a 3 -d - "${target}"_job/
done

for target in ad; do
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

tsv-append -H ad_split/*.result.tsv >6_training/ad.result.tsv

rm -fr ./*_split ./*_job

for target in ad; do
  train_upper=$(cut -f 6 6_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 6_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.2 -va3="${train_upper}" -va4=0.2 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <6_training/${target}.result.tsv \
    >6_training/${target}.result.filter.tsv
done

#0.79285 0.78308

for target in ad; do
  bash result_stat.sh 6_training/${target}.result.filter.tsv
done

| #Item | Value |
| --- | --- |
| 6_training/ad.result.filter.tsv | 28224 |
| count | 3302 |
| rocauc_min | 0.79288 |
| rocauc_median | 0.79535 |
| rocauc_max | 0.82177 |
| testauc_min | 0.79293 |
| testauc_median | 0.79594 |
| testauc_max | 0.83920 |

for target in ad; do
  mkdir -p 6_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 422 6_${target}_bootstrap
  bmr split 6_training/${target}.result.filter.tsv \
    -c 6000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh 6_${target}_bootstrap ${f} ${target} 0.79 0.2 6_${target}_bootstrap
  done
done

rm ./output.*
rm -fr ./*_split

for target in ad; do
  tsv-append -H 6_${target}_bootstrap/*.count.tsv \
    >6_${target}_bootstrap/${target}.result.count
  tsv-append -H 6_${target}_bootstrap/*.bootstrap.tsv \
    >6_${target}_bootstrap/result.tsv
done

for target in ad; do
  bash result_stat.sh -b 6_${target}_bootstrap/result.tsv
done

BS_PASS=70
for target in ad; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <6_"${target}"_bootstrap/result.tsv >6_training/"${target}".bs.tsv
done

---

mkdir -p 7_training

for target in ad; do
  sed -i "s:ProbeID:#marker:g" 6_training/"${target}".bs.tsv
  bmr nextstep 6_training/"${target}".bs.tsv 1_training/"${target}".bs.tsv |
    tsv-uniq >7_training/"${target}".formula.tsv
done

for target in ad; do
  bmr split 7_training/"${target}".formula.tsv \
    -c 20000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 120 -a 3 -d - "${target}"_job/
done

for target in ad; do
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

tsv-append -H ad_split/*.result.tsv >7_training/ad.result.tsv

rm -fr ./*_split ./*_job

for target in ad; do
  train_upper=$(cut -f 6 7_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 7_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.15 -va3="${train_upper}" -va4=0.15 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <7_training/${target}.result.tsv \
    >7_training/${target}.result.filter.tsv
done

#0.80647 0.80651

for target in ad; do
  bash result_stat.sh 7_training/${target}.result.filter.tsv
done

| #Item | Value |
| --- | --- |
| 7_training/ad.result.filter.tsv | 581937 |
| count | 10512 |
| rocauc_min | 0.80650 |
| rocauc_median | 0.80842 |
| rocauc_max | 0.83333 |
| testauc_min | 0.80651 |
| testauc_median | 0.80988 |
| testauc_max | 0.85458 |

for target in ad; do
  mkdir -p 7_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 422 7_${target}_bootstrap
  bmr split 7_training/${target}.result.filter.tsv \
    -c 6000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh 7_${target}_bootstrap ${f} ${target} 0.80 0.15 7_${target}_bootstrap
  done
done

rm ./output.*
rm -fr ./*_split

for target in ad; do
  tsv-append -H 7_${target}_bootstrap/*.count.tsv \
    >7_${target}_bootstrap/${target}.result.count
  tsv-append -H 7_${target}_bootstrap/*.bootstrap.tsv \
    >7_${target}_bootstrap/result.tsv
done

for target in ad; do
  bash result_stat.sh -b 7_${target}_bootstrap/result.tsv
done

BS_PASS=95
for target in ad; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <7_"${target}"_bootstrap/result.tsv >7_training/"${target}".result1.tsv
done