```shell
mkdir -p 1_training split job

bmr split training.tsv.gz -c 1000 --mode column --rc 3 -o split | bash

find split -type f -name "*[0-9]" |
  sed "s:split/training.::g" | sort |
  split -l 120 -a 3 -d - job/

bmr split testing.tsv.gz -c 1000 --mode column --rc 3 -o split | bash

mkdir -p 1_ad_training 1_mci_training 1_hc_training
for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 24 -J "training-${f}" "
    parallel --no-run-if-empty --line-buffer -k -j 24 '
      echo '\''==> Processing {}'\''
      Rscript univariate.R ad split/training.{} split/testing.{} 1_ad_training/{}.tsv
    ' <${f} "
done
for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 24 -J "training-${f}" "
    parallel --no-run-if-empty --line-buffer -k -j 24 '
      echo '\''==> Processing {}'\''
      Rscript univariate.R mci split/training.{} split/testing.{} 1_mci_training/{}.tsv
    ' <${f} "
done
for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 24 -J "training-${f}" "
    parallel --no-run-if-empty --line-buffer -k -j 24 '
      echo '\''==> Processing {}'\''
      Rscript univariate.R hc split/training.{} split/testing.{} 1_hc_training/{}.tsv
    ' <${f} "
done
```

```shell
for target in ad mci hc; do
  tsv-append -H 1_${target}_training/*.tsv >1_training/${target}.result.tsv
done

for target in ad mci hc; do
  train_upper=$(cut -f 6 1_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 1_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.45 -va3="${test_upper}" -va4=0.45 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <1_training/${target}.result.tsv \
    >1_training/${target}.result.filter.tsv
done

for target in ad mci hc; do
  bash result_stat.sh 1_training/${target}.result.tsv
done

for target in ad mci hc; do
  bash result_stat.sh 1_training/${target}.result.filter.tsv
done
```

| #Item                    | Value        |
|--------------------------|--------------|
| 1_training/ad.result.tsv | 837611       |
| count                    | 837610       |
| reg_p_median             | 0.4319586861 |
| reg_p_min                | 4e-10        |
| rocauc_min               | 0.37984      |
| rocauc_max               | 0.71212      |
| testauc_min              | 0.32815      |
| testauc_max              | 0.73506      |

| #Item                     | Value        |
|---------------------------|--------------|
| 1_training/mci.result.tsv | 837611       |
| count                     | 837610       |
| reg_p_median              | 0.5259647186 |
| reg_p_min                 | 3.7237e-06   |
| rocauc_min                | 0.35929      |
| rocauc_max                | 0.70320      |
| testauc_min               | 0.30665      |
| testauc_max               | 0.68835      |

| #Item                    | Value        |
|--------------------------|--------------|
| 1_training/hc.result.tsv | 837611       |
| count                    | 837610       |
| reg_p_median             | 0.4178692309 |
| reg_p_min                | 3.26e-08     |
| rocauc_min               | 0.38884      |
| rocauc_max               | 0.67869      |
| testauc_min              | 0.34686      |
| testauc_max              | 0.68239      |

| #Item                           | Value         |
|---------------------------------|---------------|
| 1_training/ad.result.filter.tsv | 8833          |
| count                           | 8832          |
| reg_p_median                    | 0.00466343225 |
| reg_p_min                       | 3.4e-09       |
| rocauc_min                      | 0.39476       |
| rocauc_max                      | 0.71064       |
| testauc_min                     | 0.34619       |
| testauc_max                     | 0.73506       |

| #Item                            | Value        |
|----------------------------------|--------------|
| 1_training/mci.result.filter.tsv | 2903         |
| count                            | 2902         |
| reg_p_median                     | 0.0442313075 |
| reg_p_min                        | 6.327e-06    |
| rocauc_min                       | 0.36479      |
| rocauc_max                       | 0.70320      |
| testauc_min                      | 0.37427      |
| testauc_max                      | 0.68788      |

| #Item                           | Value        |
|---------------------------------|--------------|
| 1_training/hc.result.filter.tsv | 10774        |
| count                           | 10773        |
| reg_p_median                    | 0.0027098496 |
| reg_p_min                       | 5.64e-08     |
| rocauc_min                      | 0.39883      |
| rocauc_max                      | 0.67750      |
| testauc_min                     | 0.38077      |
| testauc_max                     | 0.68239      |

```shell
for target in ad mci hc; do
  bash select_col.sh \
    -f 1-3 training.tsv.gz 1_training/${target}.result.filter.tsv \
    >1_training/${target}.data.tsv
done
```

```shell
for target in ad mci hc; do
  bash BS.sh 1_training/${target}.data.tsv 365 1_${target}_bootstrap
done

for target in ad mci hc; do
  bsub -n 24 -J "bs_${target}" \
    bash unibootstrap.sh 1_${target}_bootstrap \
    1_training/${target}.result.filter.tsv \
    ${target} 0.55 0.45 \
    1_${target}_bootstrap
done

for target in ad mci hc; do
  mv 1_${target}_bootstrap/${target}.result.filter.tsv.bootstrap.tsv \
    1_${target}_bootstrap/result.tsv
done

for target in ad mci hc; do
  bash result_stat.sh -b 1_${target}_bootstrap/result.tsv
done
```

| #Item                     | Value        |
|---------------------------|--------------|
| 1_ad_bootstrap/result.tsv | 8774         |
| count                     | 8773         |
| reg_p_median              | 0.0046161486 |
| reg_p_min                 | 3.4e-09      |
| rocauc_min                | 0.39476      |
| rocauc_max                | 0.71064      |
| testauc_min               | 0.34619      |
| testauc_max               | 0.73506      |
| bin(BS)                   | count        |
| 100                       | 185          |
| 95                        | 2017         |
| 90                        | 2961         |
| 85                        | 2519         |
| 80                        | 849          |
| 75                        | 78           |
| 70                        | 16           |
| 65                        | 31           |
| 60                        | 34           |
| 55                        | 44           |
| 50                        | 31           |
| 45                        | 8            |

| #Item                      | Value         |
|----------------------------|---------------|
| 1_mci_bootstrap/result.tsv | 2813          |
| count                      | 2812          |
| reg_p_median               | 0.04211826805 |
| reg_p_min                  | 6.327e-06     |
| rocauc_min                 | 0.36479       |
| rocauc_max                 | 0.70320       |
| testauc_min                | 0.37427       |
| testauc_max                | 0.68788       |
| bin(BS)                    | count         |
| 100                        | 32            |
| 95                         | 395           |
| 90                         | 544           |
| 85                         | 679           |
| 80                         | 704           |
| 75                         | 303           |
| 70                         | 88            |
| 65                         | 14            |
| 60                         | 18            |
| 55                         | 7             |
| 50                         | 20            |
| 45                         | 8             |

| #Item                     | Value        |
|---------------------------|--------------|
| 1_hc_bootstrap/result.tsv | 10774        |
| count                     | 10773        |
| reg_p_median              | 0.0027098496 |
| reg_p_min                 | 5.64e-08     |
| rocauc_min                | 0.39883      |
| rocauc_max                | 0.67750      |
| testauc_min               | 0.38077      |
| testauc_max               | 0.68239      |
| bin(BS)                   | count        |
| 100                       | 550          |
| 95                        | 5252         |
| 90                        | 3446         |
| 85                        | 1244         |
| 80                        | 200          |
| 75                        | 18           |
| 70                        | 12           |
| 65                        | 9            |
| 60                        | 12           |
| 55                        | 13           |
| 50                        | 10           |
| 45                        | 5            |
| 40                        | 2            |

```shell
BS_PASS=95
for target in ad mci hc; do
    tsv-filter -H --ge 8:${BS_PASS} \
      <1_"${target}"_bootstrap/result.tsv \
      >1_training/"${target}".bs.tsv
done

for target in ad mci hc; do
  bash select_col.sh -f 1-3 \
    training.tsv.gz 1_training/"${target}".bs.tsv \
    >1_"${target}".training.tsv
  bash select_col.sh -f 1-3 \
    testing.tsv.gz 1_training/"${target}".bs.tsv \
    >1_"${target}".testing.tsv
  sed -i "s:Row.names:#sample:g" 1_"${target}".training.tsv
  sed -i "s:Row.names:#sample:g" 1_"${target}".testing.tsv
done
```

```shell
mkdir -p 2_training

for target in ad mci hc; do
  sed -i "s:ProbeID:#marker:g" 1_training/"${target}".bs.tsv
  bmr nextstep 1_training/"${target}".bs.tsv |
    tsv-uniq >2_training/"${target}".formula.tsv
done

for target in ad mci hc; do
  bmr split 2_training/"${target}".formula.tsv \
    -c 6000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad mci hc; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 120 -a 3 -d - "${target}"_job/
done

for target in ad mci hc; do
  for f in $(find "${target}"_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo "${f}"
    bsub -n 24 -J "training-${f}" \
      "
        parallel --no-run-if-empty --line-buffer -k -j 24 '
        echo '\''==> Processing {}'\''
        Rscript multivariate.R ${target} 1_"${target}".training.tsv 1_"${target}".testing.tsv {} {}.result.tsv
    ' < ${f}
    "
  done
done

parallel --xapply -j 3 '
  tsv-append -H {}_split/*.result.tsv >2_training/{}.result.tsv
' ::: ad mci hc
rm -fr *_split *_job
rm output.*

for target in ad mci hc; do
  train_upper=$(cut -f 6 2_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 2_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.4 -va3="${test_upper}" -va4=0.4 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <2_training/${target}.result.tsv \
    >2_training/${target}.result.filter.tsv
done

#0.66623 0.63604
#0.69595 0.64931
#0.63966 0.60898

## reset filter for mci
  keep-header -- awk \
    -va1=0.66623 -va2=0.4 -va3=0.63604 -va4=0.4 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <2_training/mci.result.tsv \
    >2_training/mci.result.filter.tsv

for target in ad mci hc; do
  bash result_stat.sh 2_training/${target}.result.filter.tsv
done
```

| #Item                           | Value        |
|---------------------------------|--------------|
| 2_training/ad.result.filter.tsv | 13886        |
| count                           | 2014         |
| reg_p_median                    | 9.01e-07     |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.66627      |
| rocauc_max                      | 0.73982      |
| testauc_min                     | 0.63608      |
| testauc_max                     | 0.72083      |

| #Item                            | Value       |
|----------------------------------|-------------|
| 2_training/mci.result.filter.tsv | 6740        |
| count                            | 419         |
| reg_p_median                     | 7.04197e-05 |
| reg_p_min                        | 1.53e-08    |
| rocauc_min                       | 0.36560     |
| rocauc_max                       | 0.73871     |
| testauc_min                      | 0.39311     |
| testauc_max                      | 0.74421     |

| #Item                           | Value     |
|---------------------------------|-----------|
| 2_training/hc.result.filter.tsv | 135900    |
| count                           | 5801      |
| reg_p_median                    | 2.115e-06 |
| reg_p_min                       | 1e-10     |
| rocauc_min                      | 0.63970   |
| rocauc_max                      | 0.70250   |
| testauc_min                     | 0.60901   |
| testauc_max                     | 0.69717   |

```shell
for target in ad mci hc; do
  mkdir -p 2_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 365 2_${target}_bootstrap
  bmr split 2_training/${target}.result.filter.tsv \
    -c 6000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad mci hc; do
for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
  echo ${f}
  bsub -n 24 -J "bs-${f}" \
    bash multibootstrap.sh 2_${target}_bootstrap ${f} ${target} 0.6 0.4 2_${target}_bootstrap
done
done
## ......
## ......
## ......

rm output.*
rm -fr ad_split mci_split hc_split

for target in ad mci hc; do
  tsv-append -H 2_${target}_bootstrap/*.count.tsv \
    >2_${target}_bootstrap/${target}.result.count
  tsv-append -H 2_${target}_bootstrap/*.bootstrap.tsv \
    >2_${target}_bootstrap/result.tsv
done

for target in ad mci hc; do
  bash result_stat.sh -b 2_${target}_bootstrap/result.tsv
done
```

| #Item                     | Value        |
|---------------------------|--------------|
| 2_ad_bootstrap/result.tsv | 1623749      |
| count                     | 3711         |
| reg_p_median              | 4.64683e-05  |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.60001      |
| rocauc_max                | 0.76508      |
| testauc_min               | 0.60004      |
| testauc_max               | 0.72992      |
| bin(BS)                   | count        |
| 100                       | 77346        |
| 95                        | 526142       |
| 90                        | 486077       |
| 85                        | 321524       |
| 80                        | 149817       |
| 75                        | 48528        |
| 70                        | 10942        |
| 65                        | 2510         |
| 60                        | 641          |
| 55                        | 189          |
| 50                        | 30           |
| 45                        | 2            |

| #Item                      | Value        |
|----------------------------|--------------|
| 2_mci_bootstrap/result.tsv | 117399       |
| count                      | 763          |
| reg_p_median               | 0.0002388197 |
| reg_p_min                  | 4e-10        |
| rocauc_min                 | 0.60058      |
| rocauc_max                 | 0.75878      |
| testauc_min                | 0.60007      |
| testauc_max                | 0.73671      |
| bin(BS)                    | count        |
| 100                        | 11950        |
| 95                         | 61938        |
| 90                         | 29161        |
| 85                         | 9938         |
| 80                         | 3054         |
| 75                         | 943          |
| 70                         | 263          |
| 65                         | 117          |
| 60                         | 27           |
| 55                         | 6            |
| 50                         | 1            |

| #Item                     | Value        |
|---------------------------|--------------|
| 2_hc_bootstrap/result.tsv | 1011718      |
| count                     | 4501         |
| reg_p_median              | 2.67413e-05  |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.60001      |
| rocauc_max                | 0.71269      |
| testauc_min               | 0.60003      |
| testauc_max               | 0.70028      |
| bin(BS)                   | count        |
| 100                       | 12311        |
| 95                        | 153737       |
| 90                        | 207380       |
| 85                        | 201987       |
| 80                        | 165595       |
| 75                        | 123291       |
| 70                        | 78712        |
| 65                        | 41085        |
| 60                        | 18820        |
| 55                        | 7142         |
| 50                        | 1505         |
| 45                        | 148          |
| 40                        | 4            |

```shell
BS_PASS=95
for target in ad mci hc; do
    tsv-filter -H --ge 8:${BS_PASS} \
      <2_"${target}"_bootstrap/result.tsv >2_training/"${target}".bs.tsv
done

#for target in ad mci hc; do
#  bash select_col.sh -f 1-3 \
#    training.tsv.gz 2_training/"${target}".bs.tsv \
#    >2_"${target}".training.tsv
#  bash select_col.sh -f 1-3 \
#    testing.tsv.gz 2_training/"${target}".bs.tsv \
#    >2_"${target}".testing.tsv
#  sed -i "s:Row.names:#sample:g" 2_"${target}".training.tsv
#  sed -i "s:Row.names:#sample:g" 2_"${target}".testing.tsv
#done
```

```shell
mkdir -p 3_training

for target in ad mci hc; do
  bmr nextstep 2_training/"${target}".bs.tsv 2_training/"${target}".data.replace.tsv |
    tsv-uniq >3_training/"${target}".formula.tsv
done

for target in ad mci hc; do
  bmr split 3_training/"${target}".formula.tsv \
    -c 18000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad mci hc; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 400 -a 3 -d - "${target}"_job/
done

for target in ad mci hc; do
  for f in $(find "${target}"_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo "${f}"
    bsub -n 24 -J "training-${f}" \
      "
        parallel --no-run-if-empty --line-buffer -k -j 24 '
        echo '\''==> Processing {}'\''
        Rscript multivariate.R ${target} 2_training/${target}.data.tsv 2_testing/${target}.data.tsv {} {}.result.tsv
    ' < ${f}
    "
  done
done

for target in ad mci hc; do
  tsv-append -H "${target}"_split/*.result.tsv >2_training/"${target}".result.tsv
done
rm -fr ad_split mci_split hc_split

keep-header -- awk '($6>0.67&&$7>0.67)||($6<0.33&&$7<0.33)' \
  <3_training/ad.result.tsv \
  >3_training/ad.result.filter.tsv
keep-header -- awk '($6>0.67&&$7>0.67)||($6<0.33&&$7<0.33)' \
  <3_training/mci.result.tsv \
  >3_training/mci.result.filter.tsv
keep-header -- awk '($6>0.67&&$7>0.67)||($6<0.33&&$7<0.33)' \
  <3_training/mci.result.tsv \
  >3_training/mci.result.filter.tsv
  
for target in ad mci hc; do
  bash result_stat.sh 3_training/${target}.result.filter.tsv
done
```

```shell
for target in ad mci hc; do
  mkdir -p 3_${target}_bootstrap
  bash BS.sh 2_training/${target}.data.tsv 365 3_${target}_bootstrap
  bmr split 3_training/${target}.result.filter.tsv \
    -c 15000 --mode row --rr 1 \
    -o ${target}_split | bash
done

for target in ad mci hc; do
for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
  echo ${f}
  bsub -n 24 -J "bs-${f}" \
    bash multibootstrap.sh 3_${target}_bootstrap ${f} ${target} 0.67 0.33 3_${target}_bootstrap
done
done

for target in ad mci hc; do
  tsv-append -H 3_${target}_bootstrap/*.count.tsv \
    >3_${target}_bootstrap/${target}.result.count
  tsv-append -H 3_${target}_bootstrap/*.bootstrap.tsv \
    >3_${target}_bootstrap/result.tsv
done

for target in ad mci hc; do
  bash result_stat.sh -b 3_${target}_bootstrap/result.tsv
done


```

| #Item                      | Value        |
|----------------------------|--------------|
| 3_mci_bootstrap/result.tsv | 709384       |
| count                      | 763          |
| reg_p_median               | 5.2193e-06   |
| reg_p_min                  | 0.0000000000 |
| rocauc_min                 | 0.67004      |
| rocauc_max                 | 0.79368      |
| testauc_min                | 0.67004      |
| testauc_max                | 0.76927      |
| bin(BS)                    | count        |
| 100                        | 6376         |
| 95                         | 89783        |
| 90                         | 137486       |
| 85                         | 143375       |
| 80                         | 125610       |
| 75                         | 94597        |
| 70                         | 61446        |
| 65                         | 33093        |
| 60                         | 13405        |
| 55                         | 3641         |
| 50                         | 525          |
| 45                         | 43           |
| 40                         | 3            |