```shell
mkdir -p 1_training

bmr split training.tsv.gz -c 2000 --mode column --rc 3 -o split | bash
bmr split testing.tsv.gz -c 2000 --mode column --rc 3 -o split | bash

mkdir -p job

find split -type f -name "*[0-9]" |
  sed "s:split/training.::g" |
  sed "s:split/testing.::g" |
  sort | uniq |
  split -l 240 -a 3 -d - job/

mkdir -p 1_ad_training 1_mci_training 1_hc_training

for target in ad mci hc; do
  for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo "${target}-${f}"
    bsub -n 24 -J "training-${f}" "
    parallel --no-run-if-empty --line-buffer -k -j 24 '
      echo '\''==> Processing {}'\''
      Rscript univariate.R ${target} split/training.{} split/testing.{} 1_${target}_training/{}.tsv
    ' <${f} "
  done
done
```

```shell
parallel --xapply -j 3 '
  tsv-append -H 1_{}_training/*.tsv >1_training/{}.result.tsv
' ::: ad mci hc

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

#0.58587 0.56711
#0.58339 0.57096
#0.58091 0.56184

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
for target in ad mci hc; do
  bash select_col.sh \
    -f 1-3 testing.tsv.gz 1_training/${target}.result.filter.tsv \
    >1_training/${target}.testing.data.tsv
done
```

```shell
for target in ad mci hc; do
  bash BS.sh 1_training/${target}.data.tsv 485 1_${target}_bootstrap
done
for target in ad mci hc; do
  bash BS.sh 1_training/${target}.testing.data.tsv 362 1_${target}_test_bootstrap
done

for target in ad mci hc; do
  bsub -n 24 -q mpi -J "bs_${target}" \
    bash unibootstrap.sh \
    1_${target}_bootstrap \
    1_training/${target}.result.filter.tsv \
    ${target} 0.55 0.45 \
    1_${target}_bootstrap
done
for target in ad mci hc; do
  bsub -n 24 -q mpi -J "bs_${target}" \
    bash unibootstrap.sh \
    1_${target}_test_bootstrap \
    1_training/${target}.result.filter.tsv \
    ${target} 0.50 0 \
    1_${target}_test_bootstrap
done

for target in ad mci hc; do
  mv 1_${target}_bootstrap/${target}.result.filter.tsv.bootstrap.tsv \
    1_${target}_bootstrap/result.tsv
done

for target in ad mci hc; do
  bash result_stat.sh -b 1_${target}_bootstrap/result.tsv
done
```

| #Item                     | Value         |
|---------------------------|---------------|
| 1_ad_bootstrap/result.tsv | 8833          |
| count                     | 8832          |
| reg_p_median              | 0.00466343225 |
| reg_p_min                 | 3.4e-09       |
| rocauc_min                | 0.39476       |
| rocauc_max                | 0.71064       |
| testauc_min               | 0.34619       |
| testauc_max               | 0.73506       |
| bin(BS)                   | count         |
| 100                       | 221           |
| 95                        | 1882          |
| 90                        | 2595          |
| 85                        | 2548          |
| 80                        | 1241          |
| 75                        | 125           |
| 70                        | 31            |
| 65                        | 41            |
| 60                        | 57            |
| 55                        | 45            |
| 50                        | 33            |
| 45                        | 10            |
| 40                        | 3             |

| #Item                      | Value        |
|----------------------------|--------------|
| 1_mci_bootstrap/result.tsv | 2903         |
| count                      | 2902         |
| reg_p_median               | 0.0442313075 |
| reg_p_min                  | 6.327e-06    |
| rocauc_min                 | 0.36479      |
| rocauc_max                 | 0.70320      |
| testauc_min                | 0.37427      |
| testauc_max                | 0.68788      |
| bin(BS)                    | count        |
| 100                        | 34           |
| 95                         | 380          |
| 90                         | 540          |
| 85                         | 707          |
| 80                         | 623          |
| 75                         | 354          |
| 70                         | 115          |
| 65                         | 34           |
| 60                         | 35           |
| 55                         | 31           |
| 50                         | 31           |
| 45                         | 13           |
| 40                         | 4            |
| 35                         | 1            |

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
| 100                       | 627          |
| 95                        | 4259         |
| 90                        | 3904         |
| 85                        | 1708         |
| 80                        | 200          |
| 75                        | 14           |
| 70                        | 10           |
| 65                        | 15           |
| 60                        | 6            |
| 55                        | 11           |
| 50                        | 13           |
| 45                        | 6            |

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
        # shellcheck disable=SC2027
        Rscript multivariate.R ${target} 1_${target}.training.tsv 1_${target}.testing.tsv {} {}.result.tsv
    ' < ${f}
    "
  done
done

parallel --xapply -j 3 '
  tsv-append -H {}_split/*.result.tsv >2_training/{}.result.tsv
' ::: ad mci hc
rm -fr ./*_split ./*_job
rm ./output.*

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

#0.66758 0.63608
#0.69595 0.65052
#0.64218 0.61025

## reset filter for mci
keep-header -- awk \
  -va1=0.64218 -va2=0.4 -va3=0.61025 -va4=0.4 \
  '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
  <2_training/mci.result.tsv \
  >2_training/mci.result.filter.tsv

for target in ad mci hc; do
  bash result_stat.sh 2_training/${target}.result.filter.tsv
done
```

| #Item                           | Value        |
|---------------------------------|--------------|
| 2_training/ad.result.filter.tsv | 12507        |
| count                           | 1872         |
| reg_p_median                    | 7.583e-07    |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.66763      |
| rocauc_max                      | 0.73982      |
| testauc_min                     | 0.63613      |
| testauc_max                     | 0.72083      |

| #Item                            | Value        |
|----------------------------------|--------------|
| 2_training/mci.result.filter.tsv | 45476        |
| count                            | 414          |
| reg_p_median                     | 0.0002277964 |
| reg_p_min                        | 4.2e-09      |
| rocauc_min                       | 0.64219      |
| rocauc_max                       | 0.75770      |
| testauc_min                      | 0.61027      |
| testauc_max                      | 0.74421      |

| #Item                           | Value      |
|---------------------------------|------------|
| 2_training/hc.result.filter.tsv | 97149      |
| count                           | 4877       |
| reg_p_median                    | 1.5559e-06 |
| reg_p_min                       | 1e-10      |
| rocauc_min                      | 0.64221    |
| rocauc_max                      | 0.70250    |
| testauc_min                     | 0.61029    |
| testauc_max                     | 0.69717    |

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
| 2_ad_bootstrap/result.tsv | 12507        |
| count                     | 1872         |
| reg_p_median              | 7.583e-07    |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.66763      |
| rocauc_max                | 0.73982      |
| testauc_min               | 0.63613      |
| testauc_max               | 0.72083      |
| bin(BS)                   | count        |
| 100                       | 4328         |
| 95                        | 8046         |
| 90                        | 124          |
| 85                        | 6            |
| 80                        | 1            |
| 75                        | 1            |

| #Item                      | Value        |
|----------------------------|--------------|
| 2_mci_bootstrap/result.tsv | 45476        |
| count                      | 414          |
| reg_p_median               | 0.0002277964 |
| reg_p_min                  | 4.2e-09      |
| rocauc_min                 | 0.64219      |
| rocauc_max                 | 0.75770      |
| testauc_min                | 0.61027      |
| testauc_max                | 0.74421      |
| bin(BS)                    | count        |
| 100                        | 5827         |
| 95                         | 26421        |
| 90                         | 10898        |
| 85                         | 2140         |
| 80                         | 181          |
| 75                         | 6            |
| 70                         | 2            |

| #Item                     | Value      |
|---------------------------|------------|
| 2_hc_bootstrap/result.tsv | 97149      |
| count                     | 4877       |
| reg_p_median              | 1.5559e-06 |
| reg_p_min                 | 1e-10      |
| rocauc_min                | 0.64221    |
| rocauc_max                | 0.70250    |
| testauc_min               | 0.61029    |
| testauc_max               | 0.69717    |
| bin(BS)                   | count      |
| 100                       | 4973       |
| 95                        | 61845      |
| 90                        | 28562      |
| 85                        | 1743       |
| 80                        | 22         |
| 75                        | 3          |

```shell
BS_PASS=95
for target in ad mci hc; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <2_"${target}"_bootstrap/result.tsv >2_training/"${target}".bs.tsv
done

for target in ad mci hc; do
  bash select_col.sh -f 1-3 \
    training.tsv.gz 2_training/"${target}".bs.tsv \
    >2_"${target}".training.tsv
  bash select_col.sh -f 1-3 \
    testing.tsv.gz 2_training/"${target}".bs.tsv \
    >2_"${target}".testing.tsv
  sed -i "s:Row.names:#sample:g" 2_"${target}".training.tsv
  sed -i "s:Row.names:#sample:g" 2_"${target}".testing.tsv
done
```

```shell
mkdir -p 3_training

for target in ad mci hc; do
  sed -i "s:ProbeID:#marker:g" 2_training/"${target}".bs.tsv
  bmr nextstep 2_training/"${target}".bs.tsv 1_training/"${target}".bs.tsv |
    tsv-uniq >3_training/"${target}".formula.tsv
done

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

parallel --xapply -j 3 '
  tsv-append -H {}_split/*.result.tsv >3_training/{}.result.tsv
' ::: ad mci hc
rm -fr ./*_split ./*_job

for target in ad mci hc; do
  train_upper=$(cut -f 6 3_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 3_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.35 -va3="${test_upper}" -va4=0.35 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <3_training/${target}.result.tsv \
    >3_training/${target}.result.filter.tsv
done

#0.70715 0.67691
#0.73012 0.67565
#0.67515 0.63940

## reset filter for hc
keep-header -- awk \
  -va1=0.67515 -va2=0.35 -va3=0.65 -va4=0.35 \
  '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
  <3_training/hc.result.tsv \
  >3_training/hc.result.filter.tsv

for target in ad mci hc; do
  bash result_stat.sh 3_training/${target}.result.filter.tsv
done
```

| #Item                           | Value        |
|---------------------------------|--------------|
| 3_training/ad.result.filter.tsv | 61752        |
| count                           | 2077         |
| reg_p_median                    | 9e-10        |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.70719      |
| rocauc_max                      | 0.77861      |
| testauc_min                     | 0.67695      |
| testauc_max                     | 0.74078      |

| #Item                            | Value        |
|----------------------------------|--------------|
| 3_training/mci.result.filter.tsv | 20236        |
| count                            | 409          |
| reg_p_median                     | 1.267e-07    |
| reg_p_min                        | 0.0000000000 |
| rocauc_min                       | 0.73018      |
| rocauc_max                       | 0.78824      |
| testauc_min                      | 0.67572      |
| testauc_max                      | 0.75441      |

| #Item                           | Value        |
|---------------------------------|--------------|
| 3_training/hc.result.filter.tsv | 369987       |
| count                           | 4886         |
| reg_p_median                    | 6.3e-09      |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.67519      |
| rocauc_max                      | 0.72799      |
| testauc_min                     | 0.65003      |
| testauc_max                     | 0.71784      |

```shell
for target in ad mci hc; do
  mkdir -p 3_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 365 3_${target}_bootstrap
  bmr split 3_training/${target}.result.filter.tsv \
    -c 6000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad mci hc; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh 3_${target}_bootstrap ${f} ${target} 0.65 0.35 3_${target}_bootstrap
  done
done
## ......
## ......
## ......

rm ./output.*
rm -fr ./*_split

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

| #Item                     | Value        |
|---------------------------|--------------|
| 3_ad_bootstrap/result.tsv | 61752        |
| count                     | 2077         |
| reg_p_median              | 9e-10        |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.70719      |
| rocauc_max                | 0.77861      |
| testauc_min               | 0.67695      |
| testauc_max               | 0.74078      |
| bin(BS)                   | count        |
| 100                       | 11209        |
| 95                        | 43764        |
| 90                        | 6559         |
| 85                        | 213          |
| 80                        | 6            |

| #Item                      | Value        |
|----------------------------|--------------|
| 3_mci_bootstrap/result.tsv | 20236        |
| count                      | 409          |
| reg_p_median               | 1.267e-07    |
| reg_p_min                  | 0.0000000000 |
| rocauc_min                 | 0.73018      |
| rocauc_max                 | 0.78824      |
| testauc_min                | 0.67572      |
| testauc_max                | 0.75441      |
| bin(BS)                    | count        |
| 100                        | 9931         |
| 95                         | 10277        |
| 90                         | 27           |

| #Item                     | Value        |
|---------------------------|--------------|
| 3_hc_bootstrap/result.tsv | 369987       |
| count                     | 4886         |
| reg_p_median              | 6.3e-09      |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.67519      |
| rocauc_max                | 0.72799      |
| testauc_min               | 0.65003      |
| testauc_max               | 0.71784      |
| bin(BS)                   | count        |
| 100                       | 759          |
| 95                        | 42607        |
| 90                        | 134615       |
| 85                        | 146961       |
| 80                        | 41480        |
| 75                        | 3424         |
| 70                        | 140          |

```shell
BS_PASS=95
for target in ad mci hc; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <3_"${target}"_bootstrap/result.tsv >3_training/"${target}".bs.tsv
done
```

```shell
mkdir -p 4_training

for target in ad mci hc; do
  sed -i "s:ProbeID:#marker:g" 3_training/"${target}".bs.tsv
  bmr nextstep 3_training/"${target}".bs.tsv 1_training/"${target}".bs.tsv |
    tsv-uniq >4_training/"${target}".formula.tsv
done

for target in ad mci hc; do
  bmr split 4_training/"${target}".formula.tsv \
    -c 30000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad mci hc; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 300 -a 3 -d - "${target}"_job/
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

parallel --xapply -j 3 '
  tsv-append -H {}_split/*.result.tsv >4_training/{}.result.tsv
' ::: ad mci hc
rm -fr ./*_split ./*_job

for target in ad mci hc; do
  train_upper=$(cut -f 6 4_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 4_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.3 -va3="${test_upper}" -va4=0.3 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <4_training/${target}.result.tsv \
    >4_training/${target}.result.filter.tsv
done

#0.74358 0.70111
#0.77650 0.72334
#0.70845 0.67070

## reset filter for ad
keep-header -- awk \
  -va1=0.74358 -va2=0.3 -va3=0.71 -va4=0.3 \
  '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
  <4_training/ad.result.tsv \
  >4_training/ad.result.filter.tsv

for target in ad mci hc; do
  bash result_stat.sh 4_training/${target}.result.filter.tsv
done
```

| #Item                           | Value        |
|---------------------------------|--------------|
| 4_training/ad.result.filter.tsv | 129104       |
| count                           | 2095         |
| reg_p_median                    | 0            |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.74362      |
| rocauc_max                      | 0.79959      |
| testauc_min                     | 0.71002      |
| testauc_max                     | 0.75603      |

| #Item                            | Value        |
|----------------------------------|--------------|
| 4_training/mci.result.filter.tsv | 12252        |
| count                            | 398          |
| reg_p_median                     | 4e-10        |
| reg_p_min                        | 0.0000000000 |
| rocauc_min                       | 0.77656      |
| rocauc_max                       | 0.81455      |
| testauc_min                      | 0.72340      |
| testauc_max                      | 0.77231      |

| #Item                           | Value        |
|---------------------------------|--------------|
| 4_training/hc.result.filter.tsv | 889086       |
| count                           | 4886         |
| reg_p_median                    | 0            |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.70849      |
| rocauc_max                      | 0.75371      |
| testauc_min                     | 0.67073      |
| testauc_max                     | 0.72575      |

```shell
for target in ad mci hc; do
  mkdir -p 4_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 365 4_${target}_bootstrap
  bmr split 4_training/${target}.result.filter.tsv \
    -c 12000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad mci hc; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh 4_${target}_bootstrap ${f} ${target} 0.7 0.3 4_${target}_bootstrap
  done
done
## ......
## ......
## ......

rm ./output.*
rm -fr ./*_split

parallel --xapply -j 3 '
  tsv-append -H 4_{}_bootstrap/*.count.tsv \
    >4_{}_bootstrap/{}.result.count
  tsv-append -H 4_{}_bootstrap/*.bootstrap.tsv \
    >4_{}_bootstrap/result.tsv
' ::: ad mci hc

for target in ad mci hc; do
  bash result_stat.sh -b 4_${target}_bootstrap/result.tsv
done
```

| #Item                     | Value        |
|---------------------------|--------------|
| 4_ad_bootstrap/result.tsv | 129104       |
| count                     | 2095         |
| reg_p_median              | 0            |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.74362      |
| rocauc_max                | 0.79959      |
| testauc_min               | 0.71002      |
| testauc_max               | 0.75603      |
| bin(BS)                   | count        |
| 100                       | 7074         |
| 95                        | 97139        |
| 90                        | 24735        |
| 85                        | 155          |

| #Item                      | Value        |
|----------------------------|--------------|
| 4_mci_bootstrap/result.tsv | 12252        |
| count                      | 398          |
| reg_p_median               | 4e-10        |
| reg_p_min                  | 0.0000000000 |
| rocauc_min                 | 0.77656      |
| rocauc_max                 | 0.81455      |
| testauc_min                | 0.72340      |
| testauc_max                | 0.77231      |
| bin(BS)                    | count        |
| 100                        | 7393         |
| 95                         | 4855         |
| 90                         | 3            |

| #Item                     | Value        |
|---------------------------|--------------|
| 4_hc_bootstrap/result.tsv | 889086       |
| count                     | 4886         |
| reg_p_median              | 0            |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.70849      |
| rocauc_max                | 0.75371      |
| testauc_min               | 0.67073      |
| testauc_max               | 0.72575      |
| bin(BS)                   | count        |
| 100                       | 5            |
| 95                        | 635          |
| 90                        | 5937         |
| 85                        | 26477        |
| 80                        | 80332        |
| 75                        | 180455       |
| 70                        | 275350       |
| 65                        | 229219       |
| 60                        | 78806        |
| 55                        | 11230        |
| 50                        | 629          |
| 45                        | 10           |

```shell
BS_PASS=95
for target in ad mci hc; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <4_"${target}"_bootstrap/result.tsv >4_training/"${target}".bs.tsv
done

## BS_PASS=85 for hc
```

```shell
mkdir -p 5_training

for target in ad mci hc; do
  sed -i "s:ProbeID:#marker:g" 4_training/"${target}".bs.tsv
  bmr nextstep 4_training/"${target}".bs.tsv 1_training/"${target}".bs.tsv |
    tsv-uniq >5_training/"${target}".formula.tsv
done

for target in ad mci hc; do
  bmr split 5_training/"${target}".formula.tsv \
    -c 50000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad mci hc; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 240 -a 3 -d - "${target}"_job/
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

parallel --xapply -j 3 '
  tsv-append -H {}_split/*.result.tsv >5_training/{}.result.tsv
' ::: ad mci hc
rm -fr ./*_split ./*_job

for target in ad mci hc; do
  train_upper=$(cut -f 6 5_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 5_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.25 -va3="${test_upper}" -va4=0.25 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <5_training/${target}.result.tsv \
    >5_training/${target}.result.filter.tsv
done

#0.77102 0.72739
#0.81200 0.75691
#0.73822 0.68933

## reset filter for hc
keep-header -- awk \
  -va1=0.77102 -va2=0.25 -va3=0.74 -va4=0.25 \
  '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
  <5_training/ad.result.tsv \
  >5_training/ad.result.filter.tsv
  
keep-header -- awk \
  -va1=0.8 -va2=0.25 -va3=0.75691 -va4=0.25 \
  '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
  <5_training/mci.result.tsv \
  >5_training/mci.result.filter.tsv

for target in ad mci hc; do
  bash result_stat.sh 5_training/${target}.result.filter.tsv
done
```

| #Item                           | Value        |
|---------------------------------|--------------|
| 5_training/ad.result.filter.tsv | 56509        |
| count                           | 1968         |
| reg_p_median                    | 0            |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.77106      |
| rocauc_max                      | 0.81809      |
| testauc_min                     | 0.74003      |
| testauc_max                     | 0.76955      |

| #Item                            | Value        |
|----------------------------------|--------------|
| 5_training/mci.result.filter.tsv | 54589        |
| count                            | 406          |
| reg_p_median                     | 0            |
| reg_p_min                        | 0.0000000000 |
| rocauc_min                       | 0.80005      |
| rocauc_max                       | 0.84348      |
| testauc_min                      | 0.75697      |
| testauc_max                      | 0.80311      |

| #Item                           | Value        |
|---------------------------------|--------------|
| 5_training/hc.result.filter.tsv | 527603       |
| count                           | 4886         |
| reg_p_median                    | 0            |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.73825      |
| rocauc_max                      | 0.78082      |
| testauc_min                     | 0.68936      |
| testauc_max                     | 0.73711      |

```shell
for target in ad mci hc; do
  mkdir -p 5_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 365 5_${target}_bootstrap
  bmr split 5_training/${target}.result.filter.tsv \
    -c 10000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad mci; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh 5_${target}_bootstrap ${f} ${target} 0.75 0.25 5_${target}_bootstrap
  done
done
for target in hc; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh 5_${target}_bootstrap ${f} ${target} 0.73822 0.25 5_${target}_bootstrap
  done
done
## ......
## ......
## ......

rm ./output.*
rm -fr ./*_split

parallel --xapply -j 3 '
  tsv-append -H 5_{}_bootstrap/*.count.tsv \
    >5_{}_bootstrap/{}.result.count
  tsv-append -H 5_{}_bootstrap/*.bootstrap.tsv \
    >5_{}_bootstrap/result.tsv
' ::: ad mci hc

for target in ad mci hc; do
  bash result_stat.sh -b 5_${target}_bootstrap/result.tsv
done
```

| #Item                     | Value        |
|---------------------------|--------------|
| 5_ad_bootstrap/result.tsv | 56509        |
| count                     | 1968         |
| reg_p_median              | 0            |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.77106      |
| rocauc_max                | 0.81809      |
| testauc_min               | 0.74003      |
| testauc_max               | 0.76955      |
| bin(BS)                   | count        |
| 100                       | 211          |
| 95                        | 6664         |
| 90                        | 23981        |
| 85                        | 19763        |
| 80                        | 5518         |
| 75                        | 366          |
| 70                        | 5            |

| #Item                      | Value        |
|----------------------------|--------------|
| 5_mci_bootstrap/result.tsv | 54589        |
| count                      | 406          |
| reg_p_median               | 0            |
| reg_p_min                  | 0.0000000000 |
| rocauc_min                 | 0.80005      |
| rocauc_max                 | 0.84348      |
| testauc_min                | 0.75697      |
| testauc_max                | 0.80311      |
| bin(BS)                    | count        |
| 100                        | 8881         |
| 95                         | 41809        |
| 90                         | 3793         |
| 85                         | 105          |

| #Item                     | Value        |
|---------------------------|--------------|
| 5_hc_bootstrap/result.tsv | 527603       |
| count                     | 4886         |
| reg_p_median              | 0            |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.73825      |
| rocauc_max                | 0.78082      |
| testauc_min               | 0.68936      |
| testauc_max               | 0.73711      |
| bin(BS)                   | count        |
| 90                        | 75           |
| 85                        | 746          |
| 80                        | 3864         |
| 75                        | 14588        |
| 70                        | 45699        |
| 65                        | 110669       |
| 60                        | 180342       |
| 55                        | 137384       |
| 50                        | 32053        |
| 45                        | 2151         |
| 40                        | 31           |

```shell
BS_PASS=85
for target in ad mci; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <5_"${target}"_bootstrap/result.tsv >5_training/"${target}".bs.tsv
done

## BS_PASS=70 for hc
```

```shell
mkdir -p 6_training

for target in ad mci hc; do
  sed -i "s:ProbeID:#marker:g" 5_training/"${target}".bs.tsv
  bmr nextstep 5_training/"${target}".bs.tsv 1_training/"${target}".bs.tsv |
    tsv-uniq >6_training/"${target}".formula.tsv
done

for target in ad mci hc; do
  bmr split 6_training/"${target}".formula.tsv \
    -c 40000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad mci hc; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 240 -a 3 -d - "${target}"_job/
done

for target in ad mci hc; do
  for f in $(find "${target}"_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo "${f}"
    bsub -n 24 -J "training-${f}" \
      "
        parallel --no-run-if-empty --line-buffer -k -j 24 '
        echo '\''==> Processing {}'\''
        Rscript multivariate.R ${target} \
          1_${target}.training.tsv \
          1_${target}.testing.tsv {} {}.result.tsv
    ' < ${f}
    "
  done
done

parallel --xapply -j 3 '
  tsv-append -H {}_split/*.result.tsv >6_training/{}.result.tsv
' ::: ad mci hc
rm -fr ./*_split ./*_job

for target in ad mci hc; do
  train_upper=$(cut -f 6 6_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 6_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.2 -va3="${test_upper}" -va4=0.2 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <6_training/${target}.result.tsv \
    >6_training/${target}.result.filter.tsv
done

#0.77102 0.72739
#0.81200 0.75691
#0.73822 0.68933

## reset filter for hc
keep-header -- awk \
  -va1=0.77102 -va2=0.25 -va3=0.74 -va4=0.25 \
  '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
  <5_training/ad.result.tsv \
  >5_training/ad.result.filter.tsv
  
keep-header -- awk \
  -va1=0.8 -va2=0.25 -va3=0.75691 -va4=0.25 \
  '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
  <5_training/mci.result.tsv \
  >5_training/mci.result.filter.tsv

for target in ad mci hc; do
  bash result_stat.sh 6_training/${target}.result.filter.tsv
done
```

| #Item                           | Value        |
|---------------------------------|--------------|
| 6_training/ad.result.filter.tsv | 360417       |
| count                           | 2099         |
| reg_p_median                    | 0            |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.79518      |
| rocauc_max                      | 0.84051      |
| testauc_min                     | 0.75235      |
| testauc_max                     | 0.78728      |

| #Item                            | Value        |
|----------------------------------|--------------|
| 6_training/mci.result.filter.tsv | 30666        |
| count                            | 407          |
| reg_p_median                     | 0            |
| reg_p_min                        | 0.0000000000 |
| rocauc_min                       | 0.83267      |
| rocauc_max                       | 0.86711      |
| testauc_min                      | 0.78244      |
| testauc_max                      | 0.81601      |

| #Item                           | Value        |
|---------------------------------|--------------|
| 6_training/hc.result.filter.tsv | 1066678      |
| count                           | 4886         |
| reg_p_median                    | 0            |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.75897      |
| rocauc_max                      | 0.79979      |
| testauc_min                     | 0.70400      |
| testauc_max                     | 0.74153      |

```shell
for target in ad mci hc; do
  mkdir -p 6_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 365 6_${target}_bootstrap
  bmr split 6_training/${target}.result.filter.tsv \
    -c 10000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad mci; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh \
      6_${target}_bootstrap \
      ${f} ${target} 0.8 0.2 \
      6_${target}_bootstrap
  done
done
for target in hc; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -J "bs-${f}" \
      bash multibootstrap.sh \
      6_${target}_bootstrap \
      ${f} ${target} 0.75894 0.2 \
      6_${target}_bootstrap
  done
done
## ......
## ......
## ......

rm ./output.*
rm -fr ./*_split

parallel --xapply -j 3 '
  tsv-append -H 6_{}_bootstrap/*.count.tsv \
    >6_{}_bootstrap/{}.result.count
  tsv-append -H 6_{}_bootstrap/*.bootstrap.tsv \
    >6_{}_bootstrap/result.tsv
' ::: ad mci hc

for target in ad mci hc; do
  bash result_stat.sh -b 6_${target}_bootstrap/result.tsv
done
```

| #Item                     | Value        |
|---------------------------|--------------|
| 6_ad_bootstrap/result.tsv | 360417       |
| count                     | 2099         |
| reg_p_median              | 0            |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.79518      |
| rocauc_max                | 0.84051      |
| testauc_min               | 0.75235      |
| testauc_max               | 0.78728      |
| bin(BS)                   | count        |
| 95                        | 15           |
| 90                        | 301          |
| 85                        | 3122         |
| 80                        | 8438         |
| 75                        | 20660        |
| 70                        | 35671        |
| 65                        | 70283        |
| 60                        | 115720       |
| 55                        | 87667        |
| 50                        | 17399        |
| 45                        | 1086         |
| 40                        | 54           |

| #Item                      | Value        |
|----------------------------|--------------|
| 6_mci_bootstrap/result.tsv | 30666        |
| count                      | 407          |
| reg_p_median               | 0            |
| reg_p_min                  | 0.0000000000 |
| rocauc_min                 | 0.83267      |
| rocauc_max                 | 0.86711      |
| testauc_min                | 0.78244      |
| testauc_max                | 0.81601      |
| bin(BS)                    | count        |
| 100                        | 1312         |
| 95                         | 20282        |
| 90                         | 8339         |
| 85                         | 699          |
| 80                         | 32           |
| 75                         | 1            |

| #Item                     | Value        |
|---------------------------|--------------|
| 6_hc_bootstrap/result.tsv | 1066678      |
| count                     | 4886         |
| reg_p_median              | 0            |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.75897      |
| rocauc_max                | 0.79979      |
| testauc_min               | 0.70400      |
| testauc_max               | 0.74153      |
| bin(BS)                   | count        |
| 95                        | 118          |
| 90                        | 1184         |
| 85                        | 8834         |
| 80                        | 26542        |
| 75                        | 93680        |
| 70                        | 231520       |
| 65                        | 378231       |
| 60                        | 257413       |
| 55                        | 64937        |
| 50                        | 4139         |
| 45                        | 79           |

```shell
BS_PASS=75
for target in ad mci; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <6_"${target}"_bootstrap/result.tsv >6_training/"${target}".bs.tsv
done
BS_PASS=80
for target in hc; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <6_"${target}"_bootstrap/result.tsv >6_training/"${target}".bs.tsv
done
```

```shell
mkdir -p 7_training

for target in ad mci hc; do
  sed -i "s:ProbeID:#marker:g" 6_training/"${target}".bs.tsv
  bmr nextstep 6_training/"${target}".bs.tsv 1_training/"${target}".bs.tsv |
    tsv-uniq >7_training/"${target}".formula.tsv
done

for target in ad mci hc; do
  bmr split 7_training/"${target}".formula.tsv \
    -c 12000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad mci hc; do
  mkdir -p "${target}"_job
  find "${target}"_split -type f -name "*[0-9]" |
    sort |
    split -l 240 -a 3 -d - "${target}"_job/
done

for target in ad mci hc; do
  for f in $(find "${target}"_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo "${f}"
    bsub -n 24 -q mpi -J "training-${f}" \
      "
        parallel --no-run-if-empty --line-buffer -k -j 24 '
        echo '\''==> Processing {}'\''
        Rscript multivariate.R ${target} \
          1_${target}.training.tsv \
          1_${target}.testing.tsv {} {}.result.tsv
    ' < ${f}
    "
  done
done

parallel --xapply -j 3 '
  tsv-append -H {}_split/*.result.tsv >7_training/{}.result.tsv
' ::: ad mci hc
rm -fr ./*_split ./*_job

for target in ad mci hc; do
  train_upper=$(cut -f 6 7_training/${target}.result.tsv |
    grep -v "rocauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  test_upper=$(cut -f 7 7_training/${target}.result.tsv |
    grep -v "testauc" | sort -n |
    awk '{all[NR] = $0} END{print all[int(NR*0.95 - 0.5)]}')
  echo "$train_upper" "$test_upper"
  keep-header -- awk \
    -va1="${train_upper}" -va2=0.15 -va3="${test_upper}" -va4=0.15 \
    '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
    <7_training/${target}.result.tsv \
    >7_training/${target}.result.filter.tsv
done

#0.82267 0.76361
#0.85925 0.80351
#0.78175 0.71664

for target in ad mci hc; do
  bash result_stat.sh 7_training/${target}.result.filter.tsv
done
```

| #Item                           | Value        |
|---------------------------------|--------------|
| 7_training/ad.result.filter.tsv | 301065       |
| count                           | 2099         |
| reg_p_median                    | 0            |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.82271      |
| rocauc_max                      | 0.85691      |
| testauc_min                     | 0.76365      |
| testauc_max                     | 0.79397      |

| #Item                            | Value        |
|----------------------------------|--------------|
| 7_training/mci.result.filter.tsv | 17102        |
| count                            | 391          |
| reg_p_median                     | 0            |
| reg_p_min                        | 0.0000000000 |
| rocauc_min                       | 0.85932      |
| rocauc_max                       | 0.89362      |
| testauc_min                      | 0.80358      |
| testauc_max                      | 0.83195      |

| #Item                           | Value        |
|---------------------------------|--------------|
| 7_training/hc.result.filter.tsv | 391146       |
| count                           | 4881         |
| reg_p_median                    | 0            |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.78178      |
| rocauc_max                      | 0.81214      |
| testauc_min                     | 0.71667      |
| testauc_max                     | 0.74897      |

```shell
for target in ad mci; do
  mkdir -p 7_${target}_bootstrap
  bash BS.sh 1_"${target}".training.tsv 365 7_${target}_bootstrap
  bmr split 7_training/${target}.result.filter.tsv \
    -c 8000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad mci; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -q mpi -J "bs-${f}" \
      bash multibootstrap.sh \
      7_${target}_bootstrap \
      ${f} ${target} 0.82267 0.15 \
      7_${target}_bootstrap
  done
done
for target in hc; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -q mpi -J "bs-${f}" \
      bash multibootstrap.sh \
      7_${target}_bootstrap \
      ${f} ${target} 0.78175 0.15 \
      7_${target}_bootstrap
  done
done
## ......
## ......
## ......

for target in ad mci; do
  bash BS.sh 1_"${target}".testing.tsv 362 7_${target}_test_bootstrap
  bmr split 7_training/${target}.result.filter.tsv \
    -c 8000 --mode row --rr 1 \
    -o ${target}_split | bash
done

for target in hc; do
  bash BS.sh 1_"${target}".testing.tsv 362 7_${target}_test_bootstrap
  bmr split 7_training/${target}.result.filter.tsv \
    -c 12000 --mode row --rr 1 \
    -o ${target}_split | bash
done

## HPCC jobs limit is 200, this hand-work step should pay more attention.
for target in ad mci; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -q mpi -J "bs-${f}" \
      bash multibootstrap.sh \
      7_${target}_test_bootstrap \
      ${f} ${target} 0.76361 0.15 \
      7_${target}_test_bootstrap
  done
done
for target in hc; do
  for f in $(find ${target}_split -maxdepth 1 -type f -name "*[0-9]" | sort); do
    echo ${f}
    bsub -n 24 -q mpi -J "bs-${f}" \
      bash multibootstrap.sh \
      7_${target}_test_bootstrap \
      ${f} ${target} 0.71664 0.15 \
      7_${target}_test_bootstrap
  done
done
## ......
## ......
## ......

rm ./output.*
rm -fr ./*_split

parallel --xapply -j 3 '
  tsv-append -H 7_{}_bootstrap/*.count.tsv \
    >7_{}_bootstrap/{}.result.count
  tsv-append -H 7_{}_bootstrap/*.bootstrap.tsv \
    >7_{}_bootstrap/result.tsv
' ::: ad mci hc

for target in ad mci hc; do
  bash result_stat.sh -b 7_${target}_bootstrap/result.tsv
done
```

| #Item                     | Value        |
|---------------------------|--------------|
| 7_ad_bootstrap/result.tsv | 133066       |
| count                     | 2091         |
| reg_p_median              | 0            |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.82271      |
| rocauc_max                | 0.85691      |
| testauc_min               | 0.76365      |
| testauc_max               | 0.79397      |
| bin(BS)                   | count        |
| 100                       | 1            |
| 95                        | 57           |
| 90                        | 935          |
| 85                        | 5067         |
| 80                        | 17092        |
| 75                        | 35130        |
| 70                        | 42373        |
| 65                        | 24864        |
| 60                        | 7161         |
| 55                        | 381          |
| 50                        | 4            |

| #Item                      | Value        |
|----------------------------|--------------|
| 7_mci_bootstrap/result.tsv | 17102        |
| count                      | 391          |
| reg_p_median               | 0            |
| reg_p_min                  | 0.0000000000 |
| rocauc_min                 | 0.85932      |
| rocauc_max                 | 0.89362      |
| testauc_min                | 0.80358      |
| testauc_max                | 0.83195      |
| bin(BS)                    | count        |
| 100                        | 3240         |
| 95                         | 12644        |
| 90                         | 1187         |
| 85                         | 30           |

| #Item                     | Value        |
|---------------------------|--------------|
| 7_hc_bootstrap/result.tsv | 391146       |
| count                     | 4881         |
| reg_p_median              | 0            |
| reg_p_min                 | 0.0000000000 |
| rocauc_min                | 0.78178      |
| rocauc_max                | 0.81214      |
| testauc_min               | 0.71667      |
| testauc_max               | 0.74897      |
| bin(BS)                   | count        |
| 95                        | 8            |
| 90                        | 262          |
| 85                        | 2405         |
| 80                        | 15549        |
| 75                        | 47609        |
| 70                        | 119199       |
| 65                        | 133701       |
| 60                        | 60725        |
| 55                        | 11326        |
| 50                        | 358          |
| 45                        | 3            |

```shell
BS_PASS=80
for target in ad mci hc; do
  tsv-filter -H --ge 8:${BS_PASS} \
    <7_"${target}"_bootstrap/result.tsv >7_training/"${target}".bs.tsv
done
```