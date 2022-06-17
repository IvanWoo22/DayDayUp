```shell
mkdir -p 1_training split job

bmr split training.tsv.gz -c 1000 --mode column --rc 3 -o split | bash

find split -type f -name "*[0-9]" |
  sed "s:split/training.::g" | sort |
  split -l 120 -a 3 -d - job/

bmr split testing.tsv.gz -c 1000 --mode column --rc 3 -o split | bash

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

## reset filter for hc
#keep-header -- awk \
#  -va1=0.70845 -va2=0.3 -va3=0.68 -va4=0.3 \
#  '($6>a1&&$7>a3)||($6<a2&&$7<a4)' \
#  <4_training/hc.result.tsv \
#  >4_training/hc.result.filter.tsv

for target in ad mci hc; do
  bash result_stat.sh 3_training/${target}.result.filter.tsv
done
```

| #Item                           | Value        |
|---------------------------------|--------------|
| 4_training/ad.result.filter.tsv | 465420       |
| count                           | 2102         |
| reg_p_median                    | 0            |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.74362      |
| rocauc_max                      | 0.79963      |
| testauc_min                     | 0.70115      |
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
| 4_training/hc.result.filter.tsv | 192487       |
| count                           | 4820         |
| reg_p_median                    | 0            |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.70849      |
| rocauc_max                      | 0.75291      |
| testauc_min                     | 0.68001      |
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