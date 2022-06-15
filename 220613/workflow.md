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
  keep-header -- awk '($6>0.55&&$7>0.55)||($6<0.45&&$7<0.45)' \
    <1_training/${target}.result.tsv \
    >1_training/${target}.result.filter.tsv
done

for target in ad mci hc; do
  bash result_stat.sh 1_training/${target}.result.filter.tsv
done
```

| #Item                           | Value        |
|---------------------------------|--------------|
| 1_training/ad.result.filter.tsv | 40210        |
| count                           | 40209        |
| reg_p_median                    | 0.0311733151 |
| reg_p_min                       | 3.4e-09      |
| rocauc_min                      | 0.39476      |
| rocauc_max                      | 0.71064      |
| testauc_min                     | 0.34619      |
| testauc_max                     | 0.73506      |

| #Item                            | Value        |
|----------------------------------|--------------|
| 1_training/mci.result.filter.tsv | 29691        |
| count                            | 29690        |
| reg_p_median                     | 0.1627480769 |
| reg_p_min                        | 6.327e-06    |
| rocauc_min                       | 0.36479      |
| rocauc_max                       | 0.70320      |
| testauc_min                      | 0.37427      |
| testauc_max                      | 0.68788      |

| #Item                           | Value         |
|---------------------------------|---------------|
| 1_training/hc.result.filter.tsv | 36571         |
| count                           | 36570         |
| reg_p_median                    | 0.01576443285 |
| reg_p_min                       | 5.64e-08      |
| rocauc_min                      | 0.39883       |
| rocauc_max                      | 0.67750       |
| testauc_min                     | 0.38077       |
| testauc_max                     | 0.68239       |

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
    bash bootstrap.sh 1_${target}_bootstrap \
    1_training/${target}.result.filter.tsv \
    ${target} 1_${target}_bootstrap
done

for target in ad mci hc; do
  mv 1_${target}_bootstrap/${target}.result.filter.tsv.bootstrap.tsv \
    1_${target}_bootstrap/result.tsv
  bash result_stat.sh -b 1_${target}_bootstrap/result.tsv
done
```

| #Item                      | Value         |
|----------------------------|---------------|
| 1_ad_bootstrap/result.tsv  | 40210         |
| count                      | 40209         |
| reg_p_median               | 0.0311733151  |
| reg_p_min                  | 3.4e-09       |
| rocauc_min                 | 0.39476       |
| rocauc_max                 | 0.71064       |
| testauc_min                | 0.34619       |
| testauc_max                | 0.73506       |
| bin(BS)                    | count         |
| 100                        | 397           |
| 95                         | 3314          |
| 90                         | 5004          |
| 85                         | 5924          |
| 80                         | 5236          |
| 75                         | 4042          |
| 70                         | 3628          |
| 65                         | 3444          |
| 60                         | 3546          |
| 55                         | 3067          |
| 50                         | 1893          |
| 45                         | 602           |
| 40                         | 102           |
| 35                         | 10            |

| #Item                      | Value         |
|----------------------------|---------------|
| 1_mci_bootstrap/result.tsv | 29691         |
| count                      | 29690         |
| reg_p_median               | 0.1627480769  |
| reg_p_min                  | 6.327e-06     |
| rocauc_min                 | 0.36479       |
| rocauc_max                 | 0.70320       |
| testauc_min                | 0.37427       |
| testauc_max                | 0.68788       |
| bin(BS)                    | count         |
| 100                        | 52            |
| 95                         | 711           |
| 90                         | 1186          |
| 85                         | 1660          |
| 80                         | 2169          |
| 75                         | 2789          |
| 70                         | 3269          |
| 65                         | 3950          |
| 60                         | 4637          |
| 55                         | 4789          |
| 50                         | 3203          |
| 45                         | 1106          |
| 40                         | 156           |
| 35                         | 13            |

| #Item                     | Value         |
|---------------------------|---------------|
| 1_hc_bootstrap/result.tsv | 36571         |
| count                     | 36570         |
| reg_p_median              | 0.01576443285 |
| reg_p_min                 | 5.64e-08      |
| rocauc_min                | 0.39883       |
| rocauc_max                | 0.67750       |
| testauc_min               | 0.38077       |
| testauc_max               | 0.68239       |
| bin(BS)                   | count         |
| 100                       | 330           |
| 95                        | 4171          |
| 90                        | 5683          |
| 85                        | 5307          |
| 80                        | 4121          |
| 75                        | 3402          |
| 70                        | 2963          |
| 65                        | 2812          |
| 60                        | 2774          |
| 55                        | 2531          |
| 50                        | 1732          |
| 45                        | 628           |
| 40                        | 102           |
| 35                        | 13            |
| 30                        | 1             |

```shell
BS_PASS=90
for target in ad mci hc; do
    tsv-filter -H --ge 8:${BS_PASS} \
      <1_"${target}"_bootstrap/result.tsv >1_training/"${target}".bs.tsv
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
bmr replace 1_"${target}".training.tsv -o 2_training/"${target}".data.tsv
done

for target in ad mci hc; do
bmr nextstep 2_training/"${target}".data.replace.tsv |
tsv-uniq > 2_training/"${target}".formula.tsv
done

for target in ad mci hc; do
bmr split 2_training/"${target}".formula.tsv \
-c 6000 --mode row --rr 1 -o "${target}"_split | bash
done

for target in ad mci hc; do
mkdir -p "${target}"_job
find "${target}"_split -type f -name "*[0-9]" |
sort |
split -l 500 -a 3 -d - "${target}"_job/
done

for f in $(find ad_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
echo ${f}
bsub -q mpi -n 24 -J "training-${f}" \
"
parallel --no-run-if-empty --line-buffer -k -j 24 '
echo '\''==> Processing {}'\''
Rscript multivariate.R ad 2_training/ad.data.tsv 2_testing/ad.data.tsv {} {}.result.tsv
' < ${f}
"
done

for f in $(find mci_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
echo ${f}
bsub -q mpi -n 24 -J "training-${f}" \
"
parallel --no-run-if-empty --line-buffer -k -j 24 '
echo '\''==> Processing {}'\''
Rscript multivariate.R mci 2_training/mci.data.tsv 2_testing/mci.data.tsv {} {}.result.tsv
' < ${f}
"
done

for f in $(find hc_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
echo ${f}
bsub -q mpi -n 24 -J "training-${f}" \
"
parallel --no-run-if-empty --line-buffer -k -j 24 '
echo '\''==> Processing {}'\''
Rscript multivariate.R hc 2_training/hc.data.tsv 2_testing/hc.data.tsv {} {}.result.tsv
' < ${f}
"
done

for target in ad mci hc; do
tsv-append -H "${target}"_split/*.result.tsv > 2_training/"${target}".result.tsv
done

for target in ad mci hc; do
  keep-header -- awk '($6>0.6&&$7>0.6)||($6<0.4&&$7<0.4)' \
<2_training/${target}.result.tsv \
>2_training/${target}.result.filter.tsv
done

for target in ad mci hc; do
  bash result_stat.sh 2_training/${target}.result.filter.tsv
done
```

| #Item                           | Value        |
|---------------------------------|--------------|
| 2_training/ad.result.filter.tsv | 7412392      |
| count                           | 8715         |
| reg_p_median                    | 0.0001661535 |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.36932      |
| rocauc_max                      | 0.76508      |
| testauc_min                     | 0.35430      |
| testauc_max                     | 0.73435      |

| #Item                            | Value        |
|----------------------------------|--------------|
| 2_training/mci.result.filter.tsv | 750808       |
| count                            | 1949         |
| reg_p_median                     | 0.0007573514 |
| reg_p_min                        | 4e-10        |
| rocauc_min                       | 0.36560      |
| rocauc_max                       | 0.75878      |
| testauc_min                      | 0.36515      |
| testauc_max                      | 0.74421      |

| #Item                           | Value        |
|---------------------------------|--------------|
| 2_training/hc.result.filter.tsv | 4141460      |
| count                           | 10184        |
| reg_p_median                    | 8.06183e-05  |
| reg_p_min                       | 0.0000000000 |
| rocauc_min                      | 0.60001      |
| rocauc_max                      | 0.71269      |
| testauc_min                     | 0.60003      |
| testauc_max                     | 0.71288      |