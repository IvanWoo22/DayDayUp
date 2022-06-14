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
  bash bootstrap.sh ${target}_bootstrap \
  1_training/${target}.result.filter.tsv \
  ${target} 1_${target}_bootstrap
done

for target in ad mci hc; do
  mv ${target}_bootstrap/${target}.result.filter.tsv.bootstrap.tsv \
    ${target}_bootstrap/result.tsv
  bash result_stat.sh -b ${target}_bootstrap/result.tsv
done
```

| #Item                   | Value        |
|-------------------------|--------------|
| ad_bootstrap/result.tsv | 38070        |
| count                   | 38069        |
| reg_p_median            | 0.0352870629 |
| reg_p_min               | 1.26e-08     |
| rocauc_min              | 0.39240      |
| rocauc_max              | 0.70287      |
| testauc_min             | 0.34393      |
| testauc_max             | 0.72948      |
| bin(BS)                 | count        |
| 100                     | 201          |
| 95                      | 3303         |
| 90                      | 5001         |
| 85                      | 5094         |
| 80                      | 4366         |
| 75                      | 3736         |
| 70                      | 3416         |
| 65                      | 3362         |
| 60                      | 3225         |
| 55                      | 3070         |
| 50                      | 2115         |
| 45                      | 937          |
| 40                      | 222          |
| 35                      | 20           |
| 30                      | 1            |

| #Item                    | Value         |
|--------------------------|---------------|
| mci_bootstrap/result.tsv | 29607         |
| count                    | 29606         |
| reg_p_median             | 0.16383083105 |
| reg_p_min                | 7.0081e-06    |
| rocauc_min               | 0.36271       |
| rocauc_max               | 0.70139       |
| testauc_min              | 0.37373       |
| testauc_max              | 0.68788       |
| bin(BS)                  | count         |
| 100                      | 43            |
| 95                       | 671           |
| 90                       | 1117          |
| 85                       | 1558          |
| 80                       | 2177          |
| 75                       | 2724          |
| 70                       | 3281          |
| 65                       | 3984          |
| 60                       | 4727          |
| 55                       | 4742          |
| 50                       | 3203          |
| 45                       | 1177          |
| 40                       | 189           |
| 35                       | 11            |
| 30                       | 2             |

| #Item                   | Value         |
|-------------------------|---------------|
| hc_bootstrap/result.tsv | 35393         |
| count                   | 35392         |
| reg_p_median            | 0.01772092275 |
| reg_p_min               | 1.623e-07     |
| rocauc_min              | 0.39913       |
| rocauc_max              | 0.67022       |
| testauc_min             | 0.38489       |
| testauc_max             | 0.67991       |
| bin(BS)                 | count         |
| 100                     | 243           |
| 95                      | 3202          |
| 90                      | 5894          |
| 85                      | 5086          |
| 80                      | 4037          |
| 75                      | 3273          |
| 70                      | 3010          |
| 65                      | 2783          |
| 60                      | 2716          |
| 55                      | 2529          |
| 50                      | 1745          |
| 45                      | 695           |
| 40                      | 160           |
| 35                      | 18            |
| 30                      | 1             |

```shell
BS_PASS=90
for target in ad mci hc; do
    tsv-filter -H --ge 9:${BS_PASS} \
      <"${target}"_bootstrap/result.tsv >1_training/"${target}".bs.tsv
done

for target in ad mci hc; do
  bash select_col.sh -f 1-3 \
    training.tsv.gz 1_"${target}".result.tsv \
    >1_"${target}".training.tsv
  bash select_col.sh -f 1-3 \
    testing.tsv.gz 1_"${target}".result.tsv \
    >1_"${target}".testing.tsv
done
```