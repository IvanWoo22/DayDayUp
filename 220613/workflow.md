```shell
for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 24 -J "training-${f}" \
    "
    parallel --no-run-if-empty --line-buffer -k -j 24 '
      echo '\''==> Processing {}'\''
      Rscript univariate.R AD split/training.{} split/testing.{} AD_result/{}.tsv
    ' <${f}
"
done
for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 24 -J "training-${f}" \
    "
    parallel --no-run-if-empty --line-buffer -k -j 24 '
      echo '\''==> Processing {}'\''
      Rscript univariate.R MCI split/training.{} split/testing.{} MCI_result/{}.tsv
    ' <${f}
"
done
for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 24 -J "training-${f}" \
    "
    parallel --no-run-if-empty --line-buffer -k -j 24 '
      echo '\''==> Processing {}'\''
      Rscript univariate.R HC split/training.{} split/testing.{} HC_result/{}.tsv
    ' <${f}
"
done
```

```shell
for target in ad mci hc; do
  keep-header -- awk '($6>0.55&&$7>0.55)||($6<0.45&&$7<0.45)' \
    <training/${target}.result.tsv \
    >training/${target}.result.filter.tsv
done

for target in ad mci hc; do
  bash result_stat.sh training/${target}.result.filter.tsv
done
```

| #Item                         | Value        |
|-------------------------------|--------------|
| training/ad.result.filter.tsv | 38070        |
| count                         | 38069        |
| reg_p_max                     | 0.9996786376 |
| reg_p_min                     | 1.26e-08     |
| rocauc_min                    | 0.39240      |
| rocauc_median                 | 0.57754      |
| testauc_max                   | 0.72948      |
| testauc_median                | 0.56622      |

| #Item                          | Value        |
|--------------------------------|--------------|
| training/mci.result.filter.tsv | 29607        |
| count                          | 29606        |
| reg_p_max                      | 0.9999540193 |
| reg_p_min                      | 7.0081e-06   |
| rocauc_min                     | 0.36271      |
| rocauc_median                  | 0.56662      |
| testauc_max                    | 0.68788      |
| testauc_median                 | 0.56623      |

| #Item                         | Value        |
|-------------------------------|--------------|
| training/hc.result.filter.tsv | 35393        |
| count                         | 35392        |
| reg_p_max                     | 0.9999244074 |
| reg_p_min                     | 1.623e-07    |
| rocauc_min                    | 0.39913      |
| rocauc_median                 | 0.57723      |
| testauc_max                   | 0.67991      |
| testauc_median                | 0.56492      |

```shell
for target in ad mci hc; do
  bash select_col.sh \
    -f 1-3 training.tsv.gz training/${target}.result.filter.tsv \
    >training/${target}.data.tsv
done
```

```shell
for target in ad mci hc; do
  bash BS.sh training/${target}.data.tsv 365 ${target}_bootstrap
done

bsub -n 24 -J "bsvalidation_ad" \
  bash bootstrap.sh ad_bootstrap \
  training/ad.result.filter.tsv \
  AD ad_bootstrap

bsub -n 24 -J "bsvalidation_mci" \
  bash bootstrap.sh mci_bootstrap \
  training/mci.result.filter.tsv \
  MCI mci_bootstrap

bsub -n 24 -J "bsvalidation_hc" \
  bash bootstrap.sh hc_bootstrap \
  training/hc.result.filter.tsv \
  HC hc_bootstrap

for target in ad mci hc; do
  mv ${target}_bootstrap/${target}.result.filter.tsv.bootstrap.tsv \
    ${target}_bootstrap/result.tsv
  bash result_stat.sh -b ${target}_bootstrap/result.tsv
done
```

| #Item                                           | Value        |
|-------------------------------------------------|--------------|
| ad_bootstrap/ad.result.filter.tsv.bootstrap.tsv | 38069        |
| count                                           | 38068        |
| 0.0117170182_max                                | 0.9996786376 |
| 0.0117170182_min                                | 1.26e-08     |
| 0.58204_min                                     | 0.39240      |
| 0.58204_max                                     | 0.70287      |
| 0.55164_min                                     | 0.34393      |
| 0.55164_max                                     | 0.72948      |
| bin(87)                                         | count        |
| 100                                             | 201          |
| 95                                              | 3303         |
| 90                                              | 5001         |
| 85                                              | 5093         |
| 80                                              | 4366         |
| 75                                              | 3736         |
| 70                                              | 3416         |
| 65                                              | 3362         |
| 60                                              | 3225         |
| 55                                              | 3070         |
| 50                                              | 2115         |
| 45                                              | 937          |
| 40                                              | 222          |
| 35                                              | 20           |
| 30                                              | 1            |

| #Item                                             | Value        |
|---------------------------------------------------|--------------|
| mci_bootstrap/mci.result.filter.tsv.bootstrap.tsv | 29606        |
| count                                             | 29605        |
| 0.3294813666_max                                  | 0.9999540193 |
| 0.3294813666_min                                  | 7.0081e-06   |
| 0.56688_min                                       | 0.36271      |
| 0.56688_max                                       | 0.70139      |
| 0.57210_min                                       | 0.37373      |
| 0.57210_max                                       | 0.68788      |
| bin(62)                                           | count        |
| 100                                               | 43           |
| 95                                                | 671          |
| 90                                                | 1117         |
| 85                                                | 1558         |
| 80                                                | 2177         |
| 75                                                | 2724         |
| 70                                                | 3281         |
| 65                                                | 3984         |
| 60                                                | 4726         |
| 55                                                | 4742         |
| 50                                                | 3203         |
| 45                                                | 1177         |
| 40                                                | 189          |
| 35                                                | 11           |
| 30                                                | 2            |

| #Item                                           | Value        |
|-------------------------------------------------|--------------|
| hc_bootstrap/hc.result.filter.tsv.bootstrap.tsv | 35392        |
| count                                           | 35391        |
| 0.0037077965_max                                | 0.9999244074 |
| 0.0037077965_min                                | 1.623e-07    |
| 0.58918_min                                     | 0.39913      |
| 0.58918_max                                     | 0.67022      |
| 0.55316_min                                     | 0.38489      |
| 0.55316_max                                     | 0.67991      |
| bin(87)                                         | count        |
| 100                                             | 243          |
| 95                                              | 3202         |
| 90                                              | 5894         |
| 85                                              | 5085         |
| 80                                              | 4037         |
| 75                                              | 3273         |
| 70                                              | 3010         |
| 65                                              | 2783         |
| 60                                              | 2716         |
| 55                                              | 2529         |
| 50                                              | 1745         |
| 45                                              | 695          |
| 40                                              | 160          |
| 35                                              | 18           |
| 30                                              | 1            |


```shell


```