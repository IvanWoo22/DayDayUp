




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

