mkdir split
split RRBS_sites_merge.tsv -l 2500 -d -a 3 split/loci_
cd split || exit
# shellcheck disable=SC2011
ls -- loci_* | xargs -i{} mv {} {}.tmp
cd ..

mkdir -p job
basename -s split/ -s .tmp -a split/*[0-9]* |
  sort |
  split -l 128 -a 1 -d - job/

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "ttest${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript PairwiseWilcoxon.R split/{}.tmp split/{}.pw.tsv
    ' <${f}
    "
done

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "ttest${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript PairwiseWilcoxonFDR.R split/{}.tmp split/{}.pwfdr.tsv
    ' <${f}
    "
done

tsv-append split/*.pw.tsv | sed 's/\"//g' |
  awk '
    {{printf $1 "\t" $2 "\t" $3 "\t" $4}; sum=0;
    {for(i=5;i<=NF;i++){if($i<0.05){printf "\t1"; sum++};if($i>=0.05){printf "\t0"}}};
    {printf "\t" sum "\n"}}
  ' >pw005.tsv
tsv-append split/*.pw.tsv | sed 's/\"//g' |
  awk '
    {{printf $1 "\t" $2 "\t" $3 "\t" $4}; sum=0;
    {for(i=5;i<=NF;i++){if($i<0.01){printf "\t1"; sum++};if($i>=0.01){printf "\t0"}}};
    {printf "\t" sum "\n"}}
  ' >pw001.tsv
tsv-append split/*.pwfdr.tsv | sed 's/\"//g' |
  awk '
    {{printf $1 "\t" $2 "\t" $3 "\t" $4}; sum=0;
    {for(i=5;i<=NF;i++){if($i<0.05){printf "\t1"; sum++};if($i>=0.05){printf "\t0"}}};
    {printf "\t" sum "\n"}}
  ' >pwfdr005.tsv

for i in 005 001 fdr005; do
  echo "==>${i}"
  awk '{h[sprintf("%.0f",$NF)]++}END{for (i in h){print i,h[i]}}' pw${i}.tsv |
    sort -nrk 1,1
done
#==>005
#17 1
#16 7
#15 48
#14 207
#13 643
#12 1794
#11 4165
#10 9407
#9 12611
#8 18669
#7 27986
#6 44814
#5 67451
#4 98234
#3 138474
#2 189213
#1 246728
#0 670141
#==>001
#13 2
#12 14
#11 82
#10 576
#9 810
#8 1456
#7 2269
#6 4001
#5 6962
#4 12010
#3 22128
#2 46601
#1 122174
#0 1311508
#==>fdr005
#15 1
#14 4
#13 24
#12 108
#11 274
#10 2297
#9 2790
#8 2533
#7 3079
#6 4702
#5 5659
#4 6608
#3 8316
#2 12479
#1 19714
#0 1462005

for i in 005 001 fdr005; do
  awk '$NF>1{for(i=1;i<(NF-1);i++){printf $i "\t"};printf $(NF-1) "\n"}' pw${i}.tsv \
    >pw${i}.filter.tsv
  wc -l pw${i}.filter.tsv
done
#613724 pw005.filter.tsv
#96911 pw001.filter.tsv
#48874 pwfdr005.filter.tsv

for i in 005 001 fdr005; do
  echo "==> Processing ${i} ..."
  perl shuffle.pl <pw${i}.filter.tsv >pw${i}.shuffle.tsv
  echo "=> 1 done."
  # shellcheck disable=SC2034
  for k in {2..10}; do
    perl shuffle.pl <pw${i}.shuffle.tsv \
      >pw${i}.shuffle.tsv.bak
    mv pw${i}.shuffle.tsv.bak pw${i}.shuffle.tsv
    echo "=> ${k} done."
  done
  echo
done

###

mkdir split

for i in 005 001 fdr005; do
  split pw${i}.shuffle.tsv -l 2000 -d -a 3 split/pw${i}.shuffle
  cd split || exit
  # shellcheck disable=SC2011
  ls -- pw${i}.shuffle* | xargs -i{} mv {} {}.tmp
  cd ..

  mkdir -p ${i}_job
  basename -s split/pw${i}.shuffle -s .tmp -a split/pw${i}.shuffle* |
    sort |
    split -l 128 -a 1 -d - ${i}_job/
done

for i in 005 001 fdr005; do
  for f in $(find ${i}_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo "${f}"
    bsub -n 128 -q amd_milan -J "moran-${f}" \
      "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript MoranI.R split/{}.tmp split/{}.Moran.tsv
    ' <${f}
    "
  done
done

for i in 005 001 fdr005; do
  tsv-append split/pw${i}.shuffle*.Moran.tsv >pw${i}.shuffle.Moran.tsv
done

rm output.*
rm -fr ./split ./*job

###

mkdir split

for i in 005 001 fdr005; do
  split pw${i}.filter.tsv -l 2000 -d -a 3 split/pw${i}.filter
  cd split || exit
  # shellcheck disable=SC2011
  ls -- pw${i}.filter* | xargs -i{} mv {} {}.tmp
  cd ..

  mkdir -p ${i}_job
  basename -s split/pw${i}.filter -s .tmp -a split/pw${i}.filter* |
    sort |
    split -l 128 -a 1 -d - ${i}_job/
done

for i in 005 001 fdr005; do
  for f in $(find ${i}_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo "${f}"
    bsub -n 128 -q amd_milan -J "moran-${f}" \
      "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript MoranI.R split/{}.tmp split/{}.Moran.tsv
    ' <${f}
    "
  done
done

for i in 005 001 fdr005; do
  tsv-append split/pw${i}.filter*.Moran.tsv >pw${i}.filter.Moran.tsv
done

rm output.*
rm -fr ./split ./*job

###

for i in 005 001 fdr005; do
  bsub -n 1 -J "Moran-${i}" \
    "
    Rscript MoranDraw.R pw${i}.filter.Moran.tsv pw${i}.shuffle.Moran.tsv pw${i}.Moran.pdf
  "
done

###

mkdir split
split RRBS_sites_merge.tsv -l 2500 -d -a 3 split/loci_
cd split || exit
# shellcheck disable=SC2011
ls -- loci_* | xargs -i{} mv {} {}.tmp
cd ..

mkdir -p job
basename -s split/ -s .tmp -a split/*[0-9]* |
  sort |
  split -l 128 -a 1 -d - job/

mkdir png_kendall
for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "kd${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
        source /share/home/wangq/miniconda3/etc/profile.d/conda.sh
        conda activate renv
        Rscript KendallAll.R split/{}.tmp split/{}.kd.tsv
    ' <${f}
    "
done

tsv-append split/loci_*.kd.tsv |
  awk '$(NF-1)<0.05' >kd005.tsv

for i in 005 001 fdr005; do
  tsv-join pw${i}.filter.tsv -f kd005.tsv -k 1,2 >pw${i}.kd005.filter.tsv
  wc -l pw${i}.kd005.filter.tsv
done

#135473 pw005.kd005.filter.tsv
#56611 pw001.kd005.filter.tsv
#37524 pwfdr005.kd005.filter.tsv

for i in 005 001 fdr005; do
  echo "==> Processing ${i} ..."
  perl shuffle.pl <pw${i}.kd005.filter.tsv >pw${i}.kd005.shuffle.tsv
  echo "=> 1 done."
  # shellcheck disable=SC2034
  for k in {2..10}; do
    perl shuffle.pl <pw${i}.kd005.shuffle.tsv \
      >pw${i}.kd005.shuffle.tsv.bak
    mv pw${i}.kd005.shuffle.tsv.bak pw${i}.kd005.shuffle.tsv
    echo "=> ${k} done."
  done
  echo
done

###

mkdir split

for j in filter shuffle; do
  for i in 005 001 fdr005; do
    split pw${i}.kd005.${j}.tsv -l 2000 -d -a 3 split/pw${i}.kd005.${j}
    cd split || exit
    # shellcheck disable=SC2011
    ls -- pw${i}.kd005.${j}* | xargs -i{} mv {} {}.tmp
    cd ..

    mkdir -p ${i}_job
    basename -s split/pw${i}.kd005.${j} -s .tmp -a split/pw${i}.kd005.${j}* |
      sort |
      split -l 128 -a 1 -d - ${i}_job/

    for f in $(find ${i}_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
      echo "${f}"
      bsub -n 128 -q amd_milan -J "moran-${f}" \
        "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript MoranI.R split/{}.tmp split/{}.Moran.tsv
    ' <${f}
    "
    done
  done
done
for j in filter shuffle; do
  for i in 005 001 fdr005; do
    tsv-append split/pw${i}.kd005.${j}*.Moran.tsv >pw${i}.kd005.${j}.Moran.tsv
  done
done

rm output.*
rm -fr ./split ./*job

###

cut -f 1 RRBS_sites_merge.tsv | sort | uniq >job

bsub -n 128 -q amd_milan -J "region-test" \
  "
  parallel --no-run-if-empty --line-buffer -k -j 128 '
    echo '\''==> Processing {}'\''
    Rscript RegionPairwiseWilcoxon.R RRBS_sites_merge.tsv {} {}.pw.tsv
    ' <job
  "

tsv-append chr*.pw.tsv |
  awk '
    {{printf $1 "\t" $2 "\t" $3 "\t" $4}; sum=0;
    {for(i=5;i<=NF;i++){if($i<0.05){printf "\t1"; sum++};if($i>=0.05){printf "\t0"}}};
    {printf "\t" sum "\n"}}
  ' >region.pw005.tsv

awk '{h[sprintf("%.0f",$NF)]++}END{for (i in h){print i,h[i]}}' region.pw005.tsv | sort -nrk 1,1
#21 1
#20 14
#19 65
#18 241
#17 579
#16 1186
#15 1844
#14 2804
#13 3937
#12 5365
#11 7308
#10 9938
#9 10747
#8 14156
#7 18923
#6 26547
#5 34458
#4 43527
#3 54533
#2 66609
#1 79696
#0 109643
awk '{h[sprintf("%.0f",($4/10))]++}END{for (i in h){print i,h[i]}}' region.pw005.tsv | sort -nrk 1,1
#90-100 4
#80-90 9
#70-80 8
#60-70 67
#50-60 160
#40-50 611
#30-40 1934
#20-30 10045
#10-20 53617
#0-10 425666
awk '$NF>1&&$4>2{for(i=1;i<(NF-1);i++){printf $i "\t"};printf $(NF-1) "\n"}' region.pw005.tsv >region.pw005.filter.tsv
#125718 region.pw005.filter.tsv

perl shuffle.pl <region.pw005.filter.tsv >region.pw005.shuffle.tsv
echo "=> 1 done."
# shellcheck disable=SC2034
for k in {2..10}; do
  perl shuffle.pl <region.pw005.shuffle.tsv \
    >region.pw005.shuffle.tsv.bak
  mv region.pw005.shuffle.tsv.bak region.pw005.shuffle.tsv
  echo "=> ${k} done."
done

###

mkdir split

for i in filter shuffle; do
  split region.pw005.${i}.tsv -l 600 -d -a 3 split/region.pw005.${i}
  cd split || exit
  # shellcheck disable=SC2011
  ls -- region.pw005.${i}* | xargs -i{} mv {} {}.tmp
  cd ..

  mkdir -p ${i}_job
  basename -s split/region.pw005.${i} -s .tmp -a split/region.pw005.${i}* |
    sort |
    split -l 128 -a 1 -d - ${i}_job/
done

for i in filter shuffle; do
  for f in $(find ${i}_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
    echo "${f}"
    bsub -n 128 -q amd_milan -J "moran-${f}" \
      "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript MoranI.R split/{}.tmp split/{}.Moran.tsv
    ' <${f}
    "
  done
done

for i in filter shuffle; do
  tsv-append split/region.pw005.${i}*.Moran.tsv >region.pw005.${i}.Moran.tsv
done

rm output.*
rm -fr ./split ./*job

###

cut -f 1 RRBS_sites_merge.tsv | sort | uniq >job
bsub -n 128 -q amd_milan -J "region-kd" "
  parallel --no-run-if-empty --line-buffer -k -j 128 '
    echo '\''==> Processing {}'\''
    source /share/home/wangq/miniconda3/etc/profile.d/conda.sh
    conda activate renv
    Rscript RegionKendallAll.R RRBS_sites_merge.tsv {} {}.kd.tsv
    ' <job
  "

###

tsv-join RRBS_sites_merge.tsv -f kd005.tsv -k 1,2 >RRBS_sites_merge.filter.tsv
mkdir split
for i in 001 005; do
  split kd${i}.tsv -l 450 -d -a 3 split/kd${i}
  cd split || exit
  # shellcheck disable=SC2011
  ls -- kd${i}* | xargs -i{} mv {} {}.tmp
  cd ..

  mkdir -p ${i}_job
  basename -s split/kd${i} -s .tmp -a split/kd${i}* |
    sort |
    split -l 128 -a 1 -d - ${i}_job/
done

for f in $(find 005_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
        echo '\''==> Processing {}'\''
        Rscript WelcoxonSmoothSitesTest2.R RRBS_sites_merge.filter.tsv split/{}.tmp split/{}.smooth.test2.tsv
        Rscript WelcoxonFilterTest4.R RRBS_sites_merge.filter.tsv split/{}.tmp split/{}.unusual.test4.tsv
        Rscript WelcoxonFilterTest5.R RRBS_sites_merge.filter.tsv split/{}.tmp split/{}.unusual.test5.tsv
      ' <${f}
    "
done

i=005
tsv-append split/kd${i}*.smooth.test2.tsv >kd${i}.filter.smooth.test2.tsv
tsv-append split/kd${i}*.unusual.test5.tsv >kd${i}.filter.unusual.test5.tsv
tsv-join kd${i}.filter.smooth.test2.tsv -f pw005.filter.tsv -k 1,2 >pw005.kd${i}.filter.smooth.test2.tsv
tsv-join -e pw005.filter.tsv -f kd${i}.filter.unusual.test5.tsv -k 1,2 >pw005.kd${i}.filter.unusual.test5.tsv
i=001
tsv-append split/kd${i}*.unusual.test4.tsv >kd${i}.filter.unusual.test4.tsv
tsv-join -e pw005.filter.tsv -f kd${i}.filter.unusual.test4.tsv -k 1,2 >pw005.kd${i}.filter.unusual.test4.tsv

rm output.*
rm -fr ./split ./*job

tsv-join RRBS_sites_merge.tsv -f pw005.kd005.filter.smooth.test2.tsv -k 1,2 >changed_sites.f1.tsv
tsv-join RRBS_sites_merge.tsv -f pw005.kd001.filter.unusual.test4.tsv -k 1,2 >changed_sites.f2.tsv
tsv-join RRBS_sites_merge.tsv -f pw005.kd005.filter.unusual.test5.tsv -k 1,2 >changed_sites.f3.tsv
tsv-join RRBS_sites_merge.tsv -f pw005.kd001.filter.unusual.steep.tsv -k 1,2 >changed_sites.f4.tsv

for i in f1 f2 f3; do
  bsub -n 1 -J "Moran-${i}" \
    "
    source /share/home/wangq/miniconda3/etc/profile.d/conda.sh
    conda activate renv
    tsv-join pw005.filter.Moran.tsv -f changed_sites.${i}.tsv -k 1,2 >${i}.filter.Moran.tsv
    tsv-join pw005.shuffle.Moran.tsv -f changed_sites.${i}.tsv -k 1,2 >${i}.shuffle.Moran.tsv
    Rscript MoranDraw.R ${i}.filter.Moran.tsv ${i}.shuffle.Moran.tsv changed_sites.${i}.Moran.pdf
    rm ${i}.filter.Moran.tsv ${i}.shuffle.Moran.tsv
  "
done

###

tsv-join kd005.tsv -f pw005.kd005.filter.tsv -k 1,2 >pw005.kd005.kendall.tsv
tsv-join kd005.tsv -f pw005.kd001.filter.tsv -k 1,2 >pw005.kd001.kendall.tsv

for i in 001 005; do
  sort -k1,1 -k2,2n <pw005.kd${i}.kendall.tsv |
    cut -f 1,2,31 |
    perl site_merge.pl 50 1000 2 |
    sort -k1,1 -k2,2n >changed.kd${i}.50.bed
  awk '($4>2&&$5==0)' changed.kd${i}.50.bed \
    >changed.up.kd${i}.50.bed
  awk '($4>2&&$6==0)' changed.kd${i}.50.bed \
    >changed.down.kd${i}.50.bed
done

wc -l changed.*.50.bed
#  3309 changed.down.kd001.50.bed
#  7344 changed.down.kd005.50.bed
#  4560 changed.kd001.50.bed
# 10338 changed.kd005.50.bed
#  1246 changed.up.kd001.50.bed
#  2834 changed.up.kd005.50.bed

cut -f 1 RRBS_sites_merge.tsv | sort | uniq >job
bsub -n 128 -q amd_milan -J "region-test" \
  "
  parallel --no-run-if-empty --line-buffer -k -j 128 '
    echo '\''==> Processing {}'\''
    Rscript Sites2Regions.R RRBS_sites_merge.tsv {} RRBS_regions.{}.bed RRBS_regioned_sites.{}.bed 50 1000 3
    ' <job
  "
sort -k1,1 -k2,2n RRBS_regions.chr{{1..22},X,Y,M}.bed >RRBS_regions_merge.bed
sort -k1,1 -k2,2n RRBS_regioned_sites.chr{{1..22},X,Y,M}.bed >RRBS_regioned_sites_merge.bed
rm RRBS_regions.chr{{1..22},X,Y,M}.bed
rm RRBS_regioned_sites.chr{{1..22},X,Y,M}.bed

closestBed -d -a RRBS_regions_merge.bed -b Promoter_gencode.V40.bed |
  awk '$NF==0{print $8}' |
  sort | uniq | awk -F '.' '{print $1}' \
  >promoter/background.genelist.tsv
closestBed -d -a RRBS_regions_merge.bed -b Genebody_gencode.V40.bed |
  awk '$NF==0{print $8}' |
  sort | uniq | awk -F '.' '{print $1}' \
  >genebody/background.genelist.tsv

mkdir -p split job png_kendall
for i in f2 f4; do
  split changed_sites.${i}.tsv -l 300 -d -a 3 split/loci_${i}_
  cd split || exit
  # shellcheck disable=SC2011
  ls -- loci_${i}_* | xargs -i{} mv {} {}.tmp
  cd ..
done

basename -s split/ -s .tmp -a split/*[0-9]* |
  sort |
  split -l 128 -a 1 -d - job/

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "kd${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
        source /share/home/wangq/miniconda3/etc/profile.d/conda.sh
        conda activate renv
        Rscript KendallAll.R split/{}.tmp split/{}.kd.tsv
    ' <${f}
    "
done

for i in f2 f4; do
  tsv-append split/loci_${i}_*.kd.tsv |
    awk '$(NF-1)<0.05' >changed_sites.${i}.kd.tsv
done
for i in f2 f4; do
  awk '$NF>0' changed_sites.${i}.kd.tsv | wc -l
  awk '$NF<0' changed_sites.${i}.kd.tsv | wc -l
done
#5902
#11596
#6785
#13261

for i in f2 f4; do
  awk '{print $1 "\t" $2 "\t" $2}' changed_sites.${i}.kd.tsv |
    sort -k1,1 -k2,2n >changed_sites.${i}.kd.bed
done
for i in up down; do
  closestBed -d -a changed.${i}.kd001.50.bed -b changed_sites.f2.kd.bed |
    awk '$NF==0' >changed_sites.${i}.f2.region50.bed
  closestBed -d -a changed.${i}.kd001.50.bed -b changed_sites.f4.kd.bed |
    awk '$NF==0' >changed_sites.${i}.f4.region50.bed
done

wc -l changed_sites.*.region*
#  4302 changed_sites.down.f2.region50.bed
#  1384 changed_sites.up.f2.region50.bed
#  7076 changed_sites.down.f4.region50.bed
#  2417 changed_sites.up.f4.region50.bed

mkdir changed_region
for i in f2 f4; do
  echo "==> ${i}"
  for k in up down; do
    cut -f 1-6 changed_sites.${k}.${i}.region50.bed |
      sort -k1,1 -k2,2n | uniq >changed_region/${k}.${i}.50.bed
    wc -l changed_region/${k}.${i}.50.bed
  done
done
#==> f2
#722 changed_region/up.f2.50.bed
#2177 changed_region/down.f2.50.bed
#==> f4
#974 changed_region/up.f4.50.bed
#2658 changed_region/down.f4.50.bed

mkdir promoter genebody
for i in f2 f4; do
  for k in up down; do
    closestBed -d -a changed_region/${k}.${i}.50.bed -b Promoter_gencode.V40.bed |
      awk '$NF==0{print $10}' |
      sort | uniq | awk -F '.' '{print $1}' \
      >promoter/${k}.${i}.50.genelist.tsv
    closestBed -d -a changed_region/${k}.${i}.50.bed -b Genebody_gencode.V40.bed |
      awk '$NF==0{print $10}' |
      sort | uniq | awk -F '.' '{print $1}' \
      >genebody/${k}.${i}.50.genelist.tsv
  done
done

wc -l promoter/* genebody/*
# 22509 promoter/background.genelist.tsv
#  1150 promoter/down.f1.50.genelist.tsv
#   733 promoter/down.f2.50.genelist.tsv
#   840 promoter/down.f4.50.genelist.tsv
#   428 promoter/up.f1.50.genelist.tsv
#   238 promoter/up.f2.50.genelist.tsv
#   316 promoter/up.f4.50.genelist.tsv
# 26198 genebody/background.genelist.tsv
#  2268 genebody/down.f1.50.genelist.tsv
#  1531 genebody/down.f2.50.genelist.tsv
#  1739 genebody/down.f4.50.genelist.tsv
#   954 genebody/up.f1.50.genelist.tsv
#   586 genebody/up.f2.50.genelist.tsv
#   719 genebody/up.f4.50.genelist.tsv

tsv-join RRBS_sites_merge.tsv -f kd001.tsv -k 1,2 >RRBS_sites_merge.filter.tsv
mkdir split
split kd001.tsv -l 200 -d -a 3 split/kd001
cd split || exit
# shellcheck disable=SC2011
ls -- kd001* | xargs -i{} mv {} {}.tmp
cd ..
mkdir -p 001_job
basename -s split/kd001 -s .tmp -a split/kd001* |
  sort |
  split -l 128 -a 1 -d - 001_job/

for f in $(find 001_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
        echo '\''==> Processing {}'\''
        Rscript WelcoxonFilterTestSteep1.R RRBS_sites_merge.filter.tsv split/{}.tmp split/{}.unusual.steep1.tsv
      ' <${f}
    "
done

tsv-append split/kd001*.unusual.steep1.tsv >kd001.filter.unusual.steep1.tsv
tsv-join pw005.filter.tsv -f kd001.filter.unusual.steep1.tsv -k 1,2 >pw005.kd001.filter.unusual.steep1.tsv

rm output.*
rm -fr ./split ./*job

tsv-join RRBS_sites_merge.tsv -f pw005.kd001.filter.unusual.steep.tsv -k 1,2 >changed_sites.steep.tsv

mkdir -p split job png_kendall
split changed_sites.steep.tsv -l 300 -d -a 3 split/loci_steep_
cd split || exit
# shellcheck disable=SC2011
ls -- loci_steep_* | xargs -i{} mv {} {}.tmp
cd ..

basename -s split/ -s .tmp -a split/*[0-9]* |
  sort |
  split -l 128 -a 1 -d - job/

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "kd${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
        source /share/home/wangq/miniconda3/etc/profile.d/conda.sh
        conda activate renv
        Rscript KendallAll.R split/{}.tmp split/{}.kd.tsv
    ' <${f}
    "
done

tsv-append split/loci_steep_*.kd.tsv |
  awk '$(NF-1)<0.05' >changed_sites.steep.kd.tsv
awk '$NF>0' changed_sites.steep.kd.tsv | wc -l
awk '$NF<0' changed_sites.steep.kd.tsv | wc -l

mkdir split
split kd001.tsv -l 100 -d -a 3 split/kd001
cd split || exit
# shellcheck disable=SC2011
ls -- kd001* | xargs -i{} mv {} {}.tmp
cd ..
mkdir -p 001_job
basename -s split/kd001 -s .tmp -a split/kd001* |
  sort |
  split -l 128 -a 1 -d - 001_job/

for f in $(find 001_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
        echo '\''==> Processing {}'\''
        Rscript WelcoxonFilterTestEdge3.R RRBS_sites_merge.filter.tsv split/{}.tmp split/{}.unusual.edge3.tsv
        Rscript WelcoxonFilterTestEdge4.R RRBS_sites_merge.filter.tsv split/{}.tmp split/{}.unusual.edge4.tsv
        Rscript WelcoxonFilterTestEdge5.R RRBS_sites_merge.filter.tsv split/{}.tmp split/{}.unusual.edge5.tsv
      ' <${f}
    "
done

for i in edge3 edge4 edge5; do
  tsv-append split/kd001*.unusual.${i}.tsv >kd001.filter.unusual.${i}.tsv
  tsv-join pw005.filter.tsv -f kd001.filter.unusual.${i}.tsv -k 1,2 >pw005.kd001.filter.unusual.${i}.tsv
done
