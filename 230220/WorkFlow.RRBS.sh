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

for i in 005 001 fdr005; do
  split pw${i}.kd005.shuffle.tsv -l 2000 -d -a 3 split/pw${i}.kd005.shuffle
  cd split || exit
  # shellcheck disable=SC2011
  ls -- pw${i}.kd005.shuffle* | xargs -i{} mv {} {}.tmp
  cd ..

  mkdir -p ${i}_job
  basename -s split/pw${i}.kd005.shuffle -s .tmp -a split/pw${i}.kd005.shuffle* |
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
  tsv-append split/pw${i}.kd005.shuffle*.Moran.tsv >pw${i}.kd005.shuffle.Moran.tsv
done

rm output.*
rm -fr ./split ./*job

###

mkdir split

for i in 005 001 fdr005; do
  split pw${i}.kd005.filter.tsv -l 2000 -d -a 3 split/pw${i}.kd005.filter
  cd split || exit
  # shellcheck disable=SC2011
  ls -- pw${i}.kd005.filter* | xargs -i{} mv {} {}.tmp
  cd ..

  mkdir -p ${i}_job
  basename -s split/pw${i}.kd005.filter -s .tmp -a split/pw${i}.kd005.filter* |
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
  tsv-append split/pw${i}.kd005.filter*.Moran.tsv >pw${i}.kd005.filter.Moran.tsv
done

rm output.*
rm -fr ./split ./*job

###
