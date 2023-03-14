# Step1.
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

tsv-append split/*.pw.tsv | sed 's/\"//g' |
  awk '
    {{printf $1 "\t" $2 "\t" $3 "\t" $4}; sum=0;
    {for(i=5;i<=NF;i++){if($i<0.05){printf "\t1"; sum++};if($i>=0.05){printf "\t0"}}};
    {printf "\t" sum "\n"}}
  ' >pw005.tsv

awk '{h[sprintf("%.0f",$NF)]++}END{for (i in h){print i,h[i]}}' pw005.tsv |
  sort -nrk 1,1
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

awk '$NF>1{for(i=1;i<(NF-1);i++){printf $i "\t"};printf $(NF-1) "\n"}' pw005.tsv \
  >pw005.filter.tsv
wc -l pw005.filter.tsv

#613724 pw005.filter.tsv

perl shuffle.pl <pw005.filter.tsv >pw005.shuffle.tsv
echo "=> 1 done."
for k in {2..10}; do
  perl shuffle.pl <pw005.shuffle.tsv \
    >pw005.shuffle.tsv.bak
  mv pw005.shuffle.tsv.bak pw005.shuffle.tsv
  echo "=> ${k} done."
done

mkdir split
split pw005.shuffle.tsv -l 2000 -d -a 3 split/pw005.shuffle
cd split || exit
# shellcheck disable=SC2011
ls -- pw005.shuffle* | xargs -i{} mv {} {}.tmp
cd ..
mkdir -p 005_job
basename -s split/pw005.shuffle -s .tmp -a split/pw005.shuffle* |
  sort |
  split -l 128 -a 1 -d - 005_job/

for f in $(find 005_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "moran-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript MoranI.R split/{}.tmp split/{}.Moran.tsv
    ' <${f}
    "
done

for f in $(find 005_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "moran-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript MoranI.R split/{}.tmp split/{}.Moran.tsv
    ' <${f}
    "
done
tsv-append split/pw005.shuffle*.Moran.tsv >pw005.shuffle.Moran.tsv

split pw005.filter.tsv -l 2000 -d -a 3 split/pw005.filter
cd split || exit
# shellcheck disable=SC2011
ls -- pw005.filter* | xargs -i{} mv {} {}.tmp
cd ..

mkdir -p 005_job
basename -s split/pw005.filter -s .tmp -a split/pw005.filter* |
  sort |
  split -l 128 -a 1 -d - 005_job/

for f in $(find 005_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "moran-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript MoranI.R split/{}.tmp split/{}.Moran.tsv
    ' <${f}
    "
done

tsv-append split/pw005.filter*.Moran.tsv >pw005.filter.Moran.tsv

rm output.*
rm -fr ./split ./*job

# On workstation.
Rscript MoranDraw.R pw005.filter.Moran.tsv pw005.shuffle.Moran.tsv pw005.Moran.pdf

# Step 2.
tsv-join RRBS_sites_merge.tsv -f pw005.filter.tsv -k 1,2 >RRBS_sites_merge.filter1.tsv

mkdir split
split RRBS_sites_merge.filter1.tsv -l 1200 -d -a 3 split/loci_
cd split || exit
# shellcheck disable=SC2011
ls -- loci_* | xargs -i{} mv {} {}.tmp
cd ..

mkdir job
basename -s split/ -s .tmp -a split/*[0-9]* |
  sort |
  split -l 120 -a 1 -d - job/

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 120 '
        echo '\''==> Processing {}'\''
        Rscript WelcoxonFilterOut.R split/{}.tmp split/{}.filter.out.tsv
      ' <${f}
    "
done

tsv-append split/loci_*.filter.out.tsv >filter.out.tsv

rm output.*
rm -fr ./split ./*job

tsv-join -e pw005.filter.tsv -f filter.out.tsv -k 1,2 >pw005.filter.out.tsv
tsv-join -e pw005.filter.Moran.tsv -f filter.out.tsv -k 1,2 >pw005.filter.out.Moran.tsv
tsv-join -e pw005.shuffle.Moran.tsv -f filter.out.tsv -k 1,2 >pw005.shuffle.out.Moran.tsv

wc -l pw005.filter.out.tsv

# On workstation.
Rscript MoranDraw.R pw005.filter.out.Moran.tsv pw005.shuffle.out.Moran.tsv pw005.filter.out.Moran.pdf

