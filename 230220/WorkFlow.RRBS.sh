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
      Rscript PairwiseWilcoxon.R split/{}.tmp split/{}_pww.tsv
    ' <${f}
    "
done

tsv-append split/*_pww.tsv | sed 's/\"//g' |
  awk '
    {{printf $1 "\t" $2 "\t" $3 "\t" $4}; sum=0;
    {for(i=5;i<=NF;i++){if($i<0.05){printf "\t1"; sum++};if($i>=0.05){printf "\t0"}}};
    {printf "\t" sum "\n"}}
  ' >pwwfdr.tsv

awk '{h[sprintf("%.2f",$NF)]++}END{for (i in h){print i,h[i]}}' pww.tsv |
  sort -nrk 1,1
#17.00 1
#16.00 7
#15.00 48
#14.00 207
#13.00 643
#12.00 1794
#11.00 4165
#10.00 9407
#9.00 12611
#8.00 18669
#7.00 27986
#6.00 44814
#5.00 67451
#4.00 98234
#3.00 138474
#2.00 189213
#1.00 246728
#0.00 670141

awk '$NF>0{for(i=1;i<(NF-1);i++){printf $i "\t"};printf $(NF-1) "\n"}' pww.tsv \
  >pww_wc005.tsv

perl shuffle.pl <pww_wc005.tsv >pww_wc005.shuffle.tsv
# shellcheck disable=SC2034
for i in {1..9}; do
  perl shuffle.pl <pww_wc005.shuffle.tsv \
    >pww_wc005.shuffle.tsv.bak
  mv pww_wc005.shuffle.tsv.bak pww_wc005.shuffle.tsv
done

mkdir split
split pww_wc005.shuffle.tsv -l 2000 -d -a 3 split/loci_
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
  bsub -n 128 -q amd_milan -J "moran-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript MoranI.R split/{}.tmp split/{}_Moran.tsv
    ' <${f}
    "
done
tsv-append split/*_Moran.tsv >pww_wc005_Moran.shuffle.tsv

rm output.*
rm -fr ./split ./job

###

mkdir split
split pww_wc005.tsv -l 2000 -d -a 3 split/loci_
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
  bsub -n 128 -q amd_milan -J "moran-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
      echo '\''==> Processing {}'\''
      Rscript MoranI.R split/{}.tmp split/{}_Moran.tsv
    ' <${f}
    "
done
tsv-append split/*_Moran.tsv >pww_wc005_Moran.tsv

rm output.*
rm -fr ./split ./job
