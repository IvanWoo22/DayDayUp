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
tsv-join RRBS_sites_merge.tsv -f pw005.filter.tsv -k 1,2 >RRBS_sites_merge.filter.tsv
mkdir split
split RRBS_sites_merge.filter.tsv -l 1000 -d -a 4 split/loci_
cd split || exit
# shellcheck disable=SC2011
ls -- loci_* | xargs -i{} mv {} {}.tmp
cd ..
mkdir job
basename -s split/ -s .tmp -a split/*[0-9]* |
  sort |
  split -l 120 -a 2 -d - job/

for f in $(find job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 120 '
        echo '\''==> Processing {}'\''
        Rscript WelcoxonFilterSteep.R RRBS_sites_merge.filter.tsv split/{}.tmp split/{}.steep.tsv
      ' <${f}
    "
done

tsv-append split/loci_*.steep.tsv >steep.filter.tsv
tsv-join pw005.filter.tsv -f steep.filter.tsv -k 1,2 >pw005.filter.steep.tsv
tsv-join pw005.filter.Moran.tsv -f steep.filter.tsv -k 1,2 >pw005.filter.steep.Moran.tsv
tsv-join pw005.shuffle.Moran.tsv -f steep.filter.tsv -k 1,2 >pw005.shuffle.steep.Moran.tsv
rm output.*
rm -fr ./split ./*job

# On workstation.
Rscript MoranDraw.R pw005.filter.steep.Moran.tsv pw005.shuffle.steep.Moran.tsv pw005.filter.steep.Moran.pdf

tsv-join -e pw005.filter.tsv -f steep.filter.tsv -k 1,2 >pw005.filter.step2.tsv
tsv-join pw005.filter.step2.tsv -f kd001.tsv -k 1,2 >pw005.kd001.filter.step2.tsv
tsv-join pw005.filter.Moran.tsv -f pw005.kd001.filter.step2.tsv -k 1,2 > pw005.kd001.filter.step2.Moran.tsv
tsv-join RRBS_sites_merge.tsv -f pw005.kd001.filter.step2.tsv -k 1,2 >RRBS_sites_merge.filter.tsv

mkdir split
split pw005.kd001.filter.step2.tsv -l 100 -d -a 3 split/kd001
cd split || exit
# shellcheck disable=SC2011
ls -- kd001* | xargs -i{} mv {} {}.tmp
cd ..
mkdir -p kd001_job
basename -s split/kd001 -s .tmp -a split/kd001* |
  sort |
  split -l 128 -a 1 -d - kd001_job/
for f in $(find kd001_job -maxdepth 1 -type f -name "[0-9]*" | sort); do
  echo "${f}"
  bsub -n 128 -q amd_milan -J "-${f}" \
    "
      parallel --no-run-if-empty --line-buffer -k -j 128 '
        echo '\''==> Processing {}'\''
        Rscript WelcoxonFilterSmooth.R RRBS_sites_merge.filter.tsv split/{}.tmp split/{}.smooth.tsv
      ' <${f}
    "
done

tsv-append split/kd001*.smooth.tsv >pw005.kd001.filter.smooth.tsv
tsv-join pw005.filter.Moran.tsv -f pw005.kd001.filter.smooth.tsv -k 1,2 >pw005.filter.smooth.Moran.tsv
tsv-join pw005.shuffle.Moran.tsv -f pw005.kd001.filter.smooth.tsv -k 1,2 >pw005.shuffle.smooth.Moran.tsv

tsv-join RRBS_sites_merge.tsv -f pw005.kd001.filter.smooth.tsv -k 1,2 >changed_sites.f2.tsv
tsv-join RRBS_sites_merge.tsv -f pw005.filter.steep.tsv -k 1,2 >changed_sites.f4.tsv

rm -rf png_kendall
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
  tsv-append split/loci_${i}_*.kd.tsv >changed_sites.${i}.kd.tsv
  awk '$NF>0' changed_sites.${i}.kd.tsv | wc -l
  awk '$NF<0' changed_sites.${i}.kd.tsv | wc -l
done

####

sort -k1,1 -k2,2n changed_sites.f{2,4}.kd.tsv >changed_sites.kd.tsv
tsv-join RRBS_sites_merge.tsv -f changed_sites.kd.tsv -k 1,2 >RRBS_sites_merge.changed.tsv
awk '{print $1 "\t" $2 "\t" ($2+1) "\t" $NF "\t" 2}' changed_sites.f2.kd.tsv >changed_sites.f2.kd.bed
awk '{print $1 "\t" $2 "\t" ($2+1) "\t" $NF "\t" 1}' changed_sites.f4.kd.tsv >changed_sites.f4.kd.bed
sort -k1,1 -k2,2n changed_sites.f{2,4}.kd.bed >changed_sites.bed
cut -f 1 changed_sites.f2.kd.tsv | sort | uniq >job
parallel --no-run-if-empty --line-buffer -k -j 24 '
    echo '\''==> Processing {}'\''
    /usr/bin/Rscript Sites2RegionsInChanging.R changed_sites.bed changed_region.{}.bed {} 250 100000 3
    ' <job

sort -k1,1 -k2,2n changed_region.chr{{1..22},X}.bed |
  awk '$6*$7==0' >changed_regions.bed
rm changed_region.chr{{1..22},X}.bed job

awk '$8==0||($8==1&&$9>8)' changed_regions.bed >changed_regions.smooth.bed
awk '$9==0||($9==1&&$8>8)' changed_regions.bed >changed_regions.steep.bed
awk '!(($9==0||($9==1&&$8>8))||($8==0||($8==1&&$9>8)))' changed_regions.bed >changed_regions.mix.bed
rm changed_regions.*.{1..4}.beta.tsv

for fea in steep smooth mix; do
  awk '$7==0' changed_regions.${fea}.bed >changed_regions.${fea}.up.bed
done
for fea in steep smooth mix; do
  awk '$6==0' changed_regions.${fea}.bed >changed_regions.${fea}.down.bed
done

for dir in up down; do
  for fea in steep smooth mix; do
    closestBed -d -a changed_regions.${fea}.${dir}.bed -b Promoter_gencode.V40.bed |
      awk '$NF==0{print $13}' |
      sort | uniq | awk -F '.' '{print $1}' \
      >promoter/${fea}.${dir}.50.genelist.tsv
    closestBed -d -a changed_regions.${fea}.${dir}.bed -b Genebody_gencode.V40.bed |
      awk '$NF==0{print $13}' |
      sort | uniq | awk -F '.' '{print $1}' \
      >genebody/${fea}.${dir}.50.genelist.tsv
  done
done

wc -l promoter/* genebody/*
#  624 promoter/mix.down.50.genelist.tsv
#  170 promoter/mix.up.50.genelist.tsv
#   86 promoter/smooth.down.50.genelist.tsv
#   37 promoter/smooth.up.50.genelist.tsv
#  226 promoter/steep.down.50.genelist.tsv
#  106 promoter/steep.up.50.genelist.tsv
# 1281 genebody/mix.down.50.genelist.tsv
#  412 genebody/mix.up.50.genelist.tsv
#  181 genebody/smooth.down.50.genelist.tsv
#   78 genebody/smooth.up.50.genelist.tsv
#  516 genebody/steep.down.50.genelist.tsv
#  258 genebody/steep.up.50.genelist.tsv