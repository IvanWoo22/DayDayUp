echo "Pseudom_aeru_PAO1
Pseudom_puti_KT2440_GCF_000007565_2
Pseudom_chl_aureofaciens_30_84_GCF_000281915_1
Pseudom_entomophi_L48_GCF_000026105_1
Pseudom_fluo_SBW25_GCF_000009225_2
Pseudom_prot_Pf_5_GCF_000012265_1
Pseudom_sav_pv_phaseolicola_1448A_GCF_000012205_1
Pseudom_stu_A1501_GCF_000013785_1
Pseudom_syr_pv_syringae_B728a_GCF_000012245_1
Pseudom_aeru_UCBPP_PA14_GCF_000014625_1
Pseudom_aeru_PA7_GCF_000017205_1
Pseudom_aeru_LESB58_GCF_000026645_1
" >typical.lst

PREFIX=/share/home/wangq/data/Pseudomonas

# Group by order
sed -e '1d' ${PREFIX}/ASSEMBLY/Pseudomonas.assembly.pass.csv |
  cut -d , -f 3 |
  tsv-uniq |
  nwr append stdin -r order |
  cut -f 2 |
  tsv-uniq \
    >${PREFIX}/order.lst

grep -h "IPR014311" ${PREFIX}/STRAINS/Pseudom_aeru_PAO1/*.fa.tsv | cut -f 1,4,5

NP_248824.1 CDD cd01303
NP_248824.1 PANTHER PTHR11271:SF6
NP_248824.1 TIGRFAM TIGR02967
NP_250212.1 CDD cd01303
NP_250212.1 TIGRFAM TIGR02967
NP_250212.1 PANTHER PTHR11271:SF6

NP_248824.1 GDEase00
NP_250212.1 GDEase02

NAME="GDE"

mkdir -p wyf/ATHGT
mkdir -p wyf/ATHGT/hmm
cd wyf/ATHGT || return
curl -L http://www.pantherdb.org/panther/exportHmm.jsp\?acc\=PTHR11271:SF6 >hmm/${NAME}.panther.hmm
aria2c -x 12 https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR02967.1.HMM -o hmm/${NAME}.tigrfam.hmm

E_VALUE="1e-20"
for DB in panther tigrfam; do
  while IFS= read -r ORDER; do
    echo >&2 "==> ORDER [${ORDER}]"
    parallel --no-run-if-empty --linebuffer -k -j 8 "
      gzip -dcf ${PREFIX}/ASSEMBLY/{}/*_protein.faa.gz |
        hmmsearch -E ${E_VALUE} --domE ${E_VALUE} --noali --notextw hmm/${NAME}.${DB}.hmm - |
        grep '>>' | perl -nl -e '
          m{>>\s+(\S+)} or next;
          \$n = \$1;
          \$s = \$n;
          \$s =~ s/\.\d+//;
          printf qq{%s\t%s_%s\n}, \$n, {}, \$s;'
      " <${PREFIX}/taxon/"${ORDER}"
  done <${PREFIX}/order.lst >${DB}.replace.tsv
  echo >&2
done

tsv-join panther.replace.tsv \
  -f tigrfam.replace.tsv \
  >GDE.replace.tsv

wc -l ./*.replace.tsv
#4032 panther.replace.tsv
#3979 tigrfam.replace.tsv
#8011 total

faops some ${PREFIX}/PROTEINS/all.replace.fa \
  <(cut -f 2 GDE.replace.tsv) \
  GDE.replace.fa

cd hmm || return
aria2c -x 12 https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz
cat hmm_PGAP/*.HMM >PGAP.hmm

cd ..

for DB in PFAM PGAP TIGR; do
  hmmpress hmm/${DB}/HMM
done

for DB in PANTHER PFAM PGAP TIGR; do
  bsub -n 24 -J "${DB}" "
  E_VALUE=\"1e-10\"
  NAME=\"GDE\"
  hmmscan --cpu 20 -E \${E_VALUE} --domE \${E_VALUE} --noali \
    -o \${NAME}.${DB}.progress.txt --tblout \${NAME}.${DB}.tbl \
    hmm/${DB}/HMM \${NAME}.replace.fa
    "
done

for DB in PANTHER PFAM PGAP; do
  NAME="GDE"
  grep -v "#" ${NAME}.${DB}.tbl |
    awk '{printf $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $8 "\t" $9 "\t" $19;
    for(i=20;i<=NF;i++){gsub(/[[:blank:]]/, "", $i); printf "_" $i};printf "\n"}' \
      >${NAME}.abstract.${DB}.tsv
done

tsv-filter --le 4:1e-50 ${NAME}.abstract.${DB}.tsv \
  >${NAME}.cutoff.${DB}.tsv

cut -f 8 ${NAME}.abstract.${DB}.tsv | sort | uniq -c
#Amidohydrolase_family	7380

perl compare1.pl branched-chain/branched-chain.cutoff.pfam.tsv \
  >branched-chain/branched-chain_minevalue.pfam.tsv

cat branched-chain/branched-chain_minevalue.pfam.tsv | tsv-select -f 1,3 |
  tsv-filter --str-in-fld 2:"Branched-chain" |
  tsv-join -d 1 \
    -f PROTEINS/all.strain.tsv -k 1 \
    --append-fields 2 |
  tsv-join -d 3 \
    -f strains.taxon.tsv -k 1 \
    --append-fields 4 |
  tsv-summarize -g 3,4 --count |
  keep-header -- tsv-sort -k3,3n >branched-chain/branched-chain_hmmscan_copy.pfam.tsv

#统计拷贝数的分布390
tsv-summarize -g 3,2 --count branched-chain/branched-chain_hmmscan_copy.pfam.tsv \
  >branched-chain/branched-chain_hmmscan_GCF_copy.pfam.tsv
sed -i '1icopy\tgenus\tGCF' branched-chain/branched-chain_hmmscan_GCF_copy.pfam.tsv

#查看pfam和tigerfam结果的交集2140
cat branched-chain/branched-chain_minevalue.pfam.tsv |
  grep -f <(cut -f 1 branched-chain/branched-chain_tigerfam.minevalue.tsv) | wc -l
