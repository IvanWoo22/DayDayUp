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

tsv-join GDE.panther.replace.tsv \
  -f GDE.tigrfam.replace.tsv \
  >GDE.replace.tsv

wc -l ./*.replace.tsv
#  3974 ./GDE.replace.tsv
#  4032 ./panther.replace.tsv
#  3979 ./tigrfam.replace.tsv
# 11985 total

faops some ${PREFIX}/PROTEINS/all.replace.fa \
  <(cut -f 2 GDE.replace.tsv) \
  GDE.replace.fa

cd hmm || return
aria2c -x 12 https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM.tgz
cat hmm_PGAP/*.HMM >PGAP.hmm

cd ..

for DB in PFAM TIGR; do
  hmmpress hmm/${DB}/HMM
done

for DB in PANTHER PFAM TIGR; do
  bsub -n 24 -J "${DB}" "
  E_VALUE=\"1e-10\"
  NAME=\"GDE\"
  hmmscan --cpu 20 -E \${E_VALUE} --domE \${E_VALUE} --noali \
    -o \${NAME}.${DB}.progress.txt --tblout \${NAME}.${DB}.tbl \
    hmm/${DB}/HMM \${NAME}.replace.fa
    "
done

for DB in PANTHER PFAM TIGR; do
  NAME="GDE"
  grep -v "#" ${NAME}.${DB}.tbl |
    awk '{printf $1 "\t" $2 "\t" $3 "\t" $5 "\t" $6 "\t" $8 "\t" $9 "\t" $19;
    for(i=20;i<=NF;i++){gsub(/[[:blank:]]/, "", $i); printf "_" $i};printf "\n"}' \
      >${NAME}.abstract.${DB}.tsv
done

for DB in PANTHER PFAM TIGR; do
  NAME="GDE"
  cut -f 8 ${NAME}.abstract.${DB}.tsv | sort | uniq -c
  echo
done
#   7380 Amidohydrolase_family
#
#     59 allantoinase:_allantoinase
#      8 D-hydantoinase:_dihydropyrimidinase
#     12 EF_0837:_putative_amidohydrolase,_EF_0837/AHA_3915_family
#   3973 guan_deamin:_guanine_deaminase
#   3948 hutF:_formiminoglutamate_deiminase
#   2676 hutI:_imidazolonepropionase
#      1 isoAsp_dipep:_beta-aspartyl_peptidase
#     21 nagA:_N-acetylglucosamine-6-phosphate_deacetylase
#     10 one_C_dehyd_A:_formylmethanofuran_dehydrogenase_subunit_A
#     33 phosphono_phnM:_phosphonate_metabolism_protein_PhnM
#   3968 Se_ssnA:_putative_selenium_metabolism_protein_SsnA

for DB in PANTHER PFAM TIGR; do
  NAME="GDE"
  perl compare1.pl ${NAME}.abstract.${DB}.tsv \
    >${NAME}.minie.${DB}.tsv
  echo
done

cut -f 1,3 GDE.minie.PFAM.tsv | tsv-filter --str-in-fld 2:"Amidohydrolase_family" >GDE.kw_filter.PFAM.tsv
cut -f 1,3 GDE.minie.TIGR.tsv | tsv-filter --str-in-fld 2:"guan_deamin:_guanine_deaminase" >GDE.kw_filter.TIGR.tsv

for DB in PANTHER PFAM TIGR; do
  NAME="GDE"
  PREFIX=/share/home/wangq/data/Pseudomonas
  tsv-join -d 1 \
    -f ${PREFIX}/PROTEINS/all.strain.tsv -k 1 \
    --append-fields 2 \
    <${NAME}.kw_filter.${DB}.tsv |
    tsv-join -d 3 \
      -f strains.taxon.tsv -k 1 \
      --append-fields 4 |
    tsv-summarize -g 3,4 --count |
    keep-header -- tsv-sort -k3,3n >${NAME}.copy.${DB}.tsv
done

for DB in PANTHER PFAM TIGR; do
  NAME="GDE"
  tsv-summarize -g 3,2 --count ${NAME}.copy.${DB}.tsv \
    >${NAME}.GCFcopy.${DB}.tsv
  sed -i '1icopy\tgenus\tGCF' ${NAME}.GCFcopy.${DB}.tsv
done

cut -f 1 GDE.kw_filter.*.tsv | sort | uniq -c | awk '$1==2' | wc -l
#3472

mkdir -p blastp
PREFIX=/share/home/wangq/data/Pseudomonas
makeblastdb -in ${PREFIX}/PROTEINS/all.replace.fa -dbtype prot -out blastp/all_protein

tsv-join YggL.panther.replace.tsv \
  -f YggL.tigrfam.replace.tsv \
  >YggL.replace.tsv
faops some ${PREFIX}/PROTEINS/all.replace.fa \
  <(cut -f 2 YggL.replace.tsv) \
  YggL.replace.fa

for NAME in GDE YggL; do
  bsub -n 24 -q largemem -J "${NAME}" "
    blastp -db blastp/all_protein \
      -outfmt 6 -evalue 1e-5 -num_threads 16 \
      -query ${NAME}.replace.fa \
      -out blastp/${NAME}.result1.tsv"
done

for NAME in GDE YggL; do
  PREFIX=/share/home/wangq/data/Pseudomonas
  faops some ${PREFIX}/PROTEINS/all.replace.fa \
    <(cut -f 2 blastp/${NAME}.result1.tsv | sort | uniq) \
    ${NAME}.blastp1.fa
done

for NAME in GDE YggL; do
  bsub -n 24 -q largemem -J "${NAME}" "
    blastp -db blastp/all_protein \
      -outfmt 6 -evalue 1e-5 -num_threads 16 \
      -query ${NAME}.blastp1.fa \
      -out blastp/${NAME}.result2.tsv"
done

for NAME in GDE YggL; do
  PREFIX=/share/home/wangq/data/Pseudomonas
  faops some ${PREFIX}/PROTEINS/all.replace.fa \
    <(cut -f 2 blastp/${NAME}.result2.tsv | sort | uniq) \
    ${NAME}.blastp2.fa
done

for NAME in GDE YggL; do
  bsub -n 24 -q largemem -J "${NAME}" "
    blastp -db blastp/all_protein \
      -outfmt 6 -evalue 1e-5 -num_threads 16 \
      -query ${NAME}.blastp2.fa \
      -out blastp/${NAME}.result3.tsv"
done

GDE.replace.fa
