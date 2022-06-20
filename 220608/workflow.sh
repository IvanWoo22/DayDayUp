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

grep -h "IPR004361" ${PREFIX}/STRAINS/Pseudom_aeru_PAO1/*.fa.tsv | cut -f 1,4,5

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

#1564
#4641

for NAME in GDE YggL; do
  PREFIX=/share/home/wangq/data/Pseudomonas
  faops some ${PREFIX}/PROTEINS/all.replace.fa \
    <(cut -f 2 blastp/${NAME}.result1.tsv | sort | uniq) \
    ${NAME}.blastp1.fa
  bsub -n 24 -q largemem -J "${NAME}" "
    blastp -db blastp/all_protein \
      -outfmt 6 -evalue 1e-5 -num_threads 16 \
      -query ${NAME}.blastp1.fa \
      -out blastp/${NAME}.result2.tsv"
done

#1565
#4641

for NAME in GDE YggL; do
  PREFIX=/share/home/wangq/data/Pseudomonas
  faops some ${PREFIX}/PROTEINS/all.replace.fa \
    <(cut -f 2 blastp/${NAME}.result2.tsv | sort | uniq) \
    ${NAME}.blastp2.fa
  bsub -n 24 -q largemem -J "${NAME}" "
    blastp -db blastp/all_protein \
      -outfmt 6 -evalue 1e-5 -num_threads 16 \
      -query ${NAME}.blastp2.fa \
      -out blastp/${NAME}.result3.tsv"
done

cut -f 1 blastp/YggL.result3.tsv | sort | uniq >YggL.blastp.tsv

cut -f 1 blastp/${NAME}.result3.tsv | sort | uniq -c | awk '$1==2' | wc -l

GDE.replace.fa

E_VALUE=1e-20

cd ~/data/Pseudomonas || return

# Find all genes
sed '1d' <~/Scripts/withncbi/hmm/bac120.tsv | cut -f 1 | cut -d . -f 1 |
  parallel --xapply -j 22 "
     mkdir -p PROTEINS/{}
     while IFS= read -r GENUS; do
       while IFS= read -r STRAIN; do
         gzip -dcf /share/home/wangq/data/Pseudomonas/ASSEMBLY/\${STRAIN}/*_protein.faa.gz |
           hmmsearch -E 1e-20 --domE 1e-20 \
             --noali --notextw hmm/HMM/{}.HMM - |
           grep '>>' |
           (STRAIN=\${STRAIN}; perl -nl -e '
                     />>\\s+(\\S+)/ and printf qq{%s\\t%s\\n}, \$1, \$ENV{STRAIN};
                 ')
       done <taxon/\${GENUS} \
         >PROTEINS/{}/\${GENUS}.replace.tsv
     done <genus.lst
     echo
 "

sed '1d' <~/Scripts/withncbi/hmm/bac120.tsv | cut -f 1 | cut -d . -f 1 |
  parallel --no-run-if-empty --linebuffer -k -j 8 "
        while IFS= read -r GENUS; do
          cat PROTEINS/{}/\${GENUS}.replace.tsv
        done <genus.lst >PROTEINS/{}/{}.replace.tsv
"

sed '1d' <~/Scripts/withncbi/hmm/bac120.tsv | cut -f 1 | cut -d . -f 1 |
  parallel --no-run-if-empty --linebuffer -k -j 8 '
        cat PROTEINS/{}/{}.replace.tsv |
            grep -f model.lst |
            grep -v "GCF" > PROTEINS/{}/{}.model.tsv
       faops some /share/home/wangq/data/Pseudomonas/PROTEINS/all.uniq.fa <(
            cat PROTEINS/{}/{}.model.tsv |
                cut -f 1 |
                tsv-uniq
            ) stdout > PROTEINS/{}/{}.model.fa
    '

sed '1d' <~/Scripts/withncbi/hmm/bac120.tsv | cut -f 1 | cut -d . -f 1 |
  parallel --no-run-if-empty --linebuffer -k -j 8 '
        >&2 echo "==> marker [{}]"
        mafft --auto PROTEINS/{}/{}.model.fa >PROTEINS/{}/{}.model.aln.fa
    '

sed '1d' <~/Scripts/withncbi/hmm/bac120.tsv | cut -f 1 | cut -d . -f 1 |
  while IFS= read -r marker; do
    echo >&2 "==> marker [${marker}]"
    parallel --no-run-if-empty --linebuffer -k -j 4 "
            faops replace -s PROTEINS/${marker}/${marker}.model.aln.fa <(echo {}) stdout
        " <PROTEINS/"${marker}"/"${marker}".model.tsv \
      >PROTEINS/"${marker}"/"${marker}".model.replace.fa
  done

sed '1d' <~/Scripts/withncbi/hmm/bac120.tsv | cut -f 1 | cut -d . -f 1 |
  while IFS= read -r marker; do
    faops filter -l 0 PROTEINS/"${marker}"/"${marker}".model.replace.fa stdout
    echo
  done >PROTEINS/bac120.model.aln.fas

fasops concat PROTEINS/bac120.model.aln.fas model/model.lst -o PROTEINS/bac120.model.aln.fa

trimal -in PROTEINS/bac120.model.aln.fa -out PROTEINS/bac120.model.trim.fa -automated1

sed -e "s/Cam_jej_jejuni_NCTC_11168_ATCC_700819/Cam_jej_jejuni_NCTC_11168/" PROTEINS/bac120.model.trim.fa >PROTEINS/bac120.model.rename.fa

bsub -q largemem -n 24 -J "iq" ../iqtree-2.2.0-Linux/bin/iqtree2 -s PROTEINS/bac120.model.rename.fa \
  --prefix branched-chain.Pseudom_aeru -T AUTO -b 100

mafft --auto PROTEINS/bac120.model.rename.fa >PROTEINS/bac120.model.mafft.fa

tsv-join -d 1 -f strains.taxon.tsv -k 1 --append-fields 4 <representative.tsv |
  tsv-select -f 2,1 |
  nwr append stdin -r class |
  tsv-filter --str-in-fld 3:"Gammaproteobacteria" |
  sort | uniq >Gammaproteobacteria.strain.tsv
grep -f <(cut -f 2 Gammaproteobacteria.strain.tsv) YggL/YggL.replace.cross.tsv |
  cut -f 2 >YggL.Gammaproteobacteria.protein.tsv
faops some ${PREFIX}/PROTEINS/all.replace.fa \
  YggL.Gammaproteobacteria.protein.fa <YggL.Gammaproteobacteria.protein.tsv

mafft --auto YggL.Gammaproteobacteria.protein.fa \
  >YggL.Gammaproteobacteria.protein.mafft.fa

bsub -n 24 -J "iq8" ../iqtree-2.2.0-Linux/bin/iqtree2 -s PROTEINS/YggL.Gammaproteobacteria.protein.mafft.fa --prefix test8 -T AUTO -B 1000 -bnni -m Q.pfam+I+I+R5

perl proteinid_extract_taxoninfo.pl \
  strains.taxon.tsv \
  <YggL.Gammaproteobacteria.protein.tsv |
  sort | uniq >YggL.Gammaproteobacteria.taxon_info.tsv

grep "Pseudom_aeru" YggL.replace.cross.tsv |
  cut -f 2 >YggL.Pseudom_aeru.protein.tsv

printf "She_balt_OS678_GCF_000178875_2_WP_006082529\nAerom_jan_GCF_016127195_1_WP_033112740\n" \
  >>YggL.Pseudom_aeru.protein.tsv

faops some ${PREFIX}/PROTEINS/all.replace.fa \
  <(cat YggL.Pseudom_aeru.protein.tsv) \
  YggL.Pseudom_aeru.protein.fa

mafft --auto YggL.Pseudom_aeru.protein.fa \
  >YggL.Pseudom_aeru.protein.mafft.fa

bsub -n 24 -J "iq9" ../iqtree-2.2.0-Linux/bin/iqtree2 \
  -s PROTEINS/YggL.Pseudom_aeru.protein.mafft.fa \
  --prefix test9 -T 5 -B 1000 -bnni -m Q.pfam+G4

bsub -n 24 -J "Pseudom_aeru" ../iqtree-2.2.0-Linux/bin/iqtree2 \
  -s PROTEINS/YggL.Pseudom_aeru.protein.mafft.fa \
  --prefix Pseudom_aeru -T 5 -b 100 -m Q.pfam+G4

grep -f <(cut -f 2 YggL.Pseudomonas.strains.tsv) \
  YggL.replace.cross.tsv |
  grep -v "Pseudom_aeru_PAO1_GCF_013001005_1" |
  cut -f 2 \
    >YggL.Pseudomonas.protein.tsv

printf "She_balt_OS678_GCF_000178875_2_WP_006082529\nAerom_jan_GCF_016127195_1_WP_033112740\n" \
  >>YggL.Pseudomonas.protein.tsv

faops some ${PREFIX}/PROTEINS/all.replace.fa \
  <(cat YggL.Pseudomonas.protein.tsv) \
  PROTEINS/YggL.Pseudomonas.protein.fa

mafft --auto PROTEINS/YggL.Pseudomonas.protein.fa \
  >PROTEINS/YggL.Pseudomonas.protein.mafft.fa

bsub -n 24 -J "iq10" ../iqtree-2.2.0-Linux/bin/iqtree2 \
  -s PROTEINS/YggL.Pseudomonas.protein.mafft.fa \
  --prefix test10 -T 5 -B 1000 -bnni -m Q.pfam+G4

bsub -n 24 -J "Pseudomonas" ../iqtree-2.2.0-Linux/bin/iqtree2 \
  -s PROTEINS/YggL.Pseudomonas.protein.mafft.fa \
  --prefix Pseudomonas -T 5 -b 100 -m Q.pfam+G4

for S in \
  Pseudom_aeru_PAO1 \
  Pseudom_aeru_PA7_GCF_000017205_1 \
  Pseudom_aeru_UCBPP_PA14_GCF_000014625_1 \
  Pseudom_aeru_PAK_GCF_902172305_2 \
  Pseudom_aeru_LESB58_GCF_000026645_1; do
  echo ${S}
done >typical.lst

while IFS= read -r name; do
  echo "${name}"
  cp ${PREFIX}/ASSEMBLY/"${name}"/*_genomic.gbff.gz ASSEMBLY/"${name}".gbff.gz
  mv ASSEMBLY/"${name}".gbff.gz genebank/.
done <typical.lst

gunzip genebank/*.gz

for i in LESB58 PA7 PAK PAO1 UCBPP_PA14; do
  python fetch_gbk.py -i genebank/*${i}*.gbff \
    -e 5000 -l YggL \
    -o ./${i}_YggL
done

grep -v "Pseudom_aeru_PAO1_GCF_013001005_1" \
  YggL.replace.cross.tsv |
  grep -f <(cut -f 2 YggL.Pseudomonas.protein.tsv) |
  cut -f 2 | tsv-join -d 1 -f ${PREFIX}/PROTEINS/all.strain.tsv -k 1 --append-fields 2 |
  cut -f 2 | sort | uniq \
  >YggL.Pseudomonas.bac120.species.tsv

faops some ${PREFIX}/PROTEINS/bac120.trim.fa \
  YggL.Pseudomonas.bac120.species.tsv \
  YggL.Pseudomonas.bac120.fa

mafft --auto YggL.Pseudomonas.bac120.fa \
  >YggL.Pseudomonas.bac120.aln.mafft.fa

bsub -q mpi -n 24 -J "iq" ./iqtree2 -s branched-chain.Pseudomonas.bac120.aln.mafft.fa -m MFP --prefix branched-chain.Pseudomonas.bac120 -T 20 -B 1000 --bnni -o She_balt_GCF_003030925_1,Vi_cho_GCF_008369605_1

mv branched-chain.Pseudomonas.bac120.treefile branched-chain.Pseudomonas.bac120.newick

python fetch_gbk.py -i genebank2/Pseudom_aeru_PAO1*.gbff \
  -e 5000 -l PA3046 \
  -o ./${i}_YggL

while IFS= read -r name; do
  echo "${name}"
  cp ${PREFIX}/ASSEMBLY/"${name}"/*_genomic.gbff.gz ASSEMBLY/"${name}".gbff.gz
  mv ASSEMBLY/"${name}".gbff.gz genebank3/.
done <gamma.lst

gunzip genebank3/*.gz

## Gamma tree
grep -f <(cut -f 2 Gammaproteobacteria.strain.tsv) \
  YggL.replace.cross.tsv |
  cut -f 2 | tsv-join -d 1 -f ${PREFIX}/PROTEINS/all.strain.tsv -k 1 --append-fields 2 |
  cut -f 2 | sort | uniq \
  >YggL.Gammaproteobacteria.bac120.tsv

printf "Bac_subti_subtilis_168\nSta_aure_aureus_NCTC_8325\n" \
  >>YggL.Gammaproteobacteria.bac120.tsv
# 317

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  parallel --no-run-if-empty --linebuffer -k -j 8 '
    grep PROTEINS/{}/{}.replace.tsv ' \
    '  -f YggL.Gammaproteobacteria.bac120.tsv ' \
    '  >PROTEINS/{}/{}.Gammaproteobacteria.tsv
    faops some /share/home/wangq/data/Pseudomonas/PROTEINS/all.uniq.fa ' \
    '  <(cut -f 1 PROTEINS/{}/{}.Gammaproteobacteria.tsv |
        tsv-uniq) stdout ' \
    '  >PROTEINS/{}/{}.Gammaproteobacteria.fa '

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  parallel --no-run-if-empty --linebuffer -k -j 8 '
    mafft --auto PROTEINS/{}/{}.Gammaproteobacteria.fa >PROTEINS/{}/{}.Gammaproteobacteria.aln.fa
    '

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  while IFS= read -r marker; do
    parallel --no-run-if-empty --linebuffer -k -j 8 "
      faops replace -s PROTEINS/${marker}/${marker}.Gammaproteobacteria.aln.fa \
        <(echo {}) stdout
        " <PROTEINS/"${marker}"/"${marker}".Gammaproteobacteria.tsv \
      >PROTEINS/"${marker}"/"${marker}".Gammaproteobacteria.replace.fa
  done

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  while IFS= read -r marker; do
    faops filter -l 0 PROTEINS/"${marker}"/"${marker}".Gammaproteobacteria.replace.fa stdout
    echo
  done >PROTEINS/bac120.Gammaproteobacteria.aln.fas

fasops concat PROTEINS/bac120.Gammaproteobacteria.aln.fas \
  YggL.Gammaproteobacteria.bac120.tsv \
  -o PROTEINS/bac120.Gammaproteobacteria.aln.fa

trimal -automated1 \
  -in PROTEINS/bac120.Gammaproteobacteria.aln.fa \
  -out PROTEINS/bac120.Gammaproteobacteria.trim.fa

bsub -q fat_384 -n 80 -J "GammaTree" ../iqtree-2.2.0-Linux/bin/iqtree2 \
  -s PROTEINS/bac120.Gammaproteobacteria.trim.fa \
  --prefix Gammaproteobacteria -T 48 -B 1000 -bnni

## Gamma tree
grep -f <(cut -f 2 Gammaproteobacteria.strain.tsv) \
  YggL.replace.cross.tsv |
  cut -f 2 | tsv-join -d 1 -f ${PREFIX}/PROTEINS/all.strain.tsv -k 1 --append-fields 2 |
  cut -f 2 | sort | uniq \
  >YggL.Gammaproteobacteria.bac120.tsv

printf "Bac_subti_subtilis_168\nSta_aure_aureus_NCTC_8325\n" \
  >>YggL.Gammaproteobacteria.bac120.tsv
# 317

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  parallel --no-run-if-empty --linebuffer -k -j 8 '
    grep PROTEINS/{}/{}.replace.tsv ' \
    '  -f YggL.Gammaproteobacteria.bac120.tsv ' \
    '  >PROTEINS/{}/{}.Gammaproteobacteria.tsv
    faops some /share/home/wangq/data/Pseudomonas/PROTEINS/all.uniq.fa ' \
    '  <(cut -f 1 PROTEINS/{}/{}.Gammaproteobacteria.tsv |
        tsv-uniq) stdout ' \
    '  >PROTEINS/{}/{}.Gammaproteobacteria.fa '

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  parallel --no-run-if-empty --linebuffer -k -j 8 '
    mafft --auto PROTEINS/{}/{}.Gammaproteobacteria.fa >PROTEINS/{}/{}.Gammaproteobacteria.aln.fa
    '

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  while IFS= read -r marker; do
    parallel --no-run-if-empty --linebuffer -k -j 8 "
      faops replace -s PROTEINS/${marker}/${marker}.Gammaproteobacteria.aln.fa \
        <(echo {}) stdout
        " <PROTEINS/"${marker}"/"${marker}".Gammaproteobacteria.tsv \
      >PROTEINS/"${marker}"/"${marker}".Gammaproteobacteria.replace.fa
  done

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  while IFS= read -r marker; do
    faops filter -l 0 PROTEINS/"${marker}"/"${marker}".Gammaproteobacteria.replace.fa stdout
    echo
  done >PROTEINS/bac120.Gammaproteobacteria.aln.fas

fasops concat PROTEINS/bac120.Gammaproteobacteria.aln.fas \
  YggL.Gammaproteobacteria.bac120.tsv \
  -o PROTEINS/bac120.Gammaproteobacteria.aln.fa

trimal -automated1 \
  -in PROTEINS/bac120.Gammaproteobacteria.aln.fa \
  -out PROTEINS/bac120.Gammaproteobacteria.trim.fa

bsub -q fat_384 -n 80 -J "GammaTree" ../iqtree-2.2.0-Linux/bin/iqtree2 \
  -s PROTEINS/bac120.Gammaproteobacteria.trim.fa \
  --prefix Gammaproteobacteria -T 48 -B 1000 -bnni

## Pseudomonas tree
grep -v "Pseudom_aeru_PAO1_GCF_013001005_1" \
  YggL.replace.cross.tsv |
  grep -f <(cut -f 2 YggL.Pseudomonas.protein.tsv) |
  cut -f 2 | tsv-join -d 1 -f ${PREFIX}/PROTEINS/all.strain.tsv -k 1 --append-fields 2 |
  cut -f 2 | sort | uniq \
  >YggL.Pseudomonas.bac120.species.tsv

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  parallel --no-run-if-empty --linebuffer -k -j 8 '
    grep PROTEINS/{}/{}.replace.tsv ' \
    '  -f YggL.Pseudomonas.bac120.species.tsv ' \
    '  >PROTEINS/{}/{}.Pseudomonas.tsv
    faops some /share/home/wangq/data/Pseudomonas/PROTEINS/all.uniq.fa ' \
    '  <(cut -f 1 PROTEINS/{}/{}.Pseudomonas.tsv |
        tsv-uniq) stdout ' \
    '  >PROTEINS/{}/{}.Pseudomonas.fa '

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  parallel --no-run-if-empty --linebuffer -k -j 8 '
    mafft --auto PROTEINS/{}/{}.Pseudomonas.fa >PROTEINS/{}/{}.Pseudomonas.aln.fa '

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  while IFS= read -r marker; do
    parallel --no-run-if-empty --linebuffer -k -j 8 "
      faops replace -s PROTEINS/${marker}/${marker}.Pseudomonas.aln.fa \
        <(echo {}) stdout
        " <PROTEINS/"${marker}"/"${marker}".Pseudomonas.tsv \
      >PROTEINS/"${marker}"/"${marker}".Pseudomonas.replace.fa
  done

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  while IFS= read -r marker; do
    faops filter -l 0 PROTEINS/"${marker}"/"${marker}".Pseudomonas.replace.fa stdout
    echo
  done >PROTEINS/bac120.Pseudomonas.aln.fas

fasops concat PROTEINS/bac120.Pseudomonas.aln.fas \
  YggL.Pseudomonas.bac120.species.tsv \
  -o PROTEINS/bac120.Pseudomonas.aln.fa

trimal -automated1 \
  -in PROTEINS/bac120.Pseudomonas.aln.fa \
  -out PROTEINS/bac120.Pseudomonas.trim.fa

bsub -q largemem -n 24 -J "PseudomonasTree" ../iqtree-2.2.0-Linux/bin/iqtree2 \
  -s PROTEINS/bac120.Pseudomonas.trim.fa \
  --prefix Pseudomonas -T 22 -B 1000 -bnni

## Pseudom_aeru tree
grep -v "Pseudom_aeru_PAO1_GCF_013001005_1" \
  YggL.replace.cross.tsv |
  grep -f <(cut -f 2 YggL.Pseudom_aeru.protein.tsv) |
  cut -f 2 | tsv-join -d 1 -f ${PREFIX}/PROTEINS/all.strain.tsv -k 1 --append-fields 2 |
  cut -f 2 | sort | uniq \
  >YggL.Pseudom_aeru.bac120.species.tsv

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  parallel --no-run-if-empty --linebuffer -k -j 8 '
    grep PROTEINS/{}/{}.replace.tsv ' \
    '  -f YggL.Pseudom_aeru.bac120.species.tsv ' \
    '  >PROTEINS/{}/{}.Pseudom_aeru.tsv
    faops some /share/home/wangq/data/Pseudomonas/PROTEINS/all.uniq.fa ' \
    '  <(cut -f 1 PROTEINS/{}/{}.Pseudom_aeru.tsv |
        tsv-uniq) stdout ' \
    '  >PROTEINS/{}/{}.Pseudom_aeru.fa '

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  parallel --no-run-if-empty --linebuffer -k -j 8 '
    mafft --auto PROTEINS/{}/{}.Pseudom_aeru.fa >PROTEINS/{}/{}.Pseudom_aeru.aln.fa '

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  while IFS= read -r marker; do
    parallel --no-run-if-empty --linebuffer -k -j 8 "
      faops replace -s PROTEINS/${marker}/${marker}.Pseudom_aeru.aln.fa \
        <(echo {}) stdout
        " <PROTEINS/"${marker}"/"${marker}".Pseudom_aeru.tsv \
      >PROTEINS/"${marker}"/"${marker}".Pseudom_aeru.replace.fa
  done

sed '1d' </share/home/wangq/Scripts/withncbi/hmm/bac120.tsv |
  cut -f 1 | cut -d . -f 1 |
  while IFS= read -r marker; do
    faops filter -l 0 PROTEINS/"${marker}"/"${marker}".Pseudom_aeru.replace.fa stdout
    echo
  done >PROTEINS/bac120.Pseudom_aeru.aln.fas

fasops concat PROTEINS/bac120.Pseudom_aeru.aln.fas \
  YggL.Pseudom_aeru.bac120.species.tsv \
  -o PROTEINS/bac120.Pseudom_aeru.aln.fa

trimal -automated1 \
  -in PROTEINS/bac120.Pseudom_aeru.aln.fa \
  -out PROTEINS/bac120.Pseudom_aeru.trim.fa

bsub -q fat_384 -n 80 -J "Pseudom_aeruTree" ../iqtree-2.2.0-Linux/bin/iqtree2 \
  -s PROTEINS/bac120.Pseudom_aeru.trim.fa \
  --prefix Pseudom_aeru -T AUTO -B 1000 -bnni

YggL.Pseudom_aeru.protein.tsv
