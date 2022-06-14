PREFIX=/share/home/wangq/data/Pseudomonas
ID=IPR004361
NAME="GlyI"

grep -h "${ID}" ${PREFIX}/STRAINS/Pseudom_aeru_PAO1/*.fa.tsv | cut -f 1,4,5

## Local run and upload to HPCC.
aria2c -x 12 https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR00068.2.HMM -o hmm/${NAME}.tigrfam.hmm

E_VALUE="1e-20"
for DB in tigrfam; do
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

faops some ${PREFIX}/PROTEINS/all.replace.fa \
  <(cut -f 2 ${NAME}.replace.tsv) \
  ${NAME}.replace.fa

for DB in TIGR; do
  bsub -n 24 -J "${DB}" "
  E_VALUE=\"1e-10\"
  NAME=\"GlyI\"
  hmmscan --cpu 20 -E \${E_VALUE} --domE \${E_VALUE} --noali \
    -o \${NAME}.${DB}.progress.txt --tblout \${NAME}.${DB}.tbl \
    hmm/${DB}/HMM \${NAME}.replace.fa
    "
done

for target in ad mci hc; do
  tsv-append -H ${target}_training/*.tsv >1_training/${target}.result.tsv
done
