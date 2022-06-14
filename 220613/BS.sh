#!/usr/bin/env bash

DATA_FILE=$1
SAMPLE=$2
BS_DIR=${3:-BS}

if [ ! -e "${DATA_FILE}" ]; then
  echo "[${DATA_FILE}] is not a file" 1>&2
  exit
fi

mkdir -p "${BS_DIR}"

echo "==> ${BS_DIR}/data.tsv" 1>&2
parallel --no-run-if-empty --linebuffer -k -j 2 "
  echo '==> Bootstrap #{}' 1>&2
  tsv-sample -H -r -n ${SAMPLE} \
    <${DATA_FILE} |
    keep-header -- sort -k1,1 --buffer-size=2G \
      >${BS_DIR}/data.tsv.{}
    " ::: $(printf "%03d " {0..99})

