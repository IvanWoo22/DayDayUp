#!/usr/bin/env bash
BS_DIR=$1
RESULT_FILE=$2
TARGET=$3
MIN=$4
MAX=$5
OUTPUT_DIR=${6:-.}

if [ ! -d "${BS_DIR}" ]; then
  echo "[${BS_DIR}] is not a dir" 1>&2
  exit
fi

if [ ! -e "${RESULT_FILE}" ]; then
  echo "[${RESULT_FILE}] is not a file" 1>&2
  exit
fi

BASENAME=$(basename "${RESULT_FILE}")
parallel --no-run-if-empty --linebuffer -k -j 22 "
  echo '==> Bootstrap #{}' 1>&2
  Rscript multivalidate.R ${TARGET} ${BS_DIR}/data.tsv.{} ${RESULT_FILE} ${RESULT_FILE}.{}
  keep-header -- awk '\$6>${MIN}||\$6<${MAX}' \
    <${RESULT_FILE}.{} |
    cut -f 1 \
      >${OUTPUT_DIR}/${BASENAME}.{}
    " ::: $(printf "%03d " {0..99})
rm "${RESULT_FILE}".*[0-9]

echo '==> Outputs' 1>&2
parallel --no-run-if-empty --linebuffer -k -j 22 "
    cat ${OUTPUT_DIR}/${BASENAME}.{}
    rm  ${OUTPUT_DIR}/${BASENAME}.{}
    " ::: $(printf "%03d " {0..99}) |
  grep -v "^ProbeID" |
  sort --buffer-size=2G |
  tsv-summarize --group-by 1 --count |
  sort -k2,2nr --buffer-size=2G |
  (echo -e "ProbeID\tBS" && cat) \
    >"${OUTPUT_DIR}"/"${BASENAME}".count.tsv

tsv-join \
  --filter-file "${OUTPUT_DIR}"/"${BASENAME}".count.tsv \
  --key-fields 1 --append-fields 2 \
  <"${RESULT_FILE}" \
  >"${OUTPUT_DIR}"/"${BASENAME}".bootstrap.tsv
