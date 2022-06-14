#!/usr/bin/env bash

usage() {
  echo "bash result_stat.sh [-b] <dir>" 1>&2
  exit
}
[ $# -eq 0 ] && usage

while getopts ":b" opt; do
  case ${opt} in
  b)
    bootstrap=1
    ;;
  *)
    usage
    ;;
  esac
done

shift $((OPTIND - 1))

if [ ! -e "$*" ]; then
  echo "[$*] is not a file" 1>&2
  exit
fi
target="$*"

echo -e "#Item\tValue" >stat.tmp

{
  wc -l "${target}" |
    datamash reverse --whitespace |
    grep -v "^total"

  cut -f 1 "${target}" |
    perl -nl -e 'print for split /\+/' |
    tsv-uniq |
    tsv-summarize -H --count |
    datamash transpose

  tsv-summarize -H --max 5 --min 5 "${target}" |
    datamash transpose

  tsv-summarize -H --min 6 --max 6 "${target}" |
    datamash transpose

  tsv-summarize -H --min 7 --max 7 "${target}" |
    datamash transpose
} >>stat.tmp

if [ ${bootstrap} ]; then
  tsv-select -H --fields 8 <"${target}" |
    keep-header -- sort -k1,1nr --buffer-size=2G |
    datamash -H --full bin:5 1 |
    tsv-summarize -H --group-by 2 --count \
      >>stat.tmp
fi

mlr --itsv --omd cat stat.tmp
rm stat.tmp
