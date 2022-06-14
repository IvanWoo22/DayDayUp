#!/usr/bin/env bash

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

  tsv-summarize -H --min 6 --median 6 "${target}" |
    datamash transpose

  tsv-summarize -H --max 7 --median 7 "${target}" |
    datamash transpose
} >>stat.tmp

mlr --itsv --omd cat stat.tmp
rm stat.tmp
