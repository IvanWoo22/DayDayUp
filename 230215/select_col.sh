#!/usr/bin/env bash

usage() {
  echo "bash select_col.sh [-f 1-3] <infile> [select]" 1>&2
  echo "bash select_col.sh -f 1,5 tests/SKCM/simple.data.tsv | datamash check" 1>&2
  echo "# 6 lines, 2 fields" 1>&2
  echo "bash select_col.sh -f 1-3 tests/SKCM/meth.tsv.gz tests/SKCM/ucox.05.result.tsv | datamash check" 1>&2
  echo "# 303 lines, 78 fields" 1>&2
  exit 1
}
[ $# -eq 0 ] && usage

while getopts ":f:" opt; do
  case ${opt} in
  f)
    opt_field=$OPTARG
    ;;
  \?)
    echo "Invalid option: $OPTARG" 1>&2
    usage
    ;;
  :)
    echo "Invalid option: $OPTARG requires an argument" 1>&2
    usage
    ;;
  esac
done
shift $((OPTIND - 1))

if [ ! -f "$1" ]; then
  echo "[$1] is not a file" 1>&2
  exit
fi
opt_infile="$1"

if [ "$2" != "" ]; then
  if [ ! -f "$2" ]; then
    echo "[$2] is not a file" 1>&2
    exit
  fi
  opt_select="$2"
fi

#----------------------------#
# run
#----------------------------#
# tmpdir
mytmpdir=$(mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir')

touch "${mytmpdir}"/list

if [ "${opt_field}" != "" ]; then
  echo "${opt_field}" >>"${mytmpdir}"/list
fi

if [ "${opt_select}" != "" ]; then
  gzip -dcf "${opt_infile}" |
    head -n 1 |
    tr $'\t' '\n' |
    nl |
    perl -nl -e '
            my @fields = grep {/\w+/} split /\s+/, $_;
            next unless @fields == 2;
            print join qq{\t}, $fields[1], $fields[0];
        ' \
      >"${mytmpdir}"/fields.index

  tsv-select -f 1 \
    <"${opt_select}" \
    >"${mytmpdir}"/fields.select

  grep -Fw -f "${mytmpdir}"/fields.select \
    "${mytmpdir}"/fields.index |
    cut -f 2 \
      >>"${mytmpdir}"/list
fi

perl -nl -MAlignDB::IntSpan -e '
        BEGIN {
            our $set = AlignDB::IntSpan->new();
        }
        $set->add($_);
        END {
            my @elements = $set->elements;
            while ( scalar @elements ) {
                my @batching = splice @elements, 0, 5000;
                my $batching_set = AlignDB::IntSpan->new;
                $batching_set->add(@batching);
                print $batching_set->runlist();
            }
        }
    ' \
  <"${mytmpdir}"/list >"${mytmpdir}"/runlist

count=$(wc -l <"${mytmpdir}"/runlist)

if [ "$count" -eq "0" ]; then
  echo >&2 "No fields"
elif [ "$count" -eq "1" ]; then
  echo >&2 "Writing fields..."
  gzip -dcf "${opt_infile}" |
    tsv-select -f $(cat "${mytmpdir}"/runlist)
else
  parallel --no-run-if-empty --line-buffer -k -j 1 --seqreplace ,, "
            echo >&2 'Writing fields ,,...'
            gzip -dcf ${opt_infile} |
                tsv-select -f {} |
                datamash transpose
        " \
    <"${mytmpdir}"/runlist \
    >"${mytmpdir}"/selected

  echo >&2 "Writing all fields..."
  datamash transpose <"${mytmpdir}"/selected
fi

# clean
rm -fr "${mytmpdir}"

