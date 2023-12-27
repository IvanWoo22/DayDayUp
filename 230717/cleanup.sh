perl basic_info_extract.pl <cases.tsv >basic_info.tmp

sed -i 's/－/-/g' basic_info.tmp
sed -i 's/―/-/g' basic_info.tmp
sed -i 's/--/-/g' basic_info.tmp
sed -i 's/：/:/g' basic_info.tmp
sed -i 's/不祥/不详/g' basic_info.tmp
sed -i 's/\([^\t]*\t[^\t]*\t\)个体经营\t/\1个体经营者\t/g' basic_info.tmp
sed -i 's/\([^\t]*\t[^\t]*\t\)离休\t/\1退(离)休人员\t/g' basic_info.tmp
sed -i 's/\([^\t]*\t[^\t]*\t\)退休\t/\1退(离)休人员\t/g' basic_info.tmp
sed -i 's/\([^\t]*\t[^\t]*\t\)离退休人员\t/\1退(离)休人员\t/g' basic_info.tmp
sed -i 's/\([^\t]*\t[^\t]*\t\)离退休\t/\1退(离)休人员\t/g' basic_info.tmp
sed -i 's/\([^\t]*\t[^\t]*\t\)退（离）休人员\t/\1退(离)休人员\t/g' basic_info.tmp
sed -i 's/\([^\t]*\t[^\t]*\t\)自由职业\t/\1自由职业者\t/g' basic_info.tmp
sed -i 's/\(\([^\t]*\t\)\{5\}[[:digit:]]\+\)岁/\1/g' basic_info.tmp
perl -pe 's/(\d{4})年(\d{2})月(\d{2})日(\d{2})时(\d{2})分/$1-$2-$3 $4:$5:00/g' basic_info.tmp \
	| perl fullconvert.pl | awk '{ gsub(/\xef\xbb\xbf/,""); print }' | sort -t$'\t' -k1,1 -k15,15n -k14,14 | uniq >case_info.tsv
# vim case_info.tsv

###
# en.raw.tsv 管饲营养信息
awk -F "\t" 'NR!=1{print $2}' en.raw.tsv | sort -n | uniq >id1.list
# cases.tsv 病例信息
cut -f 1 cases.tsv | sort -n | uniq >id2.list
# tests.tsv 检验信息
awk -F "\t" 'NR!=1{print $1}' tests.tsv | sort -n | uniq >id3.list
