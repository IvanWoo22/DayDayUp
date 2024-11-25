mv ICU医嘱信息.csv ICU_DO.csv
sed ':a;N;$!ba;s/\r\n)/)/g' ICU_DO.csv \
	| sed 's/\r//g' >ICU_DO1.csv
bash fix_csv.sh ICU_DO1.csv ICU_DO2.csv
vim ICU_DO2.csv
#Change a line manually.

awk -F'\t' '{sub($2, "", $1); print $2"\t"$1"\t"$3"\t"$4"\t"$5}' ICU_Cases_raw.tsv >ICU_Cases_tmp1.tsv

sed -i 's/℃/°C/g' ICU_Cases_tmp1.tsv
sed -i 's/：/:/g' ICU_Cases_tmp1.tsv
sed -i 's/（/(/g' ICU_Cases_tmp1.tsv
sed -i 's/）/)/g' ICU_Cases_tmp1.tsv
vim ICU_Cases_tmp1.tsv
#Change many lines manually.

perl basic_info_extract.pl <ICU_Cases_tmp1.tsv >ICU_Cases_tmp2.tsv

sed -i 's/－/-/g' ICU_Cases_tmp2.tsv
sed -i 's/―/-/g' ICU_Cases_tmp2.tsv
sed -i 's/--/-/g' ICU_Cases_tmp2.tsv
sed -i 's/、//g' ICU_Cases_tmp2.tsv
sed -i 's/|//g' ICU_Cases_tmp2.tsv
sed -i 's///g' ICU_Cases_tmp2.tsv
sed -i 's/不祥/不详/g' ICU_Cases_tmp2.tsv
sed -i 's/不详/未提供/g' ICU_Cases_tmp2.tsv
sed -i 's/为提供/未提供/g' ICU_Cases_tmp2.tsv
sed -i 's/未知/未提供/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)未通过\t/\1未提供\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)未提交\t/\1未提供\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)-\t/\1未提供\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)\t/\1未提供\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)—\t/\1未提供\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)不详\t/\1未提供\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)其他\t/\1未提供\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)无\t/\1无业人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)无业\t/\1无业人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)无职业\t/\1无业人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)个体经营\t/\1个体经营者\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)个体\t/\1个体经营者\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)职工\t/\1职员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)员工\t/\1职员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)离休\t/\1退(离)休人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)退休\t/\1退(离)休人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)离退休人员\t/\1退(离)休人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)离退休\t/\1退(离)休人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)退（离）休人员\t/\1退(离)休人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)已退休\t/\1退(离)休人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)退休职工\t/\1退(离)休人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)退休职员\t/\1退(离)休人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)退休人员\t/\1退(离)休人员\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)医师\t/\1医生\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)务农\t/\1农民\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\([^\t]*\t[^\t]*\t\)自由职业\t/\1自由职业者\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\(\([^\t]*\t\)\{4\}\)\t/\1未提供\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\(\([^\t]*\t\)\{5\}[[:digit:]]\+\)岁/\1/g' ICU_Cases_tmp2.tsv
sed -i 's/\(\([^\t]*\t\)\{7\}\)已\t/\1已婚\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\(\([^\t]*\t\)\{7\}\)\t/\1未提供\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\(\([^\t]*\t\)\{7\}\)离异\t/\1离婚\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\(\([^\t]*\t\)\{8\}\)汉\t/\1汉族\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\(\([^\t]*\t\)\{8\}\)\t/\1未提供\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\(\([^\t]*\t\)\{9\}\)\t/\1未提供\t/g' ICU_Cases_tmp2.tsv
sed -i 's/\(\([^\t]*\t\)\{10\}\)\t/\1未提供\t/g' ICU_Cases_tmp2.tsv

perl -pe 's/(\d{4})年(\d{2})月(\d{2})日(\d{2})时(\d{2})分/$1-$2-$3 $4:$5:00/g' ICU_Cases_tmp2.tsv \
	| perl fullconvert.pl | awk '{ gsub(/\xef\xbb\xbf/,""); print }' \
	| sort -t$'\t' -k1,1 -k15,15n -k14,14 | uniq >ICU_Cases_tmp3.tsv
awk -F "\t" '{print $1 "\t" $2 "\t" $4 "\t" $6 "\t" $8 "\t" $11 "\t" $3 "\t" $10 "\t" $7 "\t" $5}' ICU_Cases_tmp3.tsv \
	>ICU_Cases_tmp4.tsv
cat case_info_clean.tsv >>ICU_Cases_tmp4.tsv
perl basic_info_filter.pl <ICU_Cases_tmp4.tsv >ICU_Cases_clean.tsv

awk -F '\t' '{sub($2, "", $1); print $2"\t"$1}' ICU_Checks_tmp.tsv \
	>ICU_Checks_tmp1.tsv
