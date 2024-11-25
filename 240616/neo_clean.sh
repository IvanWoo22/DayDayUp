TEST_FILE='重症医学科患者数据-检验.csv'
sed -i 's/－/-/g' ${TEST_FILE}
sed -i 's/―/-/g' ${TEST_FILE}
sed -i 's/--/-/g' ${TEST_FILE}
sed -i 's/、/;/g' ${TEST_FILE}
sed -i 's/，/;/g' ${TEST_FILE}
sed -i 's/。/./g' ${TEST_FILE}
sed -i 's/｜/|/g' ${TEST_FILE}
sed -i 's/℃/°C/g' ${TEST_FILE}
sed -i 's/：/:/g' ${TEST_FILE}
sed -i 's/（/(/g' ${TEST_FILE}
sed -i 's/）/)/g' ${TEST_FILE}

sed -i 's/A,B1,B2,B6,叶酸,B12,C,D,E/A;B1;B2;B6;叶酸;B12;C;D;E/g' ${TEST_FILE}
sed -i 's/B1,B2,B6,叶酸,B12,C/B1;B2;B6;叶酸;B12;C/g' ${TEST_FILE}
sed -i 's/B1,叶酸,B12,D/B1;叶酸;B12;D/g' ${TEST_FILE}
sed -i 's/ADP,血栓弹力图法/ADP;血栓弹力图法/g' ${TEST_FILE}
sed -i 's/纤溶亢进预测值 /纤溶亢进预测值/g' ${TEST_FILE}
sed -i 's/抗IgG,C3d/抗IgG;C3d/g' ${TEST_FILE}
sed -i 's/Ⅷ因子,Ⅸ因子/Ⅷ因子;Ⅸ因子/g' ${TEST_FILE}
sed -i 's/浑浊,有凝块/浑浊;有凝块/g' ${TEST_FILE}
sed -i 's/>=0.10,阳性/>=0.10:阳性/g' ${TEST_FILE}
sed -i 's/<0.10,阴性/<0.10:阴性/g' ${TEST_FILE}
sed -i 's/0.73,无反应性/0.73;无反应性/g' ${TEST_FILE}
sed -i 's/5.09,无反应性/5.09;无反应性/g' ${TEST_FILE}
sed -i 's/IgG,κ型单克隆免疫球蛋白/IgG κ型单克隆免疫球蛋白/g' ${TEST_FILE}
sed -i 's/IgA,λ型单克隆免疫球蛋白/IgA λ型单克隆免疫球蛋白/g' ${TEST_FILE}
sed -i 's/IgG,λ型单克隆免疫球蛋白/IgG λ型单克隆免疫球蛋白/g' ${TEST_FILE}
sed -i 's/较稀,带粘液,NULL/较稀;带粘液,NULL/g' ${TEST_FILE}
sed -i 's/软便,沾有粘液,NULL/软便;沾有粘液,NULL/g' ${TEST_FILE}
sed -i 's/软便,有水,NULL/软便;有水,NULL/g' ${TEST_FILE}
sed -i 's/硬便,沾有粘液,NULL/硬便;沾有粘液,NULL/g' ${TEST_FILE}
sed -i 's/较清亮,有凝块,NULL/较清亮;有凝块,NULL/g' ${TEST_FILE}
sed -i 's/缺乏;<20,不足:20-30,正常/缺乏:<20;不足:20-30;正常/g' ${TEST_FILE}
sed -i 's/0.14,无反应性/0.14;无反应性/g' ${TEST_FILE}
sed -i 's/4.98,无反应性/4.98;无反应性/g' ${TEST_FILE}
sed -i 's/脓细胞2-6\/Hp,红细胞(+)/脓细胞2-6\/Hp;红细胞(+)/g' ${TEST_FILE}
sed -i 's/细胞少, 不宜分析/细胞少;不宜分析/g' ${TEST_FILE}
sed -i 's/细胞少, 不宜分类/细胞少;不宜分类/g' ${TEST_FILE}
sed -i 's/絮状凝结,无法计数/絮状凝结;无法计数/g' ${TEST_FILE}
sed -i 's/肾移植配型PRA-Ⅰ,PRA-Ⅱ/肾移植配型PRA-Ⅰ;PRA-Ⅱ/g' ${TEST_FILE}
sed -i 's/健康女性:停经前,>20岁:11-43|         停经后:15-46|骨质疏松患者:13-48/健康女性:停经前:>20岁:11-43|停经后:15-46;骨质疏松患者:13-48/g' ${TEST_FILE}
sed -i 's/健康女性:停经前,>20岁:11-43; 停经后:15-46|骨质疏松患者:13-48/健康女性:停经前:>20岁:11-43|停经后:15-46;骨质疏松患者:13-48/g' ${TEST_FILE}
sed -i 's/血清中检出抗-E抗体,不排除/血清中检出抗-E抗体;不排除/g' ${TEST_FILE}

tr -d '\r' <${TEST_FILE} \
	| awk -F, 'NR>1&&length($1)>=7{print substr($1,1,7), substr($1,8),$2,$3,$4,$5,$6,$7,$8}' \
		OFS=, >ICU_test.csv

tr -d '\r' <重症医学科患者数据-文书.csv \
	| awk -F',' 'NR>1&&length($2)==7{
			# 提取第 1, 2, 3 列
			col1=$1
			col3=$3
			gsub(/ /, "", col3)
			detail_info=""
			if (match($0, /(<DetailInfo>.+<\/DetailInfo>)/, arr)) {
					detail_info=arr[1]
					gsub(/,/, "，", detail_info)
					gsub(/\xE3\x80\x80/, " ", detail_info)
					gsub(/ +/, " ", detail_info)
			}
			print substr(col1,1,7), substr(col1,8), col3, detail_info
	}' OFS=',' >ICU_text.csv
TEST_FILE='ICU_text.csv'
sed -i 's/－/-/g' ${TEST_FILE}
sed -i 's/―/-/g' ${TEST_FILE}
sed -i 's/--/-/g' ${TEST_FILE}
sed -i 's/。/./g' ${TEST_FILE}
sed -i 's/｜/|/g' ${TEST_FILE}
sed -i 's/℃/°C/g' ${TEST_FILE}
sed -i 's/：/:/g' ${TEST_FILE}
sed -i 's/（/(/g' ${TEST_FILE}
sed -i 's/）/)/g' ${TEST_FILE}
sed -i 's///g' ${TEST_FILE}
sed -i 's/不祥/不详/g' ${TEST_FILE}

TEST_FILE='重症医学科患者数据-检查.csv'
sed -i 's/－/-/g' ${TEST_FILE}
sed -i 's/―/-/g' ${TEST_FILE}
sed -i 's/--/-/g' ${TEST_FILE}
sed -i 's/。/. /g' ${TEST_FILE}
sed -i 's/｜/|/g' ${TEST_FILE}
sed -i 's/℃/°C/g' ${TEST_FILE}
sed -i 's/[[:space:]]*：/: /g' ${TEST_FILE}
sed -i 's/（/(/g' ${TEST_FILE}
sed -i 's/）/)/g' ${TEST_FILE}
sed -i 's///g' ${TEST_FILE}
sed -i 's/不祥/不详/g' ${TEST_FILE}
sed -i 's/,常规]/(常规)]/g' ${TEST_FILE}
sed -i 's/[[:space:]]\+/ /g' ${TEST_FILE}
sed -i 's/\(\[[^]]*\),\([^]]*\]\)/\1，\2/g' ${TEST_FILE}
sed -i 's/\r//g' ${TEST_FILE}
sed -i 's/胆囊,含穿刺活检/胆囊，含穿刺活检/g' ${TEST_FILE}
sed -i 's/血管评估,含门静脉/血管评估，含门静脉/g' ${TEST_FILE}
sed -i 's/NULL//g' ${TEST_FILE}

TEST_FILE='重症医学科患者数据-医嘱.csv'
sed -i 's/\r//g' ${TEST_FILE}
sed -E -i '/\(st$/ {N; s/\n//;}' ${TEST_FILE}
sed -E -i '/4mg\s*$/ {N; s/\n//;}' ${TEST_FILE}
sed -E -i '/00\s*$/ {N; s/\n//;}' ${TEST_FILE}
sed -E -i '/明抽\s*$/ {N; s/\n//;}' ${TEST_FILE}
sed -E -i ':a; s/\(([^\(\)]*),([^\(\)]*)\)/\(\1，\2\)/g; ta' ${TEST_FILE}
sed -i 's/（/(/g' ${TEST_FILE}
sed -i 's/）/)/g' ${TEST_FILE}
sed -i 's/标本检查与诊断,手术标本组织/标本检查与诊断，手术标本组织/g' ${TEST_FILE}
sed -i 's/检查(组套),手术标本/检查(组套)，手术标本/g' ${TEST_FILE}
sed -i 's/AA+ADP,血栓弹力图法/AA+ADP，血栓弹力图法/g' ${TEST_FILE}
sed -i 's/医生,/医生，/g' ${TEST_FILE}
sed -i 's/片,晚/片，晚/g' ${TEST_FILE}
sed -i 's/一次,分次推注,间隔/一次，分次推注，间隔/g' ${TEST_FILE}
sed -i 's/取石术,经输尿/取石术，经输尿/g' ${TEST_FILE}
sed -i 's/Ⅷ因子,Ⅸ因子/Ⅷ因子，Ⅸ因子/g' ${TEST_FILE}
sed -i 's/U,血浆/U，血浆/g' ${TEST_FILE}
sed -i 's/,停血滤每日减一次,/，停血滤每日减一次，/g' ${TEST_FILE}
sed -i 's/mg,qd/mg，qd/g' ${TEST_FILE}
sed -i 's/[[:space:]]\+/ /g' ${TEST_FILE}

TEST_FILE='重症医学科患者数据-基本信息.tsv'
sed -i 's/\r//g' ${TEST_FILE}
sed -i 's/（/(/g' ${TEST_FILE}
sed -i 's/）/)/g' ${TEST_FILE}
sed -i 's/NULL//g' ${TEST_FILE}
sed -i 's/,/，/g' ${TEST_FILE}
sed -i 's/\t/,/g' ${TEST_FILE}
sed -i 's/|,/,/g' ${TEST_FILE}
awk -F, 'NR>1&&length($1)>=7{print $2,substr($1,8),$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19}' \
	OFS=, ${TEST_FILE} >ICU_basic.csv

TEST_FILE='重症医学科患者数据-护理体征.csv'
sed -i 's/\r//g' ${TEST_FILE}
sed -i 's/（/(/g' ${TEST_FILE}
sed -i 's/）/)/g' ${TEST_FILE}
sed -i 's/NULL//g' ${TEST_FILE}
sed -i 's/℃/°C/g' ${TEST_FILE}
awk -F, 'NR>1&&length($1)>=7{print substr($1,1,7),substr($1,8),$2,$3,$4}' \
	OFS=, ${TEST_FILE} >ICU_nursingsigns.csv

awk -F, '{print $1,$2,$3,$4}' OFS='\t' ICU_text.csv | perl basic_info_extract.pl >ICU_text.tsv

TEST_FILE='ICU_text.tsv'
sed -i 's/－/-/g' ${TEST_FILE}
sed -i 's/―/-/g' ${TEST_FILE}
sed -i 's/--/-/g' ${TEST_FILE}
sed -i 's/、//g' ${TEST_FILE}
sed -i 's/|//g' ${TEST_FILE}
sed -i 's///g' ${TEST_FILE}
sed -i 's/不祥/不详/g' ${TEST_FILE}
sed -i 's/不详/未提供/g' ${TEST_FILE}
sed -i 's/为提供/未提供/g' ${TEST_FILE}
sed -i 's/未知/未提供/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)未通过\t/\1未提供\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)未提交\t/\1未提供\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)-\t/\1未提供\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)\t/\1未提供\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)—\t/\1未提供\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)不详\t/\1未提供\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)其他\t/\1未提供\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)无\t/\1无业人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)无业\t/\1无业人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)无职业\t/\1无业人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)个体经营\t/\1个体经营者\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)个体\t/\1个体经营者\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)职工\t/\1职员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)员工\t/\1职员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)离休\t/\1退(离)休人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)退休\t/\1退(离)休人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)离退休人员\t/\1退(离)休人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)离退休\t/\1退(离)休人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)退（离）休人员\t/\1退(离)休人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)已退休\t/\1退(离)休人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)退休职工\t/\1退(离)休人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)退休职员\t/\1退(离)休人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)退休人员\t/\1退(离)休人员\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)医师\t/\1医生\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)务农\t/\1农民\t/g' ${TEST_FILE}
sed -i 's/\([^\t]*\t[^\t]*\t\)自由职业\t/\1自由职业者\t/g' ${TEST_FILE}
sed -i 's/\(\([^\t]*\t\)\{4\}\)\t/\1未提供\t/g' ${TEST_FILE}
sed -i 's/\(\([^\t]*\t\)\{5\}[[:digit:]]\+\)岁/\1/g' ${TEST_FILE}
sed -i 's/\(\([^\t]*\t\)\{7\}\)已\t/\1已婚\t/g' ${TEST_FILE}
sed -i 's/\(\([^\t]*\t\)\{7\}\)\t/\1未提供\t/g' ${TEST_FILE}
sed -i 's/\(\([^\t]*\t\)\{7\}\)离异\t/\1离婚\t/g' ${TEST_FILE}
sed -i 's/\(\([^\t]*\t\)\{8\}\)汉\t/\1汉族\t/g' ${TEST_FILE}
sed -i 's/\(\([^\t]*\t\)\{8\}\)\t/\1未提供\t/g' ${TEST_FILE}
sed -i 's/\(\([^\t]*\t\)\{9\}\)\t/\1未提供\t/g' ${TEST_FILE}
sed -i 's/\(\([^\t]*\t\)\{10\}\)\t/\1未提供\t/g' ${TEST_FILE}

perl -pe 's/(\d{4})年(\d{2})月(\d{2})日(\d{2})时(\d{2})分/$1-$2-$3 $4:$5:00/g' ${TEST_FILE} \
	| perl fullconvert.pl | awk '{ gsub(/\xef\xbb\xbf/,""); print }' \
	| sort -t$'\t' -k1,1 -k14,14n -k13,13 | uniq | awk -F "\t" '{print $1 "\t" $2 "\t" $4 "\t" $6 "\t" $8 "\t" $11 "\t" $3 "\t" $10 "\t" $7 "\t" $5}' \
	>ICU_text_clean.tsv
cat ICU_text_clean.tsv ICU_Cases_clean.tsv \
	| perl basic_info_filter.pl \
	| awk '{ if($4 == "未提供") $4=""; print }' >ICU_texts.tsv

蒋有旺性别:男年龄:57岁出生地:江苏省镇江市句容市职业:其他民族:汉族
