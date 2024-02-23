cut -f 3,4,5,6,7,8 Atha_links_timespot.tsv | sort | uniq \
	| awk '{
    count0 = 0;
    count1 = 0;
    count2 = 0;
    count3 = 0;
    for (i = a; i <= NF; i++) {
        if ($i == 0) {
            count0++;
        } else if ($i == 1) {
            count1++;
        } else if ($i == 2) {
            count2++;
        } else if ($i == 3) {
            count3++;
        }
    }
    if ((count0 >= 0 && count1 >= 0 && count2 == 0 && count3 == 0) || (count0 >= 0 && count1 == 0 && count2 >= 0 && count3 == 0)) {
        print;
    }
}'
