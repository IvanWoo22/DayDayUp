# Divide LSC evenly into four/six/ten parts

We use these three Perl scripts to divide LSC region into four/six/ten parts. The original table looks as follows:

```
Vanilla	Van_planifolia	Target	NC_026778	148011	86359-116163,118207-148011	1-86358	116164-118206	29805	29805	86358	2043
Schizaea	Schizaea_pectinata	Target	NC_035808	156392	73223-113315,116297-156389	1-73222,156390-156392	113316-116296	40093	40093	73225	2981
Paphiopedilum	Paph_armeniacum	Target	NC_026779	162682	91735-125372,129045-162682	1-91734	125373-129044	33638	33638	91734	3672
Plantago	Plan_media	Target	NC_028520	164130	82766-121137,125751-164122	1-82765,164123-164130	121138-125750	38372	38372	82773	4613
Acacia	Aca_ligulata	Target	NC_026134	174233	92779-130995,135981-174197	1-92778,174198-174233	130996-135980	38217	38217	92814	4985
```

Just type in:

```Bash
perl split_4.pl lsc.tsv > lsc_4.tsv
perl split_6.pl lsc.tsv > lsc_6.tsv
perl split_10.pl lsc.tsv > lsc_10.tsv
```

and then you will get three new tables with divided parts.
