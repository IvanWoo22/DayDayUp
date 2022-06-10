PREFIX=/share/home/wangq/data/Pseudomonas
ID=IPR004361
NAME="GlyI"

grep -h "${ID}" ${PREFIX}/STRAINS/Pseudom_aeru_PAO1/*.fa.tsv | cut -f 1,4,5

## Local run and upload to HPCC.
aria2c -x 12 https://ftp.ncbi.nlm.nih.gov/hmm/current/hmm_PGAP.HMM/TIGR00068.1.HMM -o hmm/${NAME}.tigrfam.hmm

