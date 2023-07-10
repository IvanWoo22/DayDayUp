curl -fsSL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/309/985/GCF_000309985.2_CAAS_Brap_v3.01/GCF_000309985.2_CAAS_Brap_v3.01_genomic.gff.gz |
  pigz -dc |
  awk '$1~"^NC_024"&&$3=="gene"{if($7=="+"){print $1 "\t" $9 "\t" $4 "\t" $5}
    else{print $1 "\t" $9 "\t" $5 "\t" $4}}' |
  perl -lane '$F[1] =~ /Name=(\w+)/;
	    $F[1] = $1;
	    print join("\t", @F)' \
    >Data/Brap.gff.tsv
