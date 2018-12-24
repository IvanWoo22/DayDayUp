# Stop Sites Work Flow
## Quality Control
### Fastp
We have PE150 fastqc files, which contain a few adapters.  We use [fastp](https://github.com/OpenGene/fastp) to cut adapters and quality control.  
```bash
fastp --detect_adapter_for_pe --correction \
  -i S036_hongxu-A_10X-rRNA_BHFV2MDSXX_S20_L002_R1_001.fastq.gz \
  -I S036_hongxu-A_10X-rRNA_BHFV2MDSXX_S20_L002_R2_001.fastq.gz \
  -o 10Xr_out.R1.fq.gz \
  -O 10Xr_out.R2.fq.gz
```


### Bowtie2
We use 18s and 28s rRNA sequences to build index.
```bash
bowtie2-build 28s.fasta 28s

bowtie2 -p 2 --local -x 28s -1 10Xr_out.R1.fq.gz -2 10Xr_out.R2.fq.gz -S 10Xr_28s.sam
```
