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
