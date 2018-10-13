# fastq-join
## With usearch v11, I turnned 2 ```.fastq``` files, which are  a forward seq and a reverse seq of the same sample, into a ```.fastq``` file covering reads of the two files.  
### Basic usage  
The simplest way to use fastq_mergepairs is to specify the the forward and reverse FASTQ filenames and an output FASTQ filename.  
```
usearch -fastq_mergepairs SampleA_R1.fastq -reverse SampleA_R2.fastq -fastqout merged.fq  
```

###  Automatic R2 filename  
If the -reverse option is omitted, the reverse FASTQ filename is constructed by replacing R1 with R2. The following command line is equivalent to the example above.  
```
usearch -fastq_mergepairs SampleA_R1.fastq -fastqout merged.fq
```

### Merging multiple FASTQ file pairs in a single command
You can specify two or more FASTQ filenames following -fastq_mergepairs. In the following example, SampleA and SampleB are both merged. The R2 filenames are constructed automatically as explained above, or can be given explicitly using the -reverse option.
```
usearch -fastq_mergepairs SampleA_R1.fastq SampleB_R1.fastq -fastqout merged.fq
```
