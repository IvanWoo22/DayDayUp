Today I try to use QIIME2 for data analysis.  
Most of them are recorded in [QIIME2_learnning_log](QIIME2_learnning_log/).

***

[phred_detect.pl](phred_detect.pl) is a program to detect whether a .fastq file is use phred 33 or 64. 
Do like below:  
```Bash
phred_detect.pl fastq.file 1000
## Here 1000 is the sample-size you want to input.
```
