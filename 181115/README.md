# How to use Picard for gathering Bam files
### (Picard)[] is a A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.
java -jar picard.jar GatherBamFiles \
      I=input1.bam \ 
      I=input2.bam \ 
      O=gathered_files.bam
