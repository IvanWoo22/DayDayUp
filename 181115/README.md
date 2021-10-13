# How to use Picard for gathering Bam files

### [Picard](http://broadinstitute.github.io/picard/) is a set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.

Concatenate one or more BAM files as efficiently as possibleThis tool performs a rapid "gather" operation on BAM files
after scatter operations where the same process has been performed on different regions of a BAM file creating many
smaller BAM files that now need to be concatenated (reassembled) back together.

Assumes that the list of BAM files provided as INPUT are in the order that they should be concatenated and simply
concatenates the bodies of the BAM files while retaining the header from the first file. Operates via copying of the
gzip blocks directly for speed but also supports generation of an MD5 on the output and indexing of the output BAM file.
Only supports BAM files, does not support SAM files.

## Usage example:

```bash
java -jar picard.jar GatherBamFiles \
      I=input1.bam \ 
      I=input2.bam \ 
      O=gathered_files.bam
```
