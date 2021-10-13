# Do a microbial composition analysis with data from [EBI](https://www.ebi.ac.uk/ena/data/view/PRJEB8035)

## Use QIIME2 to analyse.

Before beginning this tutorial, create a new directory and change to that directory.

```Bash
mkdir qiime2
cd qiime2
```

### Download data about t5koday93.

Here are some of links:
> [T5KONAM1Water93](ftp.sra.ebi.ac.uk/vol1/ERA395/ERA395163/fastq/seqs_T5KONAM1Water93.fastq.gz)  
> [T5KONAM2Water93](ftp.sra.ebi.ac.uk/vol1/ERA395/ERA395163/fastq/seqs_T5KONAM2Water93.fastq.gz)  
> [T5KONAM3Water93](ftp.sra.ebi.ac.uk/vol1/ERA395/ERA395163/fastq/seqs_T5KONAM3Water93.fastq.gz)   
> ......

And more you can get on [PRJEB8035](https://www.ebi.ac.uk/ena/data/view/PRJEB8035).  
After download all files we need, remove them to ```qiime/t5koday93/```.

### Create a sample_metadata.tsv and t5koday93-manifest.

In this part we can take [“Fastq manifest” formats](https://docs.qiime2.org/2018.8/tutorials/importing/)
and [Sample metadata](https://docs.qiime2.org/2018.8/tutorials/moving-pictures/) as an example.

### Import data.

```Bash
qiime tools import\
  --type 'SampleData[SequencesWithQuality]' \
  --input-path t5koday93-manifest \
  --output-path demux.qza \
  --input-format SingleEndFastqManifestPhred33
```

### Generate a summary of the demultiplexing results.

```Bash
qiime demux summarize \
    --i-data demux.qza \
    --o-visualization demux.qzv
```

### Sequence quality control and feature table construction

#### Option 1: DADA2

[DADA2](https://www.ncbi.nlm.nih.gov/pubmed/27214047) is a pipeline for detecting and correcting (where possible)
Illumina amplicon sequence data. As implemented in the ```q2-dada2``` plugin, this quality control process will
additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the
sequencing data, and will filter chimeric sequences.  
In the ```demux.qzv``` quality plots, we see that the quality of all the bases seems to be high, so we won’t trim any
bases from the sequences. We’ll truncate our sequences at 253 bases. This next command may take up to 10 minutes to run,
and is the slowest step in this tutorial.

```Bash 
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 0 \
  --p-trunc-len 253 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats stats-dada2.qza
```

Make it visualized:

```Bash
qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv
```

If you’d like to continue the tutorial using this FeatureTable (opposed to the Deblur feature table generated in Option
2), run the following commands.

```Bash 
mv rep-seqs-dada2.qza rep-seqs.qza
mv table-dada2.qza table.qza
```

#### Option 2: Deblur

[Deblur](http://msystems.asm.org/content/2/2/e00191-16) uses sequence error profiles to associate erroneous sequence
reads with the true biological sequence from which they are derived, resulting in high quality sequence variant data.
This is applied in two steps. First, an initial quality filtering process based on quality scores is applied. This
method is an implementation of the quality filtering approach described by Bokulich et al. (2013).

```Bash  
qiime quality-filter q-score \
 --i-demux demux.qza \
 --o-filtered-sequences demux-filtered.qza \
 --o-filter-stats demux-filter-stats.qza
```

Next, the Deblur workflow is applied using the ```qiime deblur denoise-16S method```.

```Bash
qiime deblur denoise-16S \
  --i-demultiplexed-seqs demux-filtered.qza \
  --p-trim-length 253 \
  --o-representative-sequences rep-seqs-deblur.qza \
  --o-table table-deblur.qza \
  --p-sample-stats \
  --o-stats deblur-stats.qza
```

```Bash
qiime metadata tabulate \
  --m-input-file demux-filter-stats.qza \
  --o-visualization demux-filter-stats.qzv
qiime deblur visualize-stats \
  --i-deblur-stats deblur-stats.qza \
  --o-visualization deblur-stats.qzv
```  

### FeatureTable and FeatureData summaries

After the quality filtering step completes, you’ll want to explore the resulting data.

```Bash
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
```  

### Generate a tree for phylogenetic diversity analyses

QIIME supports several phylogenetic diversity metrics, including Faith’s Phylogenetic Diversity and weighted and
unweighted UniFrac.

```Bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

### Alpha and beta diversity analysis

* Alpha diversity
    * Shannon’s diversity index (a quantitative measure of community richness)
    * Observed OTUs (a qualitative measure of community richness)
    * Faith’s Phylogenetic Diversity (a qualitiative measure of community richness that incorporates phylogenetic
      relationships between the features)
    * Evenness (or Pielou’s Evenness; a measure of community evenness)
* Beta diversity
    * Jaccard distance (a qualitative measure of community dissimilarity)
    * Bray-Curtis distance (a quantitative measure of community dissimilarity)
    * unweighted UniFrac distance (a qualitative measure of community dissimilarity that incorporates phylogenetic
      relationships between the features)
    * weighted UniFrac distance (a quantitative measure of community dissimilarity that incorporates phylogenetic
      relationships between the features)

```Bash
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 2188 \
  --m-metadata-file sample-metadata.tsv \
  --output-dir core-metrics-results
```  

The sampling depth of 2188 chosen based on the DADA2 feature table summary. If you are using a Deblur feature table
rather than a DADA2 feature table, you might want to choose a different even sampling depth. Apply the logic from the
previous paragraph to help you choose an even sampling depth.

We’ll first test for associations between categorical metadata columns and alpha diversity data. We’ll do that here for
the Faith Phylogenetic Diversity (a measure of community richness) and evenness metrics.

```Bash 
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization core-metrics-results/evenness-group-significance.qzv
```

Next we’ll analyze sample composition in the context of categorical metadata using PERMANOVA (first described in
Anderson (2001)) using the ```beta-group-significance``` command.

```Bash 
qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file sample-metadata.tsv \
  --m-metadata-column treatment-group \
  --o-visualization core-metrics-results/unweighted-unifrac-treatment-significance.qzv \
  --p-pairwise
```

### Alpha rarefaction plotting

In this section we’ll explore alpha diversity as a function of sampling depth using the qiime diversity
alpha-rarefaction visualizer.

```Bash
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 3000 \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization alpha-rarefaction.qzv
```

> The value that you provide for --p-max-depth should be determined by reviewing the “Frequency per sample” information presented in the table.qzv file that was created above. In general, choosing a value that is somewhere around the median frequency seems to work well, but you may want to increase that value if the lines in the resulting rarefaction plot don’t appear to be leveling out, or decrease that value if you seem to be losing many of your samples due to low total frequencies closer to the minimum sampling depth than the maximum sampling depth.

The visualization will have two plots. The top plot is an alpha rarefaction plot, and is primarily used to determine if
the richness of the samples has been fully observed or sequenced. The bottom plot in this visualization is important
when grouping samples by metadata.

### Taxonomic analysis

In the next sections we’ll begin to explore the taxonomic composition of the samples, and again relate that to sample
metadata.   
The first step in this process is to assign taxonomy to the sequences in our ```FeatureData[Sequence]``` QIIME 2
artifact.  
We’ll do that using a pre-trained Naive Bayes classifier and the q2-feature-classifier plugin.  
This classifier was trained on the Greengenes 13_8 99% OTUs, where the sequences have been trimmed to only include 250
bases from the region of the 16S that was sequenced in this analysis (the V4 region, bound by the 515F/806R primer pair)
. We’ll apply this classifier to our sequences, and we can generate a visualization of the resulting mapping from
sequence to taxonomy.

```Bash
wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/2018.8/common/gg-13-8-99-515-806-nb-classifier.qza"
```

```Bash
qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```

Next, we can view the taxonomic composition of our samples with interactive bar plots. Generate those plots with the
following command and then open the visualization.

```Bash
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file sample-metadata.tsv \
  --o-visualization taxa-bar-plots.qzv
```

### Cluster samples using UPGMA.

In addition to using PCoA, samples were clustered using UPGMA (unweighted pair group method with arithmetic mean).

```Bash
qiime diversity beta-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-metric  unweighted_unifrac \
  --p-clustering-method  upgma \
  --m-metadata-file sample-metadata.tsv  \
  --p-sampling-depth 2188 \
  --output-dir beta_tree
  ```

## More info on [QIIME2_docs](https://docs.qiime2.org/2018.8/tutorials/moving-pictures/)
