# Trio splicing visualization and analysis

Python class for the visualization (sashimi and coverage plots) of splicing events in a trio (proband, mother, father).
Junctions that differ greatly between the proband and the mother/father will be marked.

/// ---------------------------------------- ///

## OVERVIEW:

* For each sample, aligned reads are used to find splicing junctions and their coverage. A plot of coverage and splicing events is then generated as well as a summary table of junctions.

* If a trio vcf file is provided, variants that can distinguish father- and mother-derived reads are found. These are then used for generating additional coverage and splicing plots for the proband.

* Splicing junction coverage is normalized to the library size-normalized number of reads assigned to the gene of interest (i.e. the larger the library and the more reads are assigned to the gene, the more likely a splicing junction will be detected).

* If the option is selected (default=False), a kneedle algorithm is used to define a threshold of low coverage for filtering junctions. Otherwise, thresholds are set to 0.

* Proband splicing junctions are then compared to those of the father and mother as follows:

    * Junctions that are detected above threshold in the proband but not in both parents are marked as "proband_unique".

    * Junctions that are detected above threshold in both parents but not in the proband are marked as "proband_missing".

    * A linear regression is fit on the junctions coverage of the proband and each parent (individually), under the assumption that most junctions will be used similarly, except for those possibly impacted by a variant (e.g. SNP affecting a splicing donor site). Any junction with a difference in coverage > 2 standard deviations from the fit are marked as outliers.

* A plot visualizing junctions flagged as described above will be generated. If a tsv file detailing known variants of interest is provided, these will be plotted as well for reference. This can be useful to prioritize junctions to look at in more detail.

/// ---------------------------------------- ///

## USAGE:

### OPTIONS:

For running the analysis, use the wrapper.py script, specifying the following options:

* **--genes**\
Ensembl IDs of genes of interest, comma separated.

* **--gtf**\
Path to the GTF annotation needed to identify exons for each gene of interest.\
Can be gzipped.

* **--bam**\
Paths to the aligned bam files, comma separated.\
For each file, a bai index should be present.

* **--relationships**\
Unique IDs for family members, comma separated (e.g. "proband,mother,father").\
One sample must be named "proband".\
Also, these are assumed to be in the same order as the bam files.

* **--vcf**\
Optional path to the trio vcf file, gzipped. Sample columns are assumed to be in the same order as the sample IDs above.

* **--interesting_variants**\
Optional path to a tab-separated file detailing known variants of interest to be plotted.
Must have the following columns:\

    * **contig**\
    Contig in the same format as in the GTF file.

    * **pos**\
    Position in basepairs of the variants.

    * **flag**\
    Name to be used for the variant during plotting.\
    Can be repeated (e.g. to group all variants of the same kind).

* **--low_coverage_filter**\
Optional parameter to determine a coverage threshold for splice junctions filtering.

### EXAMPLE COMMAND:

```
	python wrapper.py \
	--genes ENSG00000000123,ENSG00000000456 \
	--gtf /path/to/annotations.gtf \
	--bam /path/to/sample1.bam,/path/to/sample2.bam,/path/to/sample3.bam \
	--relationships proband,parent1,parent2 \
	--vcf /path/to/variants.vcf.gz \
	--interesting_variants /path/to/variants/of/interest.tsv \
	--low_coverage_filter
```

/// ---------------------------------------- ///

## OUTPUTS:

### FOLDER STRUCTURE:

For each gene, a folder (named as the gene) is created for the results/plots and has the following structure:

<pre>
<b>ENSG00000000123</b>
│
├── <b>proband</b>
│   │  Folder containing data for the proband.
│   │
│   ├── <b>exons.tsv</b>
│   │   Coordinates of exons extracted from the GTF file.
│   │
│   ├── <b>junctions_origin=XXX.tsv</b>
│   │   Table of detected junctions using reads of different origin ("unknown" = all reads were used, regardless of origin).
│   │   If a trio vcf was provided, files showing junctions detected using paternal- or maternal-derived reads are generated.
│   │
│   ├── <b>sashimi_origin=XXX.png</b>
│   │   Sashimi and coverage plot of junctions detected using reads of different origin ("unknown" = all reads were used, regardless of origin).
│   │   If a trio vcf was provided, files showing junctions detected using paternal- or maternal-derived reads are generated.
│   │
│   └── <b>sashimi-cleaned_origin=unknown.png</b>
│       Same as sashimi plot above, but only junctions of known exons are plotted.
│    
├── <b>parent_1_id</b>
│   Folder containing data for parent_1.
│
├── <b>parent_2_id</b>
│   Folder containing data for parent_2.
│
├── <b>family_junctions.tsv</b>
│   Summary table of all detected junctions in the trio.
│   Junctions of interest are flagged in the appropriate column.
│
├── <b>reads_origin_stats.tsv</b>
│   Number of reads detected, divided based on origin.
│
├── <b>XXX_linregress.png</b>
│   Linear regression fit of junctions coverage between two samples.
│   Outlier junctions are marked in red.
│
├── <b>proband_flagged_junctions.png</b>
│   Plot of flagged junctions.
│   If variants of interest were specified, they'll be plotted here.
│
├── <b>sashimi_all_samples.png</b>
│   Sashimi and coverage plot of junctions detected showing all trio samples.
│
└── <b>sashimi-cleaned_all_samples.png</b>
    Same as sashimi plot above, but only junctions of known exons are plotted.
</pre>

### EXAMPLE RESULTS:

#### Trio sashimi plot

<p align="center">
  <img src="https://github.com/mmarchetti90/splicing_check/blob/main/images/sashimi-cleaned_all_samples.png">
</p>

#### Linear regression and outlier junctions

<p align="center">
  <img src="https://github.com/mmarchetti90/splicing_check/blob/main/images/proband_mother_linregress.png">
</p>

#### Flagged junctions

<p align="center">
  <img src="https://github.com/mmarchetti90/splicing_check/blob/main/images/proband_flagged_junctions.png">
</p>

/// ---------------------------------------- ///

## DEPENDENCIES:

Python 3.9.16 &

	gzip
	matplotlib
	numpy
	pandas
	pysam
	seaborn
	scipy

/// ---------------------------------------- ///

## NOTES:

* Work in progress, under testing.
