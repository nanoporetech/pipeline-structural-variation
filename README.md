![.](ONT_logo.png  "Oxford Nanopore Technologies")
******************
# pipeline-structural-variation

### Overview

`pipeline-structural-variation` is a pipeline for calling structural variations on whole genome sequencing data obtained from Oxford Nanopore sequencing platforms. It accepts FASTQ files and outputs aligned reads and filtered SV calls.

### Features

The pipeline performs the following steps:
- Maps reads using minimap2
- Produces QC report using NanoPlot
- Estimates appropriate parameters for variant calling depending on read depth
- Calls variants using sniffles
- Filters variants by minimum/maximum length, read support, or type (e.g. insertion, deletion, etc.)

******************
# Getting Started

### Requirements
The following software packages must be installed prior to running:

-  [miniconda3](https://conda.io/miniconda.html) - please refer to installation [instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Installation
After installing miniconda3 and snakemake (see above), install the pipeline as follows:
```bash
# Get pipeline
$ wget -O pipeline-structural-variation.tar.gz https://github.com/nanoporetech/pipeline-structural-variation/archive/v1.6.1.tar.gz
# Unzip
$ tar xvzf pipeline-structural-variation.tar.gz
# Change to directory
$ cd pipeline-structural-variation-*
# Create conda environment with all dependencies
$ conda env create -n pipeline-structural-variation -f env.yml
# Activate environment
$ conda activate pipeline-structural-variation
# To test if the installation was successful run
$ snakemake -p all
# Deactivate environment
$ conda deactivate
```

### Input

To run the pipeline the following input files are required:

| Input | Description |
|-------|-------------|
| Reference genome | FASTA file containing the reference genome (e.g. GRCh38 for human) |
| Nanopore reads| Folder containing FASTQ files or a single concatenated FASTQ file. Reads should be **q-score filtered** (see FAQ)|

### Output

 The main output files created by the pipeline are:

| Output | Description |
|--------|-------------|
| Aligned reads | Aligned reads in indexed and sorted BAM format |
| Variant calls | Called variants in VCF format |

After the pipeline has finished you can find the aligned reads in `{output_folder}/alignment/` and the indexed and zipped VCF file in `output_folder/sv_calls/{sample_name}_sniffles_filtered.vcf`.

### Usage:

To run the pipeline with default settings invoke snakemake as follows.

```bash
$ snakemake -j 30 all --config input_fastq=/data/pass/ reference_fasta=/data/ref/hg38.fa
```

`-j` specifies how many CPU cores will be used by the pipeline. `all` is the default target (see Targets); this will run all steps required for SV calling and produce a QC report for the input reads using NanoPlot. `input_fastq` specifies the input FASTQ files or a folder containing multiple input FASTQ files (e.g. the pass folder from MinKNOW).

### Targets

|Name| Description |
|--|--|
| all | Maps reads, calls variants and produces a QC report from the input reads |
| qc | Only maps reads and produces a QC report |
| call | Maps reads and calls variants but does not produce a QC report |
| eval | Evaluates the called variants against the GIAB truth set **only** applicable when sequencing [HG002](https://www.coriell.org/0/Sections/Search/Sample_Detail.aspx?Ref=NA24385&Product=DNA) |

### Options

The pipeline accepts several input parameters. They can either be changed in the `config.yml` (see below) file or specified when running snakemake.

For example:
```bash
snakemake -j 30 eval --config input_fastq=/data/pass/ reference_fasta=/data/ref/hg38.fa min_sv_length=100
```
In the above example the minimum SV length parameter was changed to 100 bp.

##### Required parameters

These parameters have to be specified to run the pipeline.

| Parameter | Allowed | Description |
|-----------|---------|-------------|
| input_fastq | Absolute file path | FASTQ file or folder containing FASTQ files |
| reference_fasta | Absolute file path | FASTA file containing the reference genome |

##### Optional parameters

| Parameter | Allowed | Default | Description |
|-----------|---------|---------|-------------|
| sample_name | string without spaces | my_sample | Name of the sample |
| min_sv_length | Integer > 40 | 50 | Minimum SV length |
| max_sv_length | Integer | 1000000 | Maximum SV length |
| min_read_length | Integer | 1000 | Min read length. Shorter reads will be ignored |
| min_read_mapping_quality | Integer | 10 | Min mapping quality. Reads with lower mapping quality will be discarded |
| min_read_support | Integer | 'auto' | Minimum read support required to call a SV (auto for auto-detect) |

# Annotating and visualising the results

## Annotate

This step will become part of the pipeline in the future. Currently, use tools like `bedtools` or `vcfanno` to annotate the VCF file with information about genes, repeats, known SV etc. from GTF/GFF/BED files.

There are many ways to retrieve annotations. A few examples are:

-  [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgTables)
-  [ENSEMBL](https://www.ensembl.org/index.html)
-  [GENCODE](https://www.gencodegenes.org/)

## Visualise

### Integrated genome viewer

Download from [here](http://software.broadinstitute.org/software/igv/).

Supported formats: BAM, BED, VCF, WIG, …

**Recommended settings for SV analysis:**

View -> Preferences -> Alignments
- Quick consensus mode on
- Hide indels < 10 bp
- Label indels > 30 bp
- Mapping quality threshold: 20

**Works well for:**

- Deletions and insertions
- Inspect reference sequence and flanking regions
- Compare to annotations
- Inspect heterozygous variants

### Ribbon

Find [here](http://genomeribbon.com/).
Supported format: BAM, BED, VCF

**Recommended settings for SV analysis:**

- Multi-read settings -> minimum number of alignments -> 4-5
- Single-read settings -> Dot plot

**Works well for:**

- Visualizing split alignments
- Inversions, translocations, duplications

******************
# Results

We benchmarked the pipeline against the preliminary [Genome in a bottle](https://jimb.stanford.edu/giab) SV truth set for HG002.

### Precision and Recall
The pipeline was run using `auto` for determining the most suitable sniffles parameters to get a good balance between precision and recall. Depending on your application you might want to change the `min_read_support` parameter to maximize either precision or recall.

| Dataset | Pipeline | Min. read support | Coverage | Precision | Recall |
|-----------|---------|-------------|-------------|-------------|------------
| HG002 q7 filtered | v1.5.0 | auto | 60* | 95.86 | 95.42 |
| HG002 q7 filtered | v1.5.0 | auto | 45* | 96.04 | 94.71 |
| HG002 q7 filtered | v1.5.0 | auto | 30* | 96.44 | 92.93 |
| HG002 q7 filtered | v1.5.0 | auto | 15* | 96.72 | 87.88 |

\* Coverage was computed using mosdepth from the BAM file. The BAM file was neither filtered by mapping quality nor read length.


### Reproducing the results
Run the pipeline with HG002 data using the `eval` target. In the output folder you will find `eval/summary.txt` containing precision and recall values.

******************
# Help
## Licence and Copyright

(c) 2018 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

### FAQs

**I only want my results! What should I do?**
Install pipeline (see Installation) and run as follows.
```bash
$ snakemake -j 16 all --config input_fastq=/data/pass/ reference_fasta=/data/ref/hg38.fa
```
You will find your results in the following locations:

**Aligned reads:**

```results/my_samples/alignment/my_sample_minimap2.bam```

**Variant calls:**

```results/my_samples/sv_calls/my_sample_sniffles_filtered.vcf.gz```

**What kind of structural variants are supported?**

Currently the pipeline has been validated to detect insertions, deletions and duplications from whole genome sequencing data. Support for inversions and translocations will be added in the future.

**What is the minimum SV length supported?**

Currently, 50 bp.

**How can i filter my reads by q-score?**

The most recent version of MinKNOW will perform filtering automatically. When using the FASTQ files from the `pass` output folder no additional filtering is required. For unfiltered datasets use [NanoFilt](https://github.com/wdecoster/nanofilt) with a minimum q-score of 6.

**How long will it take to run the pipeline for a 30X human dataset?**

When running with 30 CPU cores this will take roughly 6-8 hours.

**How much memory will I need?**

Memory consumption is determined by minimap2. Therefore, 16 GB will be required for human datasets. For smaller genomes 8 GB will be sufficient.

**How much storage will I need?**

Unzipped FASTQ files for a human 30X human dataset will require roughly 200 GB. In addition 150 GB will be required for the mapped BAM files.  

**Are FAST5 files required to run the pipeline?**

No.

**Can I use the pipeline to detect variants in a cancer dataset or to detect somatic variants?**

The pipeline has not yet been validated on cancer samples or for somatic variant detection.

**Can I use the pipeline to detect SVs for targeted sequencing experiments (e.g. long amplicons, etc.)?**

Calling variants in targeted regions is supported but requires the user to specify a BED file containing the coordinates of the targeted regions.

For example, save the following in a file called `targets.bed`

```
chr1 6579767 6589767
```

and run pipeline as follows:

```bash
snakemake -j 30 eval --config input_fastq=/data/pass/ reference_fasta=/data/ref/hg38.fa target_bed=targets.bed
```
Make sure that the chromosome names in the BED file match the names in your reference FASTA files (e.g. chr1 vs. 1).

**Can I use this pipeline to detect gene fusions using DNA data?**

The current version does not support calling gene fusions from DNA data as translocation calling is not supported yet. However, support for translocations will be available in the future.

**Can I use this pipeline to detect gene fusion using cDNA data?**

cDNA data is not supported.

**Can I run the pipeline using docker?**

Yes. To build a docker image run:
```bash
$ make build
```
Next, run the pipeline as follows:
```bash
docker run -ti -w `pwd` -v `pwd`:`pwd` pipeline-structural-variation snakemake all
```
In the future we will provide a pre built docker image.


### Abbreviations and glossary
| Term/Abbreviation | Description |
|--|--|
| SV | Structural variation |
| WGS | Whole genome sequencing |
| Target | Targets or rules are used by snakemake to define what steps of the pipeline should be executed. Changing the target can modify the behaviour of the pipeline to fit certain applications better |
| Precision | Fraction of SV calls that are present in the truth set. Computed as TP / (TP + FP). |
| Recall | Fraction of calls from the truth set that were called correctly. Computed as TP / (TP + FN).|

### References and Supporting Information

If you use this pipeline please cite:
- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191
- Sedlazeck, F. J., Rescheneder, P., Smolka, M., Fang, H., Nattestad, M., von Haeseler, A., & Schatz, M. C. (2018). Accurate detection of complex structural variations using single-molecule sequencing. Nature Methods, 15(6), 461–468. https://doi.org/10.1038/s41592-018-0001-7
- Pedersen, B. S., & Quinlan, A. R. (2018). Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics, 34(5), 867–868. https://doi.org/10.1093/bioinformatics/btx699

When using the QC report please also cite:
- De Coster, W., D’Hert, S., Schultz, D. T., Cruts, M., & Van Broeckhoven, C. (2018). NanoPack: visualizing and processing long-read sequencing data. Bioinformatics, 34(15), 2666–2669. https://doi.org/10.1093/bioinformatics/bty149

### For additional information on SV and SV calling please see:

-  [Nanopore Structural Variation Knowledge exchange](https://nanoporetech.com/resource-centre/knowledge-exchange-exploring-structural-variation-nanopore-sequencing)
-  [100 tomato genomes in 100 days](https://vimeo.com/304389956)
