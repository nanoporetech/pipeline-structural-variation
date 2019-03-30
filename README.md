
![.](ONT_logo.png "Oxford Nanopore Technologies")

******************

# pipeline-structural-variation

### Overview


Brief description in non-slang, explicit terms that facilitate the biological understanding of the data analysis provided by the software. Few sentences. Include biological concepts if appliable. **Link out to ReadtheDocs** specific to the repository and a explanation for further information see, `here <https://readthedocs.org/>`_

### Features

The pipeline performs the following steps: 
- Maps reads using minimap2
- Produces QC report using NanoPlot
- Estimates appropriate parameters for variant calling depending on input coverage
- Calls variants using sniffles

******************

# Getting Started

### Input

To run the pipeline the following input files are required:

| Input | Description |
|-------|-------------|
| Reference genome | FASTA file containing the reference genome (e.g. GRCh38 for human) |
| Nanopore reads| Folder containing FASTQ files or a single concatenated FASTQ file. Reads should be **q-scores filtered** (see FAQ)|

### Output

The main output files created by the pipeline are:

| Output | Description |
|--------|-------------|
| Aligned reads | Aligned reads in indexed and sorted BAM format |
| Variant calls | Called variants in VCF format |

### Dependencies

To run the pipeline the following software packages have to be installed on your system:

- [miniconda3](https://conda.io/miniconda.html) - install it according to the [instructions](https://conda.io/docs/user-guide/install/index.html).
- [snakemake](https://anaconda.org/bioconda/snakemake) install using `conda`.
- The rest of the dependencies are automatically installed using the `conda` feature of `snakemake`.

### Installation

After installing miniconda3 and snakemake (see dependencies), install the pipeline as follows:
```bash
# Get pipeline
$ git clone https://github.com/nanoporetech/pipeline-structural-variation.git
# Unzip

$ conda env create -n pipeline-structural-variation -f env.yml
$ conda activate pipeline-structural-variation
# To test if the installation worked run
$ snakemake --use-conda  -p all

```

### Usage: 



```bash

```

#### Run using docker

Not available yet

#### Setup conda environment manually

If you want to setup your conda environment manually follow the instructions below, activate the "environment" before running snakemake and skip the `--use-conda` parameter.

```bash
# Get pipeline
$ git clone https://github.com/nanoporetech/pipeline-structural-variation.git
# Unzip
# Create conda environment
$ conda env create -n pipeline-structural-variation -f env.yml
# Activate environment
$ conda activate pipeline-structural-variation
# To test if the installation worked run (don't use --use-conda)
$ snakemake -p all
# Deactivate environment
# conda deactivate

```


#### Options

The pipeline accepts several input parameters. They can either be changed in the `config.yml` (see below) file ore specified when running snakemake.
For example:
```bash
snakemake -j 30 eval --config input_fastq=/data/pass/ reference_fasta=/data/ref/hg38.fa
``` 

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
| threads | Integer number | 1 | Threads to use when running analysis |
| min_sv_length | Integer > 40 | 50 | Minimum SV length | 
| max_sv_length | Integer | 1000000 | Maximum SV length | 
| min_read_length | Integer | 1000 | Min read length. Shorter reads will be ignored |
| min_read_mapping_quality | Integer | 10 | Min mapping quality. Reads will lower mapping quality will be discarded |
| min_read_support | Integer | 'auto' | Minimum read support required to call a SV (auto for auto-detect) |


# What to do with the results?

## Visualise

### IGV



******************

# Results

### Statistics and Performance

What statistics can be derived and what characterise the software i.e. false neg, false pos, what statistics models are used; boil down to hard numbers of reliability. From pacbio website: 'Performance is measured as positive predictive value (PPV); it measures TP/(TP+FP), the ratio of true positive calls over all true and false positive calls.' Describe using diagrams, tables, pictures of how to quantify a good and a bad performance...

### Data and Results

What values assess the data quality and how are these represented; explain. - **ANALYSIS** covered here

******************

# Help

## Licence and Copyright

(c) 2018 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

### FAQs

**What kind of structural variants are supported?**

Currently the pipeline has been validate for detection insertions, deletions and duplications from whole genome sequencing data. Support for inversions and translocation will be added in the future.  

**What is the minimum SV length supported?**

Currently, 50 bp.

**How can i filter my reads by q-scores?**

Most recent version of Minknow will perform filtering automatically. When using the FASTQ files from the `pass` output folder no additional filtering is required. For unfiltered datasets use [NanoFilt](https://github.com/wdecoster/nanofilt) with a minimum q-score of 6.

**How long will it take to run the pipeline for a 30X human dataset?**

When running with 30 CPU cores roughly 6-8 hours.

**How much memory will I need?**

Memory consumption is determined by minimap2. Therefore, 16 GB will be required for human datasets. For smaller genomes 8 GB will be sufficient. 

Unzipped FASTQ files for a human 30X human dataset will require roughly 200 GB. In addition 150 GB will be required for the mapped BAM files.  

**Are FAST5 files required to run the pipeline?**

No

**Can I use the pipeline to detect variants in a cancer dataset?**

The pipeline has not been validated yet on cancer samples.

**Can I use the pipeline to detect SVs for targeted sequencing experiments (e.g. long amplicons, etc.)?**

Calling variants in targeted regions is support but requires the user to specify a bed file containing the coordinates of the targeted regions.
For example save the following in a file called `targets.bed`
```
chr1    6579767 6589767
```
and run pipeline as follows:
```bash
snakemake -j 30 eval --config input_fastq=/data/pass/ reference_fasta=/data/ref/hg38.fa target_bed=targets.bed
```
Make sure that the chromosome names in the BED file match the names in your reference FASTA files (e.g. chr1 vs. 1)

**Can I use this pipeline to detect gene fusion using DNA data?**
The current version does not support calling gene fusions from DNA data. However, support for calling translocations will be available in the future. 

**Can I use this pipeline to detect gene fusion using cDNA data?**
cDNA data is currently not supported

### Abbreviations

| Abbreviation | Description |
|--|--|
| SV | Structural variation |
| WGS | Whole genome sequencing |

### References and Supporting Information

If you use this pipeline please cite:

- Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191
- Sedlazeck, F. J., Rescheneder, P., Smolka, M., Fang, H., Nattestad, M., von Haeseler, A., & Schatz, M. C. (2018). Accurate detection of complex structural variations using single-molecule sequencing. Nature Methods, 15(6), 461–468. https://doi.org/10.1038/s41592-018-0001-7
- Pedersen, B. S., & Quinlan, A. R. (2018). Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics, 34(5), 867–868. https://doi.org/10.1093/bioinformatics/btx699


##### For additional information on SV and SV calling please see:

- [Nanopore Structural Variation Knowledge exchange](https://community.nanoporetech.com/posts/4192)
- [100 tomato genomes in 100 days](https://vimeo.com/304389956)

