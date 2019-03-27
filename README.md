![.](ONT_logo.png "Oxford Nanopore Technologies")

******************

# pipeline-structural-variation

### Overview

Brief description in non-slang, explicit terms that facilitate the biological understanding of the data analysis provided by the software. Few sentences. Include biological concepts if appliable. **Link out to ReadtheDocs** specific to the repository and a explanation for further information see, `here <https://readthedocs.org/>`_

### Features

Bullet point the features that distinguish it from pre-existing or similar software, what advances have been made and how it enables the user to get what they want from the sequence data.

- Feature 1.
- Feature 2.
- Feature ...n

******************

# Getting Started

### Input

| Input | Description |
|-------|-------------|
| Reference genome | FASTA file containing the reference genome (e.g. GRCh38 for human) |
| Nanopore reads| Folder containing FASTQ files or a single concatenated FASTQ file. Reads should be **q-scores filtered** (see FAQ)|

### Output
| Output | Description |
|--------|-------------|
| Aligned reads | Aligned reads in indexed and sorted BAM format |
| Varaint calls | Called variants in VCF format |
### Dependencies

- [miniconda3](https://conda.io/miniconda.html) - install it according to the [instructions](https://conda.io/docs/user-guide/install/index.html).
- [snakemake](https://anaconda.org/bioconda/snakemake) install using `conda`.
- The rest of the dependencies are automatically installed using the `conda` feature of `snakemake`.

### Installation

After installing miniconda3 and snakemake (see dependencies), setup the pipeline as follows:
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

Any other installation method.

#### Setup conda environment manually

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
```


#### Options




******************

# Results

### Statistics and Performance

What statistics can be derived and what characterise the software i.e. false neg, false pos, what statistics models are used; boil down to hard numbers of reliability. From pacbio website: 'Performance is measured as positive predictive value (PPV); it measures TP/(TP+FP), the ratio of true positive calls over all true and false positive calls.' Describe using diagrams, tables, pictures of how to quantify a good and a bad performance...

### Data and Results

What values assess the data quality and how are these represented; explain. - **ANALYSIS** covered here

******************

# Help

### Licence and Copyright

© 2018 Oxford Nanopore Technologies Ltd.

[INSERT REPO NAME] is distributed under the terms of the Oxford Nanopore Technologies Developer licence.

### FAQs

Description of generalities of the software and user experience misunderstandings. Ref any dialogue between customer and technical services/ONT too.

### Abbreviations

Tabulate abbreviations or names or shortened text or bioinformatician terms that the biologist might not understand... remember we are trying facilitate the users' experience of the product.

### References and Supporting Information

List helpful websites, pages on the community, papers, white papers, citations and relevant internal evidence thatt might help the user.
At a very minimum a link to the BIOINFORMATICS RESOURCE IN THE COMMUNITY - IN KNOWLEDGE should be used.

#### Research Release

Research releases are provided as technology demonstrators to provide early access to features or stimulate Community development of tools. Support for this software will be minimal and is only provided directly by the developers. Feature requests, improvements, and discussions are welcome and can be implemented by forking and pull requests. However much as we would like to rectify every issue and piece of feedback users may have, the developers may have limited resource for support of this software. Research releases may be unstable and subject to rapid iteration by Oxford Nanopore Technologies.

##### Developer Release

This is a developer release provided under the terms of the Oxford Nanopore Technologies' Developer License. Developer releases are provided to allow caveated access to source code and APIs for third party tool development and exploration. Support is provided via this Github project and/or in the Community, here [https://community.nanoporetech.com/] and releases are accompanied by a release note at this website.

#### Product Release

This is a product release provided under the terms of the Oxford Nanopore Technologies' License. There is limited/no access or permission to edit the source code within this repository. Support is amply provided via this Github project and/or in the Community, here [https://community.nanoporetech.com/] and releases are accompanied by a release note at this location. We cite the Bioinformatics Resource within the Community for updates, FAQ, support, nominally.

Rev005