Boilerplate for building pipelines using snakemake
==================================================

Snakemake is a text-based workflow system which can be used to build robust bioinformatics pipelines. Similarly to GNU make a list of rules define how to create output files from input files. The dependencies between the rules are inferred automatically. The snakemake rule description language can be mixed with Python code, which makes it quite practical. To learn more about the framework check out the [snakemake wiki](https://bitbucket.org/snakemake/snakemake/wiki/browse/).

This pipeline skeleton contain a minimal layout with some utility targets. To keep things simple it does not try to package the pipeline as a python 
module but files can be included from `lib/`.

Usage
-----

Edit `config.yml` and issue snakemake <target> to invoke the target of your choice. Issue `snakemake help` to get a list of all targets and their descriptions.

Python dependencies
-------------------

See `requirements.txt` for other dependencies. Use `pip install -r requirements.txt` to install them.

Application dependencies
------------------------

Layout
------

* `README.md`
* `Snakefile`         - master snakefile
* `config.yml`        - YAML configuration file
* `snakelib/`         - snakefiles collection included by the master snakefile
* `lib/`              - python files included by analysis scripts and snakefiles
* `scripts/`          - analysis scripts
* `data/`             - input data needed by pipeline - use with caution to avoid bloated repo
* `results/`          - pipeline results to be commited - use with caution to avoid bloated repo
* `requirements.txt`  - list of python package dependencies

Useful snakemake targets
------------------------

```
help                    list all targets and descriptions
info                    print pipeline information
clean_workdir           delete working directory. WARNING: all data will be lost!
clean_resdir            delete results directory. WARNING: all data will be lost!
```
