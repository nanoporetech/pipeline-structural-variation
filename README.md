Boilerplate for building pipelines using snakemake
==================================================

Usage
-----


Python dependencies
-------------------

See `requirements.txt` for other dependencies. Use `pip install -r requirements.txt` to install them.

Application dependencies
------------------------

Layout
------

* `README.md`
* `Snakefile`         - master snakefile
* `config.mk`         - makefile for general configuration
* `analysis.mk`       - makefile for pipeline steps
* `utils.mk`          - utility makefile
* `scripts/`          - analysis scripts
* `lib/`              - python files included by analysis scripts and snakefiles
* `data/`             - input data needed by pipeline - use with caution to avoid bloated repo
* `results/`          - pipeline results to be commited - use with caution to avoid bloated repo
* `requirements.txt`  - list of python package dependencies

Useful snakemake targets
------------------------

```
help                    list all targets and descriptions
```
