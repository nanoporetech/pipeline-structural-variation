"""
HIVE setup script

Copyright (c) 2018 by Oxford Nanopore Technologies Ltd.
"""
from pipeline_structural_variation import __version__
from setuptools import setup, find_packages


setup(
    name='pipeline_structural_variation',
    version=__version__,
    author='metrichor-bio',
    # setup_requires=['pytest-runner'],
    description='SV calling workflow python glue',
    zip_safe=False,
    install_requires=[
        'pysam==0.14.1',
        'pybedtools',
        'scipy',
        'tqdm'
    ],
    packages=find_packages(exclude=("tests",)),

    entry_points={
        "console_scripts": [
            'sniffles-filter = pipeline_structural_variation.sv_filter_sniffles:main',
            'sniffles-telemetry = pipeline_structural_variation.sniffles_telemetry:main',
            'sniffles-stats = pipeline_structural_variation.sniffles_stats:main',
            'sniffles-edit = pipeline_structural_variation.sniffles_edit:main',
            'sniffles-strand = pipeline_structural_variation.sniffles_strand:main',
            'sniffles-sample = pipeline_structural_variation.vcf_downsample:main',
            'bamref2bed = pipeline_structural_variation.sniffles_bamref2bed:main',
            'cat_fastq = pipeline_structural_variation.cat_fastq:main'
        ]
    }
)
