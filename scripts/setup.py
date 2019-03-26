"""
HIVE setup script

Copyright (c) 2018 by Oxford Nanopore Technologies Ltd.
"""
from sv_calling_glue import __version__
from setuptools import setup, find_packages


setup(
    name='sv_calling_glue',
    version=__version__,
    author='metrichor-bio',
    setup_requires=['pytest-runner'],
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
            'sniffles-filter = sv_calling_glue.sv_filter_sniffles:main',
            'sniffles-telemetry = sv_calling_glue.sniffles_telemetry:main',
            'sniffles-stats = sv_calling_glue.sniffles_stats:main',
            'sniffles-edit = sv_calling_glue.sniffles_edit:main',
            'sniffles-strand = sv_calling_glue.sniffles_strand:main',
            'sniffles-sample = sv_calling_glue.vcf_downsample:main',
            'bamref2bed = sv_calling_glue.sniffles_bamref2bed:main',
            'cat_fastq = sv_calling_glue.cat_fastq:main'
        ]
    }
)
