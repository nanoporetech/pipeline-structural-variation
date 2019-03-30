import os
import sys
import logging
import base64
import zlib
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import pysam

from sv_calling_glue.sniffles_telemetry import get_summary_from_vcf, \
    compute_stats_from_bam
from sv_calling_glue.sniffles_vcf import SvFile
from sv_calling_glue.telemetry import Telemetry


def parse_args(argv):
    usage = "Simple script to filter sniffles VCF files"
    parser = ArgumentParser(description=usage,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--vcf',
                        dest='VCF',
                        type=str,
                        required=True,
                        help='Input VCF file')

    args = parser.parse_args(argv)
    return args


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to cuecat.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    vcf_file = args.VCF
    if not os.path.exists(vcf_file):
        raise OSError("Could not find {}.".format(vcf_file))

    print(get_summary_from_vcf(vcf_file))
    # print(compute_stats_from_bam(bam_file))

        # record.set('stats', compute_stats_from_bam(bam_file))
        # record.set('calls', get_calls_from_vcf(vcf_file))
        # record.set('summary', get_summary_from_vcf(vcf_file))
        # record.set('vcf', base_encode64(vcf_file))


if __name__ == '__main__':

    main()
