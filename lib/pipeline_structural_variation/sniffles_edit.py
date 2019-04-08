from __future__ import print_function

import os
import sys
import tempfile
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from pipeline_structural_variation.sniffles_vcf import SvFile


def parse_args(argv):
    usage = "Simple script to edit sniffles VCF files"
    parser = ArgumentParser(description=usage,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--vcf',
                        dest='VCF',
                        type=str,
                        required=True,
                        help='Input VCF file')
    parser.add_argument('-o', '--output',
                        dest='OUTPUT',
                        type=str,
                        default="/dev/stdout",
                        help='Output VCF file')
    parser.add_argument("--vcf-version",
                        dest="VCF_VERSION",
                        action='store_true',
                        help="Change VCF version to 4.2")
    parser.add_argument("--check",
                        dest="CHECK",
                        action='store_true',
                        help="Check file for integrity.")
    parser.add_argument("--ins-length",
                        dest="INS_LENGTH",
                        action='store_true',
                        help="Set INS length to 1")

    args = parser.parse_args(argv)
    return args


def fix_header(vcf_file):
    f = tempfile.NamedTemporaryFile("w", delete=False)
    with open(vcf_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if "##fileformat=" in line:
                print("##fileformat=VCFv4.2", file=f)
            else:
                print(line, file=f)
    f.close()

    return f.name


def check_sniffles_file(sniffles_file):
    sniffles_file.check()


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

    if args.VCF_VERSION:
        vcf_file = fix_header(vcf_file)

    sniffles_file = SvFile(vcf_file, fix_ins=args.INS_LENGTH)

    # sniffles_file.fix_overlaping_ins_dup()

    check_sniffles_file(sniffles_file)

    sniffles_file.write_vcf(args.OUTPUT)


if __name__ == '__main__':

    main()
