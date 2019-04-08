import os
import sys
import logging

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from pipeline_structural_variation.sniffles_vcf import SvFile


def parse_args(argv):
    usage = "Simple script to filter sniffles VCF files"
    parser = ArgumentParser(description=usage,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-v', '--vcf',
                        dest='VCF',
                        type=str,
                        required=True,
                        help='Input VCF file (default: stdin)')
    parser.add_argument('-o', '--output',
                        dest='OUTPUT',
                        type=str,
                        default="/dev/stdout",
                        help='Output VCF file')
    parser.add_argument('-t', "--type",
                        dest="SV_TYPE",
                        nargs="*",
                        choices=SvFile.SV_TYPES + ['all'],
                        default=None,
                        help="Only output SV of this type (default: off)")
    parser.add_argument('-m', "--min-support",
                        dest="SV_MIN_SUPPORT",
                        type=int,
                        default=0,
                        help="Min read support (default: 0)")
    parser.add_argument("--max-support",
                        dest="SV_MAX_SUPPORT",
                        type=int,
                        default=None,
                        help="Max read support (default: None)")
    parser.add_argument("--min-af",
                        dest="SV_MIN_AF",
                        type=float,
                        default=0.0,
                        help="Min SV allele frequency (default: 0)")
    parser.add_argument('-l', "--min-length",
                        dest="SV_MIN_LENGTH",
                        type=int,
                        default=50,
                        help="Min SV lengthe (default: 0)")
    parser.add_argument("--strand-support",
                        dest="STRAND_SUPPORT",
                        type=float,
                        default=0.0,
                        help="Max pvalue for strand bias filter. 0 to disable. (default: disabled)")
    parser.add_argument("--max-length",
                        dest="SV_MAX_LENGTH",
                        type=int,
                        default=None,
                        help="Max SV lengthe (default: None)")
    parser.add_argument('-c', "--chr",
                        dest="SV_CHR",
                        type=str,
                        default=None,
                        help="Only output SV only if on SV_CHR (default: off)")
    parser.add_argument("--precise",
                        dest="PRECISE",
                       	action='store_true',
                        help="Only SVs marked as precise (default: off)")
    parser.add_argument("--bed",
                        dest="BED",
                        action='store_true',
                        help="Output as BED (default: off)")
    parser.add_argument("--add-chr",
                        dest="CHR",
                        action='store_true',
                        help="Add 'chr' prefix to chromosomes (default: off)")
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

    logging.basicConfig(level=logging.INFO)

    vcf_file = args.VCF
    if not os.path.exists(vcf_file):
        raise OSError("Could not find {}.".format(vcf_file))

    sniffles_file = SvFile(vcf_file)

    sv_types = None
    if args.SV_TYPE and args.SV_TYPE[0].lower() != 'all':
        sv_types = args.SV_TYPE

    logging.info("Variants found: {}".format(sniffles_file.n_total()))
    sniffles_file.filter(min_length=args.SV_MIN_LENGTH,
                         max_length=args.SV_MAX_LENGTH,
                         min_support=args.SV_MIN_SUPPORT,
                         max_support=args.SV_MAX_SUPPORT,
                         min_af=args.SV_MIN_AF,
                         include_types=sv_types,
                         precise=args.PRECISE,
                         chr_filter=args.SV_CHR,
                         strand_support_p_value=args.STRAND_SUPPORT
                         )
    logging.info("Variants after filtering: {} ({}%)".format(sniffles_file.n_filtered(), sniffles_file.p_filtered()))

    if not args.BED:
        # sniffles_file.write_vcf("/dev/stdout")
        sniffles_file.write_vcf(args.OUTPUT)
    else:
        sniffles_file.write_bed(args.OUTPUT)


if __name__ == '__main__':

    main()
