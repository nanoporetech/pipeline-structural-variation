import os
import sys
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import pysam


def parse_args(argv):
    usage = "Extract reference BED from BAM file"
    parser = ArgumentParser(description=usage,
                            formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--bam',
                        dest='BAM',
                        type=str,
                        required=True,
                        help='Input BAM file')
    parser.add_argument('-f', '--filter',
                        dest='FILTER',
                        nargs="*",
                        help='Remove entries that contain')

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

    bam_file = args.BAM
    filter_names = args.FILTER
    if not os.path.exists(bam_file):
        raise OSError("Could not find {}.".format(bam_file))

    with pysam.AlignmentFile(bam_file, "r") as bam:
        for tid in range(0, bam.nreferences):
            ref_name = bam.get_reference_name(tid)
            if len(list(filter(lambda x: x in ref_name, filter_names))) > 0:
                continue

            print(ref_name,
                  "0",
                  str(bam.get_reference_length(ref_name)),
                  sep='\t')


if __name__ == '__main__':

    main()
