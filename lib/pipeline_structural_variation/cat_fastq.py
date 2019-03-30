from __future__ import print_function

import os
import sys
import glob
import argparse
import pysam
import logging


def parse_args(argv):
    """
    Commandline parser

    :param argv: Command line arguments
    :type argv: List
    """
    usage = "Cat long lists of FASTQ files"
    parser = argparse.ArgumentParser(
        description=usage,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-l", "--log",
                        dest="log",
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR',
                                 'CRITICAL',
                                 'debug', 'info', 'warning', 'error',
                                 'critical'],
                        default="INFO",
                        help="Print debug information")
    parser.add_argument("-o", "--output",
                        dest="OUT",
                        type=str,
                        default="/dev/stdout",
                        help="Output file. (default: stdout)")

    parser.add_argument("FASTQ",
                        nargs="+",
                        type=str,
                        help="FASTQ files or folders containing FASTQ files")

    args = parser.parse_args(argv)


    return args


def find_file_in_folder(folder, pattern="*.fastq"):
    if os.path.isfile(folder):
        return folder
    files = []
    for file in glob.glob(os.path.join(folder, pattern)):
       files.append(file)

    if len(files) == 0:
        logging.warning("Could not find {} files in {}".format(pattern, folder))

    return files


def parse_fastqs(filename, fout):
    with pysam.FastxFile(filename) as fh:
        for entry in fh:
            if entry.comment:
                entry.comment = "CO:Z:{}".format(entry.comment)

            fout.write(str(entry) + "\n")


def get_file_names(path):
    if os.path.exists(path):
        filenames = [path]
        if os.path.isdir(path):
            logging.info("Searching {} for FASTQ files".format(path))
            filenames = find_file_in_folder(path)
    else:
        logging.warning("Could not find {}".format(path))

    return filenames


def format_fq(paths, out_filename):
    """
    Concatenate FASTQ files

    :param paths: Input FASTQ files or folders containing FASTQ files
    :param out_filename: Output FASTQ file
    :return: None
    """

    with open(out_filename, mode='w') as fout:
        for path in paths:
            for filename in get_file_names(path):
                parse_fastqs(filename, fout)


def main(argv=sys.argv[1:]):
    """
    Basic command line interface to telemap.

    :param argv: Command line arguments
    :type argv: list
    :return: None
    :rtype: NoneType
    """
    args = parse_args(argv=argv)

    numeric_level = getattr(logging, args.log.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % args.log.upper())
    logging.basicConfig(level=numeric_level, format='%(message)s')

    format_fq(args.FASTQ, args.OUT)


if __name__ == '__main__':
    main()
