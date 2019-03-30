import os
import sys
import datetime
import base64
import zlib
import gzip
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import pysam

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
    parser.add_argument('-b', '--bam',
                        dest='BAM',
                        type=str,
                        required=False,
                        help='Input BAM file')
    parser.add_argument('-d', '--depth',
                        dest='DEPTH',
                        type=str,
                        required=False,
                        help='Mosdepth file')
    parser.add_argument('-s', '--read-support',
                        dest='RS',
                        type=str,
                        required=False,
                        help='Read support file')

    args = parser.parse_args(argv)
    return args


def compute_stats_from_bam(bam_file):
    stats = {'n_secondary': 0, 'n_supplementary': 0, 'n_reads': 0,
             'n_unmapped': 0, 'seqlen': 0}

    if bam_file:
        bam = pysam.AlignmentFile(bam_file, "rb")
        for query in bam:
            if query.is_secondary:
                stats['n_secondary'] += 1
                continue
            if query.is_supplementary:
                stats['n_supplementary'] += 1
                continue
            if query.is_unmapped:
                stats['n_unmapped'] += 1

            stats['n_reads'] += 1
            stats['seqlen'] += query.query_length
        bam.close()

    return stats


def get_calls_from_vcf(vcf_file):
    sniffles_file = SvFile(vcf_file)

    calls = []
    for sv in sniffles_file.filtered_variants:
        call = {'id': sv.id,
                'chr1': sv.chr1,
                'pos1': sv.pos1,
                'chr2': sv.chr2,
                'pos2': sv.pos2,
                'sv_type': SvFile.get_display_name(sv.type),
                'sv_length': sv.length,
                'read_support': sv.read_support,
                'precise': sv.precise,
                'gt': sv.gt
                }
        calls.append(call)

    return calls


def get_sv_type_summary(sniffles_file):
    sv_types = {}
    for sv in sniffles_file.filtered_variants:
        display_name = SvFile.get_base_type(sv.type)
        if display_name not in sv_types:
            sv_types[display_name] = 0
        sv_types[display_name] += 1

    return sv_types


def get_summary_from_vcf(vcf_file):
    sniffles_file = SvFile(vcf_file)

    summary = {'sv_types': get_sv_type_summary(sniffles_file)}

    return summary


def compute_depth(mosdepth_file):
    if not mosdepth_file:
        return {'avg_depth': 'NA'}
    try:
        with gzip.open(mosdepth_file, "r") as fh:
            sum_depth = 0
            count_depth = 0
            for line in fh:
                if not line:
                    continue
                cols = line.strip().split(b"\t")
                sum_depth += float(cols[3])
                count_depth += 1
            return {'avg_depth': sum_depth / count_depth }
    except Exception as e:
        return {'avg_depth': 'NA'}


def get_params(readsupport_file):
    if not readsupport_file:
        return {'min_read_support': 'NA'}
    try:
        with open(readsupport_file, "r") as fh:
            rs = int(float(fh.readline().strip()))
            return {'min_read_support': rs }
    except Exception as e:
        return {'min_read_support': 'NA'}


def base_encode64(vcf_file, compress=False):
    with open(vcf_file) as fh:
        if compress:
            encoded = base64.b64encode(zlib.compress(fh.read()).encode()).decode('ascii')
        else:
            encoded = base64.b64encode(fh.read().encode()).decode('ascii')
        return encoded


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
    bam_file = args.BAM
    mosdepth_file = args.DEPTH
    readsupport_file = args.RS
    if not os.path.exists(vcf_file):
        raise OSError("Could not find {}.".format(vcf_file))
    if bam_file and not os.path.exists(bam_file):
        raise OSError("Could not find {}.".format(bam_file))
    if mosdepth_file and not os.path.exists(mosdepth_file):
        raise OSError("Could not find {}.".format(mosdepth_file))
    if readsupport_file and not os.path.exists(readsupport_file):
        raise OSError("Could not find {}.".format(readsupport_file))

    telemetry = Telemetry()

    with telemetry as record:
        record.set('stats', compute_stats_from_bam(bam_file))
        record.set('depth', compute_depth(mosdepth_file))
        record.set('params', get_params(readsupport_file))
        record.set('calls', get_calls_from_vcf(vcf_file))
        record.set('summary', get_summary_from_vcf(vcf_file))
        record.set('vcf', base_encode64(vcf_file))
        record.set('timestamp', datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))

    for call in get_calls_from_vcf(vcf_file):
        print(call, file=sys.stderr)


if __name__ == '__main__':

    main()
