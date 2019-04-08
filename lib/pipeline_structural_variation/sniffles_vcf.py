import random
import re
import pybedtools
import logging

from pysam import VariantFile
from scipy.stats import fisher_exact
from tqdm import tqdm

from pipeline_structural_variation.sv import SvInfo


class SvFile:

    SV_TYPES = ['DEL/INV', 'DUP', 'INV', 'TRA', 'BND', 'INVDUP', 'INS', 'DEL',
                'INV/INVDUP']
    SV_TYPE_NAMES = {'DEL': 'Deletion',
                     'INS': 'Insertion',
                     'INV': 'Inversion',
                     'TRA': 'Translocation',
                     'BND': 'Translocation',
                     'DUP': 'Duplication',
                     'DEL/INV': 'Deletion/Inveserion',
                     'DUP/INS': 'Duplication',
                     'INVDUP': 'Inverted Duplication'}
    SV_BASE_TYPE =  {'DEL': 'DEL',
                     'INS': 'INS',
                     'INV': 'INV',
                     'TRA': 'TRA',
                     'BND': 'BND',
                     'DUP': 'DUP',
                     'DUP/INS': 'DUP',
                     'DEL/INV': 'INV',
                     'INVDUP': 'INV'}

    NO_LENGTH_TYPES = ['BND']

    @staticmethod
    def get_display_name(id):
        if id in SvFile.SV_TYPE_NAMES:
            return SvFile.SV_TYPE_NAMES[id]
        else:
            return id

    @staticmethod
    def get_base_type(sv_type):
        if sv_type in SvFile.SV_BASE_TYPE:
            return SvFile.SV_BASE_TYPE[sv_type]
        else:
            return sv_type

    def read_sniffles_vcf(self, vcf_path, fix_ins=False):
        variants = []
        vcf = VariantFile(vcf_path, "r")
        for rec in vcf.fetch():
            sv_info = SvInfo(rec)

            if fix_ins and sv_info.type in ['INS', 'DUP', 'DUP/INS']:
                logging.debug("Changing {}/{} to {}/{}".format(sv_info.pos1,
                                                                 sv_info.pos2,
                                                                 sv_info.pos1,
                                                                 sv_info.pos1 + 1))
                sv_info.pos2 = sv_info.pos1 + 1
                sv_info.pysam_rec.stop = sv_info.pysam_rec.start + 1

            variants.append(sv_info)

        return variants, vcf.header

    def __init__(self, vcf_path, fix_ins=False):
        self._variants, self.header = self.read_sniffles_vcf(vcf_path,
                                                             fix_ins=fix_ins)
        self.filtered_variants = self._variants

    def check(self):
        for variant in self._variants:
            if variant.pos1 >= variant.pos2 and variant.type != 'BND':
                logging.warning("POS1 >= POS2 for:")
                logging.warning("{}".format(variant))

    def sample(self, n):
        if not self.filtered_variants:
            self.filtered_variants = self._variants
        self.filtered_variants = random.sample(self.filtered_variants, int(n))

    def filter(self, min_length=0, max_length=None, min_support=0,
               max_support=None, include_types=None, chr_filter=None,
               precise=False, gt_min=0, gt_max=None, min_af=0,
               strand_support_p_value=0):
        self.filtered_variants = []

        with tqdm(total=len(self._variants)) as pbar:
            for variant in self._variants:
                pbar.update(1)
                if include_types and variant.type not in include_types:
                    continue
                if min_length and variant.type not in SvFile.NO_LENGTH_TYPES and variant.length <= min_length:
                    continue
                if max_length is not None and variant.length >= max_length:
                    continue
                if min_support and variant.read_support < min_support:
                    continue
                if min_af and variant.af < min_af:
                    continue
                if max_support is not None and variant.read_support > max_support:
                    continue
                if strand_support_p_value:
                    oddsratio, pvalue = fisher_exact(
                        [[variant.ref_support_fwd, variant.ref_support_rev],
                         [variant.read_support_fwd, variant.read_support_rev]])
                    if pvalue < strand_support_p_value:
                        continue

                if chr_filter and not (
                        re.match(chr_filter, variant.chr1, re.I) and re.match(
                        chr_filter, variant.chr2, re.I)):
                    continue
                if precise and not variant.precise:
                    continue

                if gt_max is not None and sum(variant.gt) > gt_max:
                    continue
                if gt_min and sum(variant.gt) < gt_min:
                    continue

                self.filtered_variants.append(variant)

    # def fix_overlaping_ins_dup():
    #     for variant in self._variants:

    def reset(self):
        self.filtered_variants = self._variants

    def n_total(self):
        return len(self._variants)

    def n_filtered(self):
        return len(self.filtered_variants)

    def p_filtered(self):
        if not self.n_total():
            return 0
        return self.n_filtered() * 100.0 / self.n_total()

    def bedtool(self):
        bed_string = ""

        for variant in self.filtered_variants:
            bed_string += variant.to_bed() + '\n'

        return pybedtools.BedTool(bed_string, from_string=True)

    def write_bed(self, path):
        self.bedtool().saveas(path)

    def write_vcf(self, path):
        vcf = VariantFile(path, 'w', header=self.header)
        for variant in self.filtered_variants:
            vcf.write(variant.pysam_rec)
        vcf.close()
