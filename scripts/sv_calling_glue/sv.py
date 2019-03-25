import logging
import re


class SvInfo:

    def __init__(self, pysam_vcf_rec):
        self.id = pysam_vcf_rec.id
        self.type = pysam_vcf_rec.info.get('SVTYPE', None)
        self.length = abs(pysam_vcf_rec.info.get('SVLEN', None))
        self.read_support = pysam_vcf_rec.info.get('RE', None)

        self.chr1 = pysam_vcf_rec.chrom
        self.chr2 = pysam_vcf_rec.info.get('CHR2', None)

        self.precise = pysam_vcf_rec.info.get('PRECISE')
        self.pos1 = pysam_vcf_rec.pos
        self.pos2 = pysam_vcf_rec.stop

        self.af = pysam_vcf_rec.info.get('AF', None)[0] if 'AF' in pysam_vcf_rec.info else 1.0

        if self.type.upper() == "BND":
            if len(pysam_vcf_rec.alts) > 1:
                raise ValueError("No support for multiple alt alleles: {}".format(pysam_vcf_rec))
            breakend_alt = pysam_vcf_rec.alts[0]
            m = re.search('.*[\[\]](.*):(.*)[\]\[].*', breakend_alt)
            self.chr2 = m.group(1)
            self.pos2 = m.group(2)

        if len(pysam_vcf_rec.samples.keys()) != 1:
            raise RuntimeError(
                "Currently only a single sample per file is supported")

        self.gt = pysam_vcf_rec.samples.values()[0]['GT']

        # Stranded read support
        self.read_support_fwd = None
        self.read_support_rev = None
        if 'STRANDS2' in pysam_vcf_rec.info:
            read_strands = pysam_vcf_rec.info.get('STRANDS2', None)
            # print(read_strands)
            self.read_support_fwd = read_strands[0]
            self.read_support_rev = read_strands[1]

        # Stranded ref support
        self.ref_support_fwd = None
        self.ref_support_rev = None
        ref_strands = pysam_vcf_rec.info.get('REF_strand', None)
        if ref_strands:
            # print(ref_strands)
            self.ref_support_fwd = ref_strands[0]
            self.ref_support_rev = ref_strands[1]

        self.pysam_rec = pysam_vcf_rec

    def to_bed(self):
        return "{}\t{}\t{}\t{}\t{}".format(self.chr1, self.pos1, self.pos2,
                                           self.id, self.read_support)

    def __repr__(self):
        return "{} {} {} {} {} {} {} {} {}".format(self.id, self.chr1,
                                                   self.pos1, self.chr2,
                                                   self.pos2, self.type,
                                                   self.length,
                                                   self.read_support,
                                                   self.precise)
