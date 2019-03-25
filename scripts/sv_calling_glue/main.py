from pysam import VariantFile, AlignmentFile, FastaFile
from scipy.stats import fisher_exact


def load_vcf(vcf_path):
    vcf = VariantFile(vcf_path)
    return vcf


def get_variant_candidates(vcf):
    candidates = []

    for rec in vcf.fetch():
        candidates.append(rec)

    return candidates


def load_fasta(fasta_path):
    fasta = FastaFile(fasta_path)
    return fasta


def load_sam(sam_path):
    sam = AlignmentFile(sam_path, "rb")
    return sam


class VariantCandidate:

    def __init__(self, reference, alignments, variant_call, window=10):
        # Set primary attributes
        self.reference = reference
        self.alignments = alignments
        self.variant_call = variant_call

        # Calculate context window coordinates
        self.window = window
        self.contig, self.start, self.stop = self._get_context_window()

        # Compile reference genome data for context window
        ref_data = self._compile_reference_data_in_context()
        self.reference_string = ref_data[0]
        self.reference_positions = ref_data[1]
        self.homopolymer_kmers = ref_data[2]
        self.gc_percentage = ref_data[3]

        # Compile read data for context window
        self.read_positions = self._compile_read_data_in_context()

        # Carry out aggregations
        self.average_quality_at_call = self._get_average_quality_at_call_position()
        self.percentage_positive_strand_at_call = self._get_strandedness_at_call_position()
        self.call_within_homopolymer = self._get_homopolymer_at_call_position()

    def _get_context_window(self):
        contig = self.variant_call.contig
        start = self.variant_call.pos - self.window
        stop = self.variant_call.pos + self.window

        return contig, start, stop

    # Read data methods
    def _compile_read_data_in_context(self):
        positional_data = {}

        for position in self._yield_alignments_for_positions_in_context():
            current_position = {
                'position_of_call': True if self.variant_call.pos == position.reference_pos else False,
                'reads': [],
                'reference_base': self.reference_positions[position.reference_pos],
            }

            for read in position.pileups:
                current_position['reads'].append(
                    self._compile_read_data_at_position(read)
                )

            current_position['average_quality'] = self._get_average_quality_at_position(
                current_position)
            positional_data[position.reference_pos] = current_position

        return positional_data

    def _yield_alignments_for_positions_in_context(self):
        pileup = self.alignments.pileup(
            self.contig, self.start, self.stop, stepper='all', truncate=True
        )

        for position in pileup:
            yield position

    def _compile_read_data_at_position(self, read):
        read_data = {
            'id': self._get_read_id(read),
            'base': self._get_read_base_at_position(read),
            'strand': self._get_read_strand(read),
            'quality': self._get_read_quality_at_position(read),
            'support': self._get_read_support_at_position(read)
        }

        return read_data

    def _get_read_id(self, read):
        return read.alignment.query_name

    def _get_read_strand(self, read):
        strand = '-' if read.alignment.is_reverse else '+'
        return strand

    def _get_read_base_at_position(self, read):
        base = read.alignment.query_sequence[read.query_position_or_next]
        return base

    def _get_read_quality_at_position(self, read):
        try:
            quality = read.alignment.query_qualities[read.query_position]
        except TypeError:
            quality = None

        return quality

    def _get_read_support_at_position(self, read):
        read_pos = read.alignment.query_sequence[read.query_position_or_next]
        if read_pos == self.variant_call.ref:
            support = 'ref'
        elif read_pos == self.variant_call.alts[0]:
            support = 'alt'
        else:
            support = 'other'
        #         support = 'ref' if read_pos == self.variant_call.ref else 'alt'
        return support

    def _get_average_quality_at_position(self, positional_data):
        quals = [int(read['quality']) for read in positional_data['reads'] if
                 read['quality'] not in [None]]
        if len(quals) == 0:
            return 0.0
        return sum(quals) / len(quals)

    # Reference sequence methods
    def _compile_reference_data_in_context(self):
        ref_seq = self._get_reference_sequence_in_context(self.reference)
        positions = self._get_reference_positional_data(ref_seq)
        homopolymer_kmers = self._detect_reference_homopolymer_kmers(ref_seq)
        gc_percentage = self._get_reference_gc_percentage(ref_seq)

        return ref_seq, positions, homopolymer_kmers, gc_percentage

    def _get_reference_sequence_in_context(self, reference):
        ref_seq = reference.fetch(start=self.start, end=self.stop + 1, region=self.contig)
        return ref_seq

    def _get_reference_positional_data(self, reference_sequence):
        positions = {}
        for idx, pos in enumerate(zip(reference_sequence, range(self.start, self.stop + 1))):
            positions[pos[1]] = {
                'base': pos[0],
                'idx': idx
            }

        return positions

    def _detect_reference_homopolymer_kmers(self, reference_sequence):
        kmers = []
        ex_idx = 0
        iterable = iter(enumerate(reference_sequence))

        for idx, base in iterable:
            remaining_length = len(reference_sequence[idx:-1])

            if not (remaining_length < 2):
                kmer_seed = reference_sequence[idx:idx + 3]
                kmer_bases = set(kmer_seed)

                if len(kmer_bases) <= 1:
                    for ex_idx, extension in enumerate(reference_sequence[idx + 3:]):

                        if extension == kmer_seed[0]:
                            kmer_seed = kmer_seed + extension
                            continue

                        [iterable.__next__() for i in range(2 + ex_idx)]
                        break

                    kmers.append([kmer_seed, idx, idx + 3 + ex_idx])

        return kmers

    def _detect_reference_sequence_motifs(self):
        pass

    def _get_reference_gc_percentage(self, reference_sequence):
        g = reference_sequence.count('G')
        c = reference_sequence.count('C')

        gc_percentage = int((g + c) * 100.0 / len(reference_sequence))

        return gc_percentage

    # Aggregation methods
    def _get_average_quality_at_call_position(self):
        return self.read_positions[self.variant_call.pos]['average_quality']

    def _get_strandedness_at_call_position(self):
        counts = {'ref': {'+': 0, '-': 0}, 'alt': {'+': 0, '-': 0}, 'other': {'+': 0, '-': 0}}
        for read in self.read_positions[self.variant_call.pos - 1]['reads']:
            counts[read['support']][read['strand']] += 1
        oddsratio, pvalue = fisher_exact(
            [[counts['ref']['+'], counts['ref']['-']], [counts['alt']['+'], counts['alt']['-']]])
        return pvalue

    def _get_homopolymer_at_call_position(self):
        for kmer in self.homopolymer_kmers:
            if self.reference_positions[self.variant_call.pos]['idx'] in range(kmer[1], kmer[2]):
                return True

        return False

    # Report methods
    def summarise_context(self):
        print(self.reference_string)

        # Print variant call marker
        snp_marker = ['-' for i in self.reference_string]
        snp_marker[self.reference_positions[self.variant_call.pos]['idx'] - 1] = \
        self.variant_call.alts[0]
        print(''.join(snp_marker))

        # Print homopolymeric tracts
        homopolymer_markers = [' ' for i in self.reference_string]
        for idx, base in enumerate(self.reference_string):
            for kmer in self.homopolymer_kmers:
                if idx in range(kmer[1], kmer[2]):
                    homopolymer_markers[idx] = 'H'
                    continue
        print(''.join(homopolymer_markers))

        # Print low average qual scores
        low_qual_markers = [' ' for i in self.reference_string]

        for k, v in self.read_positions.items():
            if v['average_quality'] < 5:
                low_qual_markers[self.reference_positions[k]['idx']] = 'Q'
            elif v['average_quality'] < 8:
                low_qual_markers[self.reference_positions[k]['idx']] = 'q'
        print(''.join(low_qual_markers))