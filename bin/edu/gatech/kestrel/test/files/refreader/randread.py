#!/usr/bin/python3

"""
Generate random FASTQ reads.
"""

import argparse
import random
import sys

from abc import ABCMeta, abstractmethod


class QualityGenerator(ABCMeta):

    @abstractmethod
    def get_quality_string(self):
        pass


class UniformQualityGenerator(metaclass=QualityGenerator):

    def __init__(self, params):
        tok = params.split(':')

        if len(tok) > 2:
            raise ValueError('Uniform generator: Received more than 2 parameter: {}'.format(params))

        if len(tok) == 1:
            tok[1] = tok[0]

        self.low = int(tok[0]) + 33
        self.high = int(tok[1]) + 33

    def get_quality_string(self, length):
        return ''.join(chr(random.choice(range(self.low, self.high))) for i in range(int(length)))

class FastaGenerator:

    def __init__(self, seq_generator):
        self.seq_generator = seq_generator

    def generate(self, num_sequence, read_len, start_index=1):
        for seq_index in range(start_index, start_index + num_sequence):
            yield '>Seq-{}\n{}'.format(seq_index, self.seq_generator.next_sequence(read_len))

class FastqGenerator:

    def __init__(self, seq_generator, quality_generator):
        self.seq_generator = seq_generator
        self.quality_generator = quality_generator

    def generate(self, num_sequence, read_len, start_index=1):
        for seq_index in range(start_index, start_index + num_sequence):
            yield '@Seq-{}\n{}\n+Seq-{}\n{}'.format(
                seq_index,
                self.seq_generator.next_sequence(read_len),
                seq_index,
                self.quality_generator.get_quality_string(read_len))

class RawGenerator:

    def __init__(self, seq_generator):
        self.seq_generator = seq_generator

    def generate(self, num_sequence, read_len, start_index=1):
        for seq_index in range(start_index, start_index + num_sequence):
            yield self.seq_generator.next_sequence(read_len)

class SequenceGenerator:

    def __init__(self, prob_base=1.0, seq_rna=False, all_iupac=False, mix_case=False):
        """
        Create a new sequence generator.

        :param prob_base: Probability of A, C, G, or T.
        :param n_only: Generate only N (not all IUPAC characters) for non-ACGT bases if prob_base is less than 1.0.
        """
        self.prob_base = prob_base
        self.all_iupac = all_iupac

        if self.prob_base < 0.0:
            self.prob_base = 0.0

        elif self.prob_base > 1.0:
            self.prob_base = 1.0

        # Set bases
        self.base = ['A', 'C', 'G', 'T']

        if seq_rna:
            self.base[3] = 'U'

        # Set other characters
        if all_iupac:
            self.other = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N', '.', '-']

        else:
            self.other = ['N']

        if mix_case:
            add_bases = list()

            for seq_char in self.base:
                add_bases.append(seq_char.lower())

            self.base += add_bases

            add_bases = list()

            for seq_char in self.other:
                if seq_char.isalpha():
                    add_bases.append(seq_char.lower())

            self.other += add_bases

        return

    def next_sequence(self, num_bases):
        return ''.join(self.next_base(num_bases))

    def next_base(self, num_bases):
        """
        Generate the next ``len`` bases.

        :param len: Number of bases to generate.

        :return: A generator for ``len`` bases.
        """
        for i in range(num_bases):
            rand = random.random()

            if rand < self.prob_base:
                yield self.base[int((self.prob_base - rand) / self.prob_base * len(self.base))]
            else:
                yield self.other[int((rand - self.prob_base) / (1 - self.prob_base) * len(self.other))]

def main(args=None):

    # Check arguments
    if args is None:
        args = sys.argv[1:]

    # Get parsed arguments
    parser = argparse.ArgumentParser(description='Generate random FASTQ reads.', prog='randread')

    parser.add_argument('-d', '--dist', dest='score_dist', type=str, default='u:0:93',
                        help='Distribution of FASTQ quality scores (default = u:0:93).')

    parser.add_argument('-f', '--format', dest='format', type=str, default='fastq',
                        help='Format of the file to generate (default = fastq).')

    parser.add_argument('-l', '--readlen', dest='read_len', type=int, default=100,
                        help='Read length (default = 100).')

    parser.add_argument('-m', '--mixcase', dest='mix_case', default=False, action='store_true',
                        help='Generate sequences with a mix of upper and lower case characters.')

    parser.add_argument('-n', '--numreads', dest='num_reads', type=int, default=1,
                        help='Number of reads to generate (default = 1).')

    parser.add_argument('-p', '--pbase', dest='prob_base', type=float, default=1.0,
                        help='The probability that a base in a sequence is A, C, G, or T (default = 1.0).')

    parser.add_argument('-r', '--rna', dest='seq_rna', default=False, action='store_true',
                        help='Generate U in place of T')

    parser.add_argument('-s', '--seed', dest='seed', type=int, default=None,
                        help='Set random seed, which will result in the same sequence being generated for the '
                             'same seed. By default, generated sequences are random on each run.')

    parser.add_argument('-u', '--iupac', dest='all_iupac', action='store_true', default=False,
                        help='If the probability of generating A, C, G, or T is less than 1.0, the default behavior is '
                             'to output N for all non-ACGT characters. This option will output any IUPAC nucleotide '
                             'character instead of just N for each non-ACGT character.')

    parsed_args = parser.parse_args(args=args)

    # Set base probability
    seq_generator = SequenceGenerator(parsed_args.prob_base, parsed_args.seq_rna, parsed_args.all_iupac, parsed_args.mix_case)

    # Set random seed
    random.seed(parsed_args.seed)

    # Check parsed arguments
    if parsed_args.num_reads < 1:
        raise ValueError('Cannot generate less than 1 read: {}'.format(parsed_args.num_reads))

    if parsed_args.read_len < 1:
        raise ValueError('Cannot generate reads with less than 1 base: {}'.format(parsed_args.read_len))

    # Determine distribution
    tok = parsed_args.score_dist.split(':', 1)

    tok[0] = tok[0].strip()

    if len(tok) == 1:
        tok[1] = ''
    else:
        tok[1] = tok[1].strip()

    # Get generator
    if tok[0] == 'u':
        quality_generator = UniformQualityGenerator(tok[1])

    else:
        raise ValueError('Cannot find generator for type: {}'.format(tok[0]))

    # Get sequence generator
    format = parsed_args.format.strip().lower()

    if format == 'fasta':
        seq_generator = FastaGenerator(seq_generator)
        join_str = '\n'

    elif format == 'fastq':
        seq_generator = FastqGenerator(seq_generator, quality_generator)
        join_str = '\n'

    elif format == 'raw':
        seq_generator = RawGenerator(seq_generator)
        join_str = '\n\n'

    else:
        raise ValueError('Unrecognized sequence format: {}: Must be fasta or fastq'.format(parsed_args.format))

    # Generate sequences
    print(join_str.join(seq_generator.generate(parsed_args.num_reads, parsed_args.read_len)))
    # for seq_index in range(1, parsed_args.num_reads + 1):
    #     print('@Seq-{}'.format(seq_index))
    #     print(generator.next_sequence(parsed_args.read_len))
    #     print('+')
    #     print(quality_generator.get_quality_string(parsed_args.read_len))

if __name__ == '__main__':
    main()
