import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_sequences(sequences):
    return SeqIO.parse(sequences, "fasta")


def get_intervals(intervals, sep="\t", col="chr", dtype={'start': int, 'stop': int}):
    return pd.read_csv(intervals, sep=sep, index_col=col)


def get_codon_mapping(codon_code, header=True):
    codon_mapping = {}
    with open(codon_code) as handle:

        if header:  # a header in the file will be ignored
            handle.readline()

        for line in handle:
            split = line.split()
            codon_mapping[split[0]] = split[1]

    return codon_mapping


def get_mRNA_from_intervals(sequences, intervals, index_start=0):
    mRNAs = []

    for sequence in sequences:
        interval = (intervals[intervals.index == sequence.id])

        for _, row in interval.iterrows():
            start = row['start'] + index_start
            stop = row['stop'] + 1

            mRNA = SeqRecord(
                    Seq(sequence.seq[start:stop].reverse_complement()),
                    id=sequence.id,
                    name='>' + row['id']
                    )
            mRNAs.append(mRNA)  

    return mRNAs


def main():
# if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Translate DNA to amino acid within given intervals"
    )
    parser.add_argument(
        "--sequences",
        "-s",
        dest="sequences",
        help="Path to file containing DNA sequences.",
    )
    parser.add_argument(
        "--intervals",
        "-i",
        dest="intervals",
        help="Path to file containing mRNA coordinates",
    )
    parser.add_argument(
        "--codon_code",
        "-c",
        dest="codon_code",
        help="Path to file that maps codons to amino acids.",
    )
    args = parser.parse_args()

    sequences = get_sequences(args.sequences)
    intervals = get_intervals(args.intervals)
    codon_mapping = get_codon_mapping(args.codon_code)

    mRNAs = get_mRNA_from_intervals(sequences, intervals)
