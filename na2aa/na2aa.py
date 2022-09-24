import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def get_sequences(sequences):
    return SeqIO.parse(sequences, "fasta")


def get_intervals(intervals, sep="\t", col="chr"):
    return pd.read_csv(intervals, sep=sep, index_col=col)


def get_codon_mapping(codon_code, header=True):
    """
    Create a dictionary for codon mapping from input file.
    Stop codons should be defined as '*'.
    """
    codon_mapping = {}
    with open(codon_code) as handle:

        if header:  # a header in the file will be ignored
            handle.readline()

        for line in handle:
            split = line.split()
            codon_mapping[split[0]] = split[1]

    return codon_mapping


def get_mRNA_from_intervals(sequences, intervals, index_start=0, index_stop=1):
    """
    Get the mRNA sequences with the interval file and return them as a list.
    index_start and index_stop can be used for including, excluding the start and stop
    or to take 0 and 1-based indices into account.
    """
    mRNAs = []

    for sequence in sequences:
        interval = intervals[intervals.index == sequence.id]

        for _, row in interval.iterrows():
            start = row["start"] + index_start
            stop = row["stop"] + index_stop

            mRNAs.append(
                SeqRecord(
                    Seq(sequence.seq[start:stop]), 
                    id=sequence.id, 
                    name=">" + row["id"]
                )
            )
            # add reverse complement
            mRNAs.append(
                SeqRecord(
                    Seq(sequence.seq[start:stop].reverse_complement()),
                    id=sequence.id,
                    name=">" + row["id"],
                )
            )

    return mRNAs


def get_longest_translations_from_mRNA(mRNAs, codon_mapping):
    """
    Takes a list of mRNAs and returns a dictionary with the longest
    translation for each mRNA.
    Takes into account all six reading frames
    """

    # create dict to store final results
    longest_translations = {seq.name: "" for seq in mRNAs}

    for mRNA in mRNAs:
        for i in range(0, 3):
            translated_mRNA = translate_mRNA(mRNA.seq[i:], codon_mapping)
            longest_seq = get_longest_aa_sequence(translated_mRNA)
            if len(longest_seq) > len(longest_translations[mRNA.name]):
                longest_translations[mRNA.name] = longest_seq

    return longest_translations


def translate_mRNA(sequence, codon_mapping):
    """
    Use codon_mapping to translate mRNA to amino acid sequence.
    Sequence will be returned with stop codons not breaking
    the sequence and without taking start codons into account.
    """
    aa_sequence = ""

    n = len(sequence)

    # from https://github.com/biopython/biopython/blob/master/Bio/Seq.py
    for i in range(
        0, n - n % 3, 3
    ):
        codon = sequence[i : i + 3]
        aa_sequence += codon_mapping[codon]

    return aa_sequence


def get_longest_aa_sequence(aa_sequence, start_amino_acid="M"):
    """
    Takes as input an unmodified amino acid sequence as produced by translate_mRNA()
    Amino acid sequences start with a methionine and end with a stop codon.
    """
    split = aa_sequence.split("*")

    trimmed_sequences = []
    for sequence in split:
        start = sequence.find(start_amino_acid)
        if start != -1:
            trimmed_sequences.append(sequence[start:])
    if trimmed_sequences:
        return max(trimmed_sequences, key=len)
    else:
        return ""


def main():
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

    # load files
    sequences = get_sequences(args.sequences)
    intervals = get_intervals(args.intervals)
    codon_mapping = get_codon_mapping(args.codon_code)

    # extract sequences of interest
    mRNAs = get_mRNA_from_intervals(sequences, intervals)

    # store longest translation products for each interval in a dict
    longest_translations = get_longest_translations_from_mRNA(mRNAs, codon_mapping)

    # print the results to console
    for gene in longest_translations.keys():
        print(gene)
        print(longest_translations[gene])
