import argparse
from pkg_resources import resource_stream
from Bio import SeqIO


def get_records(sequences):
    records = SeqIO.parse(sequences, "fasta")
    return records


def get_mRNA_from_intervals(sequneces, intervals, start=1):
    pass


def translate_mRNA(mRNA, codon_sun, reading_frame=(1, 2, 3, 4, 5, 6)):
    pass


def main():
    parser = argparse.ArgumentParser(
        description="Translate DNA to amino acid within given intervals"
    )
    parser.add_argument(
        "--sequences","-s", dest="sequences", help="Path to file containing DNA sequences."
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
    print(args)
    # sequences = get_records(args['sequences'])
    # args["sequences"]
    with open(args.codon_code) as handle:
        print(handle.readline())

