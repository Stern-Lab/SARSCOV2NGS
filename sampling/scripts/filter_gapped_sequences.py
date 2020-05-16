import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Filter highly gapped sequences from FASTA",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequences", required=True, help="FASTA file of sequences")
    parser.add_argument("--sequences", required=True, help="FASTA file of sequences")
    parser.add_argument("--metadata", type = int, default=50, help="upper boung percentage of gaps")
    parser.add_argument("--output", required=True, help="FASTA file of output file sequence")
    args = parser.parse_args()

    records_to_write = []
    for rec in SeqIO.parse(args.sequences, 'fasta'):
        seq = str(rec.seq)
        n = len(seq)
        seq_N = seq.count('N') / n
        if seq_N < args.precent_N:
            records_to_write.append(rec)
    cnt  = SeqIO.write(records_to_write, args.output, 'fasta')
