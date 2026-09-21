import random
import argparse
from Bio import SeqIO

BASES = ["A", "C", "G", "T"]

def mutate_base(b):
    return random.choice([x for x in BASES if x != b])

def introduce_sub(seq, error_rate):
    seq = list(seq)
    for i, b in enumerate(seq):
        if random.random() < error_rate:
            seq[i] = mutate_base(b)
    return "".join(seq)

def simulate_reads(record, read_len, n_reads, error_rate):
    genome = str(record.seq).upper()
    reads = []
    if len(genome) < read_len:
        return reads
    for i in range(n_reads):
        start = random.randint(0, len(genome) - read_len)
        read = genome[start:start + read_len]
        read = introduce_sub(read, error_rate)
        reads.append(read)
    return reads

def main():
    parser = argparse.ArgumentParser(
        description="Simulate error-controlled reads for false-positive benchmarking"
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA")
    parser.add_argument("-l", "--read_len", type=int, default=1000)
    parser.add_argument("-n", "--reads_per_genome", type=int, default=100)
    parser.add_argument("-e", "--error_rate", type=float, default=0.001)
    parser.add_argument("-s", "--seed", type=int, default=42)
    args = parser.parse_args()
    random.seed(args.seed)
    with open(args.output, "w") as out:
        for record in SeqIO.parse(args.input, "fasta"):
            reads = simulate_reads(record, args.read_len, args.reads_per_genome, args.error_rate)
            for i, read in enumerate(reads):
                header = f">{record.id}_read{i}_{args.read_len}bp"
                out.write(f"{header}\n{read}\n")

if __name__ == "__main__":
    main()
# note: error_rate 0.01 = 1%