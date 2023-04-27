import gzip
from collections import defaultdict
from itertools import product
import matplotlib.pyplot as plt
from Bio import SeqIO
import matplotlib.cm as cm

# Function to count 3-mers in a given sequence
def count_kmers(sequence, k):
    kmers = defaultdict(int)
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if "N" not in kmer:  # Ignore kmers with undefined nucleotides
            kmers[kmer] += 1
    return kmers

# Read the gzipped FASTA file and count 3-mers
filename = "GCA_000715715.1_ASM71571v1_genomic.fna.gz"
k = 3
total_kmers = defaultdict(int)

with gzip.open(filename, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        sequence = str(record.seq).upper()
        kmers = count_kmers(sequence, k)
        for kmer, count in kmers.items():
            total_kmers[kmer] += count

# Generate a histogram of 3-mer frequencies using Matplotlib
all_kmers = ["".join(x) for x in product("ACGT", repeat=k)]
counts = [total_kmers[kmer] for kmer in sorted(all_kmers)]

# Create a colormap with the number of colors equal to the number of bars
cmap = cm.get_cmap("viridis", len(all_kmers))

plt.figure(figsize=(12, 6))
plt.bar(sorted(all_kmers), counts, color=cmap.colors)
plt.xlabel('3-mers')
plt.ylabel('Frequency')
plt.title('Histogram of 3-mer Frequencies')
plt.xticks(rotation=90)
plt.show()
filename = "colorful.png"
plt.savefig(filename, dpi=300, bbox_inches="tight")
