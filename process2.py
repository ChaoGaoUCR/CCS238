import gzip
import random
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

# Read the gzipped FASTA file and get the total genome length
filename = "GCA_000715715.1_ASM71571v1_genomic.fna.gz"
total_genome_length = 0

with gzip.open(filename, "rt") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        total_genome_length += len(record.seq)

# Generate a random DNA genome of the same length
k = 3
nucleotides = ['A', 'C', 'G', 'T']
random_genome = ''.join(random.choices(nucleotides, k=total_genome_length))

# Count 3-mers in the random genome
random_kmers = count_kmers(random_genome, k)

# Generate a histogram of 3-mer frequencies using Matplotlib
all_kmers = ["".join(x) for x in product("ACGT", repeat=k)]
counts = [random_kmers[kmer] for kmer in sorted(all_kmers)]

# Create a colormap with the number of colors equal to the number of bars
cmap = cm.get_cmap("plasma", len(all_kmers))

plt.figure(figsize=(12, 6))
plt.bar(sorted(all_kmers), counts, color=cmap.colors)
plt.xlabel('3-mers')
plt.ylabel('Frequency')
plt.title('Histogram of 3-mer Frequencies in Random Genome')
plt.xticks(rotation=90)
filename = "colorful2.png"
plt.savefig(filename, dpi=300, bbox_inches="tight")
plt.show()
