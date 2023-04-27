import requests
import gzip
import io
import random
from Bio import SeqIO
from collections import Counter
import seaborn as sns
import matplotlib.pyplot as plt

# Download the genome file
url = "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/Streptomyces_californicus/latest_assembly_versions/GCA_000715715.1_ASM71571v1/GCA_000715715.1_ASM71571v1_genomic.fna.gz"
response = requests.get(url)

# Read and parse the genome file
with gzip.open(io.BytesIO(response.content), "rt") as f:
    genome_record = list(SeqIO.parse(f, "fasta"))[0]
    genome = str(genome_record.seq)

# Function to count 3-mers in a sequence
def count_kmers(sequence, k=3):
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    return Counter(kmers)

# Count 3-mers in the bacterial genome
bacterial_kmer_counts = count_kmers(genome)

# Generate a random DNA sequence of the same length
nucleotides = ["A", "C", "G", "T"]
random_genome = "".join(random.choices(nucleotides, k=len(genome)))

# Count 3-mers in the random genome
random_kmer_counts = count_kmers(random_genome)

# Function to plot the histograms and save as image files
def plot_histogram(kmer_counts, title, filename):
    sns.barplot(x=list(kmer_counts.keys()), y=list(kmer_counts.values()))
    plt.xlabel("3-mers")
    plt.ylabel("Frequency")
    plt.title(title)
    plt.xticks(rotation=90)
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.show()

# Function to print all 3-mers and their counts
def print_kmers(kmer_counts, title):
    print(title)
    for kmer, count in kmer_counts.items():
        print(f"{kmer}: {count}")

# Print all 3-mers and their counts for both genomes
print_kmers(bacterial_kmer_counts, "Bacterial Genome 3-mers:")
print_kmers(random_kmer_counts, "Random Genome 3-mers:")

# Plot histograms for both genomes and save as image files
plot_histogram(bacterial_kmer_counts, "Bacterial Genome 3-mer Distribution", "bacterial_histogram.png")
plot_histogram(random_kmer_counts, "Random Genome 3-mer Distribution", "random_histogram.png")

# Identify the most common 3-mers
most_common_bacterial_kmers = bacterial_kmer_counts.most_common(5)
print("Most common 3-mers in the bacterial genome:", most_common_bacterial_kmers)

# Biological meaning of the most common 3-mers is beyond the scope of this code solution.

