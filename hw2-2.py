import networkx as nx
import argparse

def read_fasta(file):
    """Reads a fasta file and returns a list of sequences."""
    with open(file) as f:
        reads = f.read().split('>')[1:]  # Ignore the empty string before the first '>'
        reads = [read.split('\n', 1)[1].replace('\n', '') for read in reads]
    return reads

def create_debruijn_graph(reads, k):
    G = nx.DiGraph()
    for read in reads:
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]
            G.add_edge(kmer[:-1], kmer[1:], weight=len(kmer))
    return G

def assemble_genome(G):
    genome = ""
    path = list(nx.eulerian_path(G))
    genome = path[0][0]
    for (_,v) in path:
        genome += v[-1]
    return genome

def main(input_file):
    reads = read_fasta(input_file)
    G = create_debruijn_graph(reads, k=3)
    genome = assemble_genome(G)
    print(genome)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Assemble a genome from a set of reads.')
    parser.add_argument('input_file', type=str, help='Input fasta file with the reads.')
    args = parser.parse_args()
    main(args.input_file)
