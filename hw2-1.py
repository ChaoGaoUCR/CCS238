import networkx as nx
import argparse

def read_fasta(file):
    """Reads a fasta file and returns a list of sequences."""
    with open(file) as f:
        reads = f.read().split('>')[1:]  # Ignore the empty string before the first '>'
        reads = [read.split('\n', 1)[1].replace('\n', '') for read in reads]
    return reads

def overlap(a, b, min_length=3):
    """Return length of longest suffix of 'a' overlapping
    with prefix of 'b' that is at least 'min_length'
    characters long.  If no such overlap exists,
    return 0."""
    start = 0  # start all the way at the left
    while True:
        start = a.find(b[:min_length], start)  # look for b's prefix in a
        if start == -1:  # no more occurrences to right
            return 0
        # found occurrence; check for full suffix/prefix match
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  # move just past previous match

def create_graph(reads, min_length):
    G = nx.DiGraph()
    for read in reads:
        for read2 in reads:
            if read != read2:
                len_overlap = overlap(read, read2, min_length)
                if len_overlap > 0:
                    G.add_edge(read, read2, weight=len_overlap)
    return G

def assemble_genome(G):
    genome = ""
    current_read = list(G.nodes)[0]
    while len(G) > 0:
        genome += current_read
        if len(G[current_read]) > 0:  # If there are outgoing edges from the current read
            next_read = max(G[current_read], key=lambda x: G[current_read][x]['weight'])  # Choose the neighbor with the highest 'weight'
            genome = genome[:-G[current_read][next_read]['weight']]  # Avoid duplicating the overlap
            G.remove_node(current_read)
            current_read = next_read
        else:  # If there are no outgoing edges from the current read
            G.remove_node(current_read)
            if len(G) > 0:
                current_read = list(G.nodes)[0]
    return genome

def main(input_file):
    reads = read_fasta(input_file)
    G = create_graph(reads, min_length=3)
    genome = assemble_genome(G)
    print(genome)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Assemble a genome from a set of reads.')
    parser.add_argument('input_file', type=str, help='Input fasta file with the reads.')
    args = parser.parse_args()
    main(args.input_file)
