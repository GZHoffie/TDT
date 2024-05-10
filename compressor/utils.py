
import heapdict
from tqdm import tqdm
import json
import os
from queue import PriorityQueue
import networkx as nx
from Bio import SeqIO


class Seq2Tokens:
    """
    An implementation of the Minimizer algorithm

    Args:
        vocab (:obj:`int`, `optional`):
            A dictionnary of k-mer keys and their ids :obj:`{AAA (0): 0,...}`

        unk_token (:obj:`str`, `optional`):
            The unknown token to be used by the model.
    """

    def __init__(self, seed_length, window_length, use_mask=True):
        # Initialize parameters for the minimizer.
        self._k = seed_length
        self._l = window_length

        if (self._l < self._k):
            print(f"[ERROR]\t\tthe window length (currently {self._l}) must be equal or larger than the seed length (currently {self._k}).")
            return
        else:
            print(f"[INFO]\t\tSeed length set to {self._k}, estimated vocabulary size {int(4 ** self._k / 2)}.")

        # Utils to find minimizers
        self._nucleotide_to_hash = {'A': 0, 'a': 0,
                                    'C': 1, 'c': 1,
                                    'G': 2, 'g': 2,
                                    'T': 3, 't': 3}
        self._seed_mask = int("1" * 2 * self._k, 2)
        # Random permuatation of hash values, same as SeqAn3
        self._hash_mask = 0x8F3F73B5CF1C9ADE if use_mask else 0


    def _hash(self, sequence):
        """
        Convert a DNA sequence to a sequence of numbers in 0-3.
        """
        #TODO: handle ambiguous characters 'N', 'n' properly by deleting those k-mers. 
        #TODO: Current colustions: simply ignore the ambiguous characters.
        sequence_hash = [self._nucleotide_to_hash[n] for n in sequence if n in self._nucleotide_to_hash]
        return sequence_hash
    

    def _rev_comp(self, sequence_hash):
        """
        Given a sequence of hash values returned from self._hash(sequence),
        return the reverse complement of that sequence.
        """
        rev_comp_hash = [(3 - h) for h in reversed(sequence_hash)]
        return rev_comp_hash


    def _kmers(self, sequence_hash):
        """
        Convert k-mers to their hash values.
        """
        hash_values = []
        if len(sequence_hash) - self._k + 1 <= 0:
            return hash_values
        
        # find the hash value of the first k-mer
        current_hash = 0
        for i in range(self._k):
            current_hash = current_hash << 2
            current_hash = current_hash | sequence_hash[i]
        
        hash_values.append(current_hash)
        
        # find the hash of rest of the k-mers
        for i in range(len(sequence_hash) - self._k):
            current_hash = (current_hash << 2) & self._seed_mask
            current_hash = current_hash | sequence_hash[i + self._k]
            hash_values.append(current_hash)
        
        return hash_values

    def canonical_kmers(self, sequence):
        forward_hash = self._hash(sequence)
        reverse_hash = self._rev_comp(forward_hash)

        forward_kmers = self._kmers(forward_hash)
        reverse_kmers = self._kmers(reverse_hash)

        min_kmers = [min(i, j) for i, j in zip(forward_kmers, reversed(reverse_kmers))]
        return min_kmers


    def minimizers(self, sequence):
        """
        Find the minimizers of the sequence.
        """
        minimizers = []
        if len(sequence) - self._k + 1 <= 0:
            return minimizers
        
        # Find hash values of k-mers
        forward_hash = self._hash(sequence)
        reverse_hash = self._rev_comp(forward_hash)

        masked_forward_kmers = [h ^ self._hash_mask for h in self._kmers(forward_hash)]
        masked_reverse_kmers = [h ^ self._hash_mask for h in self._kmers(reverse_hash)]

        # find canonical k-mers
        masked_min_kmers = [min(i, j) for i, j in zip(masked_forward_kmers, reversed(masked_reverse_kmers))]

        # Find min in sliding window
        sliding_window = heapdict.heapdict()
        for i in range(min(self._l - self._k + 1, len(masked_min_kmers))):
            sliding_window[i] = masked_min_kmers[i]

        if len(sliding_window) != 0:
            minimizers.append(sliding_window.peekitem()[1])

        for i in range(len(masked_min_kmers) - self._l):
            sliding_window.pop(i)
            current_index = i + self._l - self._k + 1
            sliding_window[current_index] = masked_min_kmers[current_index]
            if len(sliding_window) != 0:
                minimizer = sliding_window.peekitem()[1]
                if minimizers[-1] != minimizer:
                    minimizers.append(minimizer)

        unmasked_minimizers = [m ^ self._hash_mask for m in minimizers]
        tokens = [str(m) for m in unmasked_minimizers]
        #print(tokens)

        return tokens
    

class Seq2Graph:

    def __init__(self, condition_function, seed_length, window_length) -> None:
        self.con = condition_function
        self.tokenizer = Seq2Tokens(seed_length, window_length)
        self.graph = nx.MultiDiGraph()
    

    def _insert_sequence(self, sequence, sequence_id):
        kmers = self.tokenizer.canonical_kmers(sequence)
        filtered_kmers = [i for i in kmers if self.con(i)]
        edges = []
        for i in range(len(filtered_kmers) - 1):
            edges.append((filtered_kmers[i], filtered_kmers[i+1], dict(id=sequence_id)))
        print(len(kmers), len(edges))
        
        dup_edges = 0
        for e in edges:
            if self.graph.has_edge(e[0], e[1]):
                dup_edges += 1
            
        print(dup_edges)
        
        self.graph.add_edges_from(edges)
            
    

    def insert_all_sequences_in_file(self, file_name):

        for record in SeqIO.parse(file_name, "fasta"):
            acc = record.id.split(" ")[0]
            self._insert_sequence(str(record.seq), acc)
            print(len(sg.graph.nodes))
    

if __name__ == "__main__":
    def fracMinHash(kmer_hash):
        hash = (976369 * kmer_hash + 1982627) % 10000
        if hash <= 10:
            return True
        else:
            return False
        
    sg = Seq2Graph(fracMinHash, 12, 12)
    sg.insert_all_sequences_in_file("./data/485870.fna")
    #print(len(sg.graph.nodes))
    
    sg.insert_all_sequences_in_file("./data/2745495.fna")

    #print(len(sg.graph.nodes))


