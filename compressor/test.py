import json
import numpy as np
from Bio import SeqIO
import math
import random

# Load dictionary that maps k-mer to their corresponding index.
# A k-mer and its reverse complement are mapped to the same index.


with open("/home/zhenhao/bucket-map/seed_selection/index/9-mers.json", 'r') as dict_file:
    canonical_kmer_dict = json.load(dict_file)

# We define a utility function here that turns sequences to their 9-mer profiles.

def sequence_to_kmer_profile(sequence : str, k : int = 9, dtype=np.float32):
    """
    Return the k-mer profile of the input sequence (string)
    """
    res = np.zeros(len(set(canonical_kmer_dict.values())), dtype=dtype)
    for i in range(len(sequence) - k + 1):
        k_mer = sequence[i:i + k]
        if k_mer in canonical_kmer_dict:
            res[canonical_kmer_dict[k_mer]] = 1

    return res

def read_buckets_from_file(sequence_file_name : str, bucket_len : int, overlap_len : int):
    """
    Read the DNA sequence and store the sequences in buckets.
    """
    res = []

    for record in SeqIO.parse(sequence_file_name, "fasta"):
        record_sequence = str(record.seq)

        # Split the sequences into buckets
        num_buckets = math.ceil(len(record_sequence) / bucket_len)
        for i in range(num_buckets):
            bucket_sequence = record_sequence[i*bucket_len : (i+1)*bucket_len + overlap_len]
            res.append(bucket_sequence)

    return res

def sample_read_from_sequence(sequence : str, read_len : int):
    sample_range = max(0, len(sequence) - read_len)
    starting_pos = random.randint(0, sample_range)
    return sequence[starting_pos:starting_pos + read_len]

if __name__ == "__main__":
    buckets_1 = read_buckets_from_file("./data/485870.fna")