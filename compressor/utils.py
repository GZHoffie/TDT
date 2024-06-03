
import heapdict
from tqdm import tqdm
import json
import os
from queue import PriorityQueue
import numpy as np
from Bio import SeqIO
import random

dna_nucleotides = ['A', 'C', 'G', 'T']
nucleotide_to_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

def insert_errors(read : str, substitute_rate : float, insert_rate : float, delete_rate : float):
    """
    Randomly add errors into the read.

    Args:
        - `read`: the original DNA sequence. It should only contain 'A', 'C', 'G', 'T', or 'N'.
        - `substitution_rate`: the rate of substitutions to be added. the number of substitution errors
          would be read length * substitution rate.
        - `insert_rate`, `delete_rate`: indel error rate.
    """
    num_substitutions = int(len(read) * substitute_rate)
    num_insert = int(len(read) * insert_rate)
    num_delete = int(len(read) * delete_rate)

    res = read.upper()

    # Add deletions
    for _ in range(num_delete):
        i = random.randrange(len(res))
        res = res[:i] + res[i+1:]

    # Add insertions
    for _ in range(num_insert):
        i = random.randrange(len(res) + 1)
        new_nt = random.choice(dna_nucleotides)
        res = res[:i] + new_nt + res[i:]

    # Add substitutions
    for _ in range(num_substitutions):
        i = random.randrange(len(res))
        orig_nt = res[i]
        remaining_nt = [n for n in dna_nucleotides if n != orig_nt]
        res = res[:i] + random.choice(remaining_nt) + res[i+1:]
    
    return res
    


def create_reads_with_error(fasta_file : str, new_file : str, substitute_rate : float, insert_rate : float, delete_rate : float):
    """
    Given a file of reads, create a new file with the same reads with added errors.
    """
    with open(fasta_file, 'r') as f_in:
        lines = f_in.readlines()
    
    with open(new_file, 'w') as f_out:
        for line in lines:
            if line.startswith('>'):
                f_out.write(line)
            else:
                read = line.strip()
                new_read = insert_errors(read, substitute_rate, insert_rate, delete_rate)
                f_out.write(new_read + '\n')


def reverse_complement(read : str):
    """
    Return the reverse complement of a read.
    """
    res = "".join([nucleotide_to_complement[n] for n in reversed(read.upper()) if n in nucleotide_to_complement])
    return res
