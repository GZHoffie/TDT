import sys
sys.path.append("/home/zhenhao/virus-classification")

import torch
import json
import math
import numpy as np
from tqdm import tqdm
import random
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from compressor.utils import insert_errors, reverse_complement
from Bio import SeqIO


class RandomReadDataset(Dataset):
    def __init__(self, positive_sample_files, negative_sample_files, vectorization_method, positive_rate=0.5, read_len=300, 
                 num_reads_per_epoch=1000, random_errors=[0.01, 0.01, 0.01], random_reverse_complement=False, debug=False):
        """
        Randomly sample reads from the reference genomes.

        Args:
            - `positive_sample_files`: a list of strings, storing the filenames that contains positive samples.
            - `negative_sample_files`: a list of strings, storing the filenames that contains negative samples.
            - `vectorization_method`: an `online_vc.vectorization.Seq2Vec` class, specifying the way of transforming a read into
              feature vector.
            - `positive_rate`: the proportion of the reads from the positive samples.
            - `read_len`: the length of each read.
            - `num_reads_per_epoch`: number of reads trained per epoch.
            - `random_errors`: None (no data augmentation) or a list of 3 floats, corresponding to substituion, insertion and deletion rate
            - `random_reverse_complement`: if set to True, with probability 1/2, return the reverse complement of the read.
        """
        self.positive_samples = self._load_reference(positive_sample_files)
        self.negative_samples = self._load_reference(negative_sample_files)

        self.seq2vec = vectorization_method
        
        self.pos_rate = positive_rate
        self.read_len = read_len

        self.num_reads_per_epoch = num_reads_per_epoch

        # Data augmentation options
        self.errors = random_errors
        self.rev_comp = random_reverse_complement

        self.debug = debug
    

    def _load_reference(self, reference_list):
        res = []
        for file in reference_list:
            if file.endswith('.fna') or file.endswith('.fasta') or file.endswith('.fa'):
                file_type = "fasta"
            else:
                file_type = "fastq"
            for record in SeqIO.parse(file, file_type):
                res.append(record.seq)
        return res

    def _sample_read(self, is_positive):
        if is_positive:
            # sample from self.positive_samples
            sequence_index = random.randrange(len(self.positive_samples))
            sample_range = max(0, len(self.positive_samples[sequence_index]) - self.read_len)
            starting_pos = random.randint(0, sample_range)
            return self.positive_samples[sequence_index][starting_pos:starting_pos + self.read_len]
        else:
            # sample from self.positive_samples
            sequence_index = random.randrange(len(self.negative_samples))
            sample_range = max(0, len(self.negative_samples[sequence_index]) - self.read_len)
            starting_pos = random.randint(0, sample_range)
            return self.negative_samples[sequence_index][starting_pos:starting_pos + self.read_len]


    
    def __len__(self):
        return self.num_reads_per_epoch

    def __getitem__(self, idx):
        # Randomly choose positive/negative sample
        is_positive = (random.random() <= self.pos_rate) # if is_positive=1, then choose a positive sample; else choose a negative sample

        # Randomly choose a sequence
        sampled_sequence = str(self._sample_read(is_positive))

        # Add data augmentation
        if self.errors is not None:
            sampled_sequence = insert_errors(sampled_sequence, self.errors[0], self.errors[1], self.errors[2])
        
        if self.rev_comp and random.random() > 0.5:
            sampled_sequence = reverse_complement(sampled_sequence)


        vectorized_sequence = torch.tensor(self.seq2vec.to_signature(sampled_sequence))
        if self.debug:
            print(vectorized_sequence, is_positive)
        return vectorized_sequence, torch.tensor(is_positive, dtype=torch.float32)
    
class ReadDataset(Dataset):
    def __init__(self, positive_sample_files, negative_sample_files, vectorization_method, 
                random_errors=None, random_reverse_complement=False):
        """
        Load all reads in the fasta files.

        Args:
            - `positive_sample_files`: a list of strings, storing the filenames that contains positive samples.
            - `negative_sample_files`: a list of strings, storing the filenames that contains negative samples.
            - `vectorization_method`: an `online_vc.vectorization.Seq2Vec` class, specifying the way of transforming a read into
              feature vector.
            - `positive_rate`: the proportion of the reads from the positive samples.
            - `random_errors`: None (no data augmentation) or a list of 3 floats, corresponding to substituion, insertion and deletion rate
            - `random_reverse_complement`: if set to True, with probability 1/2, return the reverse complement of the read.
        """
        positive_samples = self._load_reference(positive_sample_files) 
        negative_samples = self._load_reference(negative_sample_files)
        self.samples = positive_samples + negative_samples
        self.labels = [True] * len(positive_samples) + [False] * len(negative_samples)

        self.seq2vec = vectorization_method

        # Data augmentation options
        self.errors = random_errors
        self.rev_comp = random_reverse_complement

    

    def _load_reference(self, reference_list):
        res = []
        for file in reference_list:
            if file.endswith('.fna') or file.endswith('.fasta') or file.endswith('.fa'):
                file_type = "fasta"
            else:
                file_type = "fastq"
            for record in SeqIO.parse(file, file_type):
                res.append(record.seq)
        return res


    
    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx):
        sequence = self.samples[idx]
        # Add data augmentation
        if self.errors is not None:
            sequence = insert_errors(sequence, self.errors[0], self.errors[1], self.errors[2])
        
        if self.rev_comp and random.random() > 0.5:
            sequence = reverse_complement(sequence)

        vectorized_sequence = torch.tensor(self.seq2vec.to_signature(str(sequence)))

        return vectorized_sequence, torch.tensor(self.labels[idx], dtype=torch.float32)
    
