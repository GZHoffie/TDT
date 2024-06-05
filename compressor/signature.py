import pickle
from compressor.vectorization import Seq2KMers
import numpy as np
from Bio import SeqIO
from tqdm import tqdm
import torch
import torch.nn as nn
from compressor.dataset import RandomReadDataset, ReadDataset
from torch.utils.data import DataLoader

class Signature:
    def __init__(self) -> None:
        pass

    def find_consensus(self, **kwargs):
        pass

    def to_signature(self, sequences):
        pass

class FMHSignature(Signature):
    def __init__(self, condition_function, seed_length, window_length=None, read_from_file=None) -> None:
        super().__init__()
        self.con = condition_function
        self.tokenizer = Seq2KMers(seed_length, window_length)
        if read_from_file is None:
            self.kmer_dict, self.kmer_num = self.tokenizer.init_vocab(self.con)
        else:
            self.load_vocab(read_from_file)
            
    

    def to_signature(self, sequences, data_type=np.float32):
        signature = np.zeros(self.kmer_num, dtype=data_type)
        for sequence in sequences:
            kmers = self.tokenizer.canonical_kmers(sequence)
            for kmer in kmers:
                if kmer in self.kmer_dict:
                    signature[self.kmer_dict[kmer]] = 1

        #print(signature, signature.sum())
        return signature
    
    def store_vocab(self, file_name):
        with open(file_name, 'wb') as f:
            pickle.dump(self.kmer_dict, f)
    

    def load_vocab(self, file_name):
        with open(file_name, 'rb') as f:
            self.kmer_dict = pickle.load(f)
        
        self.kmer_num = len(set(self.kmer_dict.values()))
    

    def find_consensus(self, **kwargs):
        """
        Assuming the name of file is the taxid.
        """
        from pathlib import Path

        #print(files)

        consensus_signature = np.zeros(self.kmer_num, dtype=np.int16)
        #signature_storage = [] # FOR DEBUG

        files = kwargs["files"]
        #print(files)

        # Attempt: simply sum up all the signatures, take
        # the ones that are shared among 80% of the species
        for f in tqdm(files):
            taxid = int(Path(f).stem)
            sequences = []
            for record in SeqIO.parse(f, "fasta"):
                sequences.append(str(record.seq))
            
            signature = self.to_signature(sequences, data_type=np.int16)
            consensus_signature += signature
            #signature_storage.append(np.array(signature, dtype=bool)) # FOR DEBUG
        
        #print(genus_signature, genus_signature.sum())
        

        return consensus_signature#, signature_storage


class LogisticRegression(nn.Module):
    # build the constructor
    def __init__(self, n_inputs, n_outputs):
        super(LogisticRegression, self).__init__()
        self.l1 = torch.nn.Linear(n_inputs, n_outputs)

    # make predictions
    def forward(self, x):
        x = self.l1(x)
        x = torch.sigmoid(x)
        #print(x)
        x = torch.flatten(x)
        #print(x)

        return x
            
class MLSignature(FMHSignature):
    def __init__(self, condition_function, seed_length, window_length=None, read_from_file=None) -> None:
        super().__init__(condition_function, seed_length, window_length, read_from_file)
    

    def _init_model(self):
        # Initialize a machine learning model. Here we use binary logistic regression
        model = LogisticRegression(self.kmer_num, 1)
        return model
    
    def find_consensus(self, **kwargs):
        positive_samples_training = kwargs["positive_samples_training"]
        negative_samples_training = kwargs["negative_samples_training"]
        positive_samples_test = kwargs["positive_samples_test"]
        negative_samples_test = kwargs["negative_samples_test"]

        # initialize model
        classifier = self._init_model()
        optimizer = torch.optim.Adam(classifier.parameters(), lr=0.001)
        criterion = torch.nn.BCELoss()

        # Load the files to the dataset
        dataset_training = RandomReadDataset(positive_samples_training, negative_samples_training, self, num_reads_per_epoch=10000, random_errors=[0.03, 0.01, 0.01])
        dataloader_training = DataLoader(dataset_training, batch_size=32, shuffle=False)

        dataset_test = RandomReadDataset(positive_samples_test, negative_samples_test, self, num_reads_per_epoch=2000, positive_rate=0.5)
        dataloader_test = DataLoader(dataset_test, batch_size=32, shuffle=False)#, collate_fn=custom_collate_fn)
        
        epochs = 50

        Loss_training = []
        Acc_training = []
        Sensitivity_training = []
        Loss_val = []
        Acc_val = []
        Sensitivity_val = []

        for epoch in range(epochs):
            classifier.train()
            #metric.reset()
            train_loss = 0.0
            train_correct_predictions = 0
            train_correct_virus_id = 0
            train_positive_samples = 0
            for i, (reads, labels) in enumerate(dataloader_training):
                
                optimizer.zero_grad()
                outputs = classifier(reads)
                #metric.update(outputs, labels)
                #print(outputs, labels)
                loss = criterion(outputs, labels)

                # Record the performance
                train_loss += loss.item()
                predicted = (outputs > 0.5).float()
                train_correct_predictions += (predicted == labels).sum()
                train_correct_virus_id += (predicted.bool() & labels.bool()).sum()
                train_positive_samples += labels.sum()

                # Training
                loss.backward()
                optimizer.step()

            Loss_training.append(train_loss / len(dataset_training))
            Acc_training.append(100 * (train_correct_predictions.item()) / len(dataset_training))
            Sensitivity_training.append(100 * train_correct_virus_id.item() / train_positive_samples)
            #AUROC_training.append(metric.compute())

            

            # validation
            classifier.eval()
            #metric.reset()
            val_loss = 0.0
            val_correct_predictions = 0
            val_correct_virus_id = 0
            val_positive_samples = 0


            with torch.no_grad():
                for i, (reads, labels) in enumerate(dataloader_test):
                    #print(reads.size(), labels.size())
                    #print(reads, labels)
                    outputs = classifier(reads)
                    #print(outputs)
                    #print(labels)
                    #print(outputs, labels)
                    loss = criterion(outputs, labels)

                    # Record the performance
                    val_loss += loss.item()
                    predicted = (outputs > 0.5).float()
                    val_correct_predictions += (predicted == labels).sum()
                    val_correct_virus_id += (predicted.bool() & labels.bool()).sum()
                    val_positive_samples += labels.sum()
            
            Loss_val.append(val_loss / len(dataset_test))
            Acc_val.append(100 * (val_correct_predictions.item()) / len(dataset_test))
            Sensitivity_val.append(100 * val_correct_virus_id.item() / val_positive_samples)
            
            print('Epoch: {}. Train Loss: {}. Accuracy: {}. Sensitivity: {}. Val Loss: {}. Accuracy: {}. Sensitivity: {}.'.format(epoch, Loss_training[-1], Acc_training[-1], Sensitivity_training[-1], Loss_val[-1], Acc_val[-1], Sensitivity_val[-1]))



            

            
    

if __name__ == "__main__":
    import glob

    def fracMinHash(kmer_hash):
        hash = (976369 * kmer_hash + 1982627) % 10000
        if hash < 1000:
            return True
        else:
            return False
    
    def all(kmer_hash):
        return True
        
    sg = MLSignature(all, 12)#, read_from_file="12-mer.pkl")
    sg.store_vocab("12-mer.pkl")
    #sg.find_consensus(positive_samples_training=glob.glob("./data/escherichia/*.fna"),
    #                  negative_samples_training=glob.glob("./data/staphylococcus/*.fna"),
    #                  positive_samples_test=glob.glob("./data/escherichia/*.fna"),
    #                  negative_samples_test=glob.glob("./data/staphylococcus/*.fna"))
    #print("562")
    #sg.insert_all_sequences_in_file("./data/562.fna")
    #print("564")
    #sg._find_consensus_in_genus(glob.glob("./data/escherichia/*.fna"))
    #print(len(sg.graph.nodes))
    #print("485870")
    #sg.insert_all_sequences_in_file("./data/485870.fna")

    #print(len(sg.graph.nodes))