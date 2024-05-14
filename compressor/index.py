import compressor.utils
from pathlib import Path
import subprocess
import os
import glob
from tqdm import tqdm
import pandas as pd


class SignatureIndex:
    def __init__(self, condition_function, k=12, p_value=0.05, kmer_file=None) -> None:
        self.signature = compressor.utils.FMHSignature(condition_function, k, read_from_file=kmer_file)
    
    def _download_reference(self, accession, directory, taxid):
        """
        Download the reference using the given accession number, and store it in
        ${directory}/${taxid}.fna. If the file already exists, we append to the file.
        """
        # mkdir if the directory doesn't exist
        Path(directory).mkdir(parents=True, exist_ok=True)

        # Change to this directory
        os.chdir(directory)

        # Download using NCBI datasets
        subprocess.run(["datasets", "download", "genome", "accession",
                         accession, "--filename", taxid + ".zip"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if os.path.exists(taxid + ".zip"):
            subprocess.run(["unzip", taxid + ".zip"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Store the results in ${directory}/${taxid}.fna
        files_name = glob.glob("ncbi_dataset/data/*/*.fna")
        with open(taxid + ".fna", "wb") as f:
            subprocess.run(["cat"] + files_name, stdout=f)
        subprocess.run(["rm", "-r", "ncbi_dataset", taxid + ".zip", "README.md"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def _remove_directory(self, directory):
        """
        Delete the directory and all the content inside.
        """
        subprocess.run(["rm", "-r", directory], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


    def find_signature(self, accession_list, taxid_list, genus_id):
        """
        Given a list of accession, determine the signature.
        """
        # Download the reference corresponding to the accession
        directory = "/home/zhenhao/TDT/data_temp/"

        for accession, taxid in tqdm(zip(accession_list, taxid_list), desc="Downloading references for genus " + str(genus_id)):
            self._download_reference(accession, directory, str(taxid))

        res = self.signature._find_consensus_in_genus(glob.glob(directory + "*.fna"))
        self._remove_directory(directory)

        return res

    def find_all_signatures(self, metadata_df, num_samples=None):
        """
        Given a metadata_df with columns "species_taxid", "genus_taxid", "ncbi_genbank_assembly_accession"
        find the signature for all genera.

        For each genera, only take `num_samples` species to speed things up.
        """
        res = {}
        genus_set = set(metadata_df.dropna(subset=['genus_taxid'])['genus_taxid'])
        for genus in genus_set:
            # sample just one genome per species
            genus_df = metadata_df[metadata_df["genus_taxid"] == genus].groupby("species_taxid").sample(1)

            # if there are more than `num_samples` species, sample only `num_samples` of them
            if num_samples is not None and len(genus_df) > num_samples:
                genus_df = genus_df.sample(num_samples)
            
            # Download the references
            signature = self.find_signature(genus_df["ncbi_genbank_assembly_accession"],
                                            genus_df["species_taxid"], genus)
            res[genus] = signature
        return res
            

        


if __name__ == "__main__":
    metadata_df = pd.read_csv("./gtdb_utils/metadata.csv")
    def fracMinHash(kmer_hash):
        hash = (976369 * kmer_hash + 1982627) % 10000
        if hash <= 10:
            return True
        else:
            return False
    
    si = SignatureIndex(fracMinHash, 12, kmer_file="/home/zhenhao/TDT/12-mer.pkl")
    res = si.find_all_signatures(metadata_df, num_samples=100)

    import json
    with open("signature.json", 'w') as f:
        json.dump(res, f)

        

