mkdir -p data
cd data
datasets download genome taxon 562 --reference --filename pseudomonas_dataset.zip
unzip pseudomonas_dataset.zip
cat ncbi_dataset/data/*/*.fna >> 562.fna
rm -r ncbi_dataset pseudomonas_dataset.zip README.md