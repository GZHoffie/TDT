mkdir -p data
cd data
taxid=561
datasets download genome taxon $taxid --reference --filename ${taxid}.zip
if [ -f ${taxid}.zip ]; then
    unzip ${taxid}.zip
    cat ncbi_dataset/data/*/*.fna > ${taxid}.fna
    rm -r ncbi_dataset ${taxid}.zip README.md
fi

accession=GCA_000190495.1
datasets download genome accession $accession --reference --filename ${accession}.zip
if [ -f ${accession}.zip ]; then
    unzip ${accession}.zip
    cat ncbi_dataset/data/*/*.fna >> ${accession}.fna
    rm -r ncbi_dataset ${accession}.zip README.md
fi
