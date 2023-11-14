# Downloading genome from the NCBI database

Refer to this [documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/) for details on downloading NCBI genome data package via command line.

```
mkdir tools
cd tools
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets
chmod +x ./datasets
echo "export PATH=\$(pwd):\${PATH}" >> ~/.bashrc
source ~/.bashrc
```

