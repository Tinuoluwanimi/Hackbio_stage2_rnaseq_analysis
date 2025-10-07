#!/bin/bash
# Create raw data directory
RAW_DIR="plants/data/raw"
mkdir -p $RAW_DIR
cd $RAW_DIR
# Download raw FASTQ files
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/097/SRR12808497/SRR12808497.fastq.gz -o SRR12808497.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/099/SRR12808499/SRR12808499.fastq.gz -o SRR12808499.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/029/SRR12808529/SRR12808529.fastq.gz -o SRR12808529.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/028/SRR12808528/SRR12808528.fastq.gz -o SRR12808528.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/027/SRR12808527/SRR12808527.fastq.gz -o SRR12808527.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/098/SRR12808498/SRR12808498.fastq.gz -o SRR12808498.fastq.gz

echo "All FASTQ files downloaded into $RAW_DIR "

