#!/bin/bash
# Download Arabidopsis genome and annotation, build STAR index

cd ~/Tinuoluwanimi/plants/references

# Download Arabidopsis thaliana genome FASTA and annotation GTF
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-59/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-59/gtf/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.59.gtf.gz

# Unzip files
gunzip -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
gunzip -f Arabidopsis_thaliana.TAIR10.59.gtf.gz

# Build STAR genome index
mkdir -p star_index

STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir ./star_index \
  --genomeFastaFiles Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  --sjdbGTFfile Arabidopsis_thaliana.TAIR10.59.gtf \
  --sjdbOverhang 100

echo "STAR genome index built successfully "
