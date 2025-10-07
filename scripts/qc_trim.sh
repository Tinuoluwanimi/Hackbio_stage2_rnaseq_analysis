#!/bin/bash
# QC + trimming for raw files

#Define directories
RAW_DIR="plants/data/raw"
QC_RAW="plants/data/qc_raw"
TRIM_DIR="plants/data/trimmed"
QC_TRIM="plants/data/qc_trimmed"
MULTIQC_RAW="plants/data/multiqc_raw"
MULTIQC_TRIM="plants/data/multiqc_trimmed"

mkdir -p $QC_RAW $TRIM_DIR $QC_TRIM $MULTIQC_RAW $MULTIQC_TRIM

# QC on raw reads
fastqc -o $QC_RAW ${RAW_DIR}/*.fastq.gz
echo "MultiQC on raw reads" 
multiqc $QC_RAW -o $MULTIQC_RAW

# Trimming with fastp 
echo "Trimming with fastp"
for f in ${RAW_DIR}/*.fastq.gz; do
    base=$(basename $f .fastq.gz)
    out=${TRIM_DIR}/${base}_trim.fastq.gz
    json=${TRIM_DIR}/${base}_fastp.json
    html=${TRIM_DIR}/${base}_fastp.html

    fastp -i $f -o $out -j $json -h $html
done

# QC on trimmed reads
fastqc -o $QC_TRIM ${TRIM_DIR}/*.fastq.gz
echo "MultiQC on trimmed reads"
multiqc $QC_TRIM -o $MULTIQC_TRIM

echo "QC + trimming complete "
