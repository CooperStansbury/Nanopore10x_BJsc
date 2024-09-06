#!/bin/bash

# Set input and output file paths
data_dir="/nfs/turbo/umms-indikar/shared/projects/cell_cycle/data/sc_cell_cycle_2024/barcodes"
biolegend_fasta="$data_dir/biolegend.fasta" 
cell_id_fasta="$data_dir/cell_ids.fasta" 
merged_barcodes_fastq="$data_dir/merged_barcodes.fastq.gz"
output_dir="$data_dir/cutadapt"
json_file="$output_dir/cutadapt.json"
info_file="$output_dir/info.tsv"
summary_file="$output_dir/summary.txt"

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# Run cutadapt
cutadapt -e 1 -j 24 --revcomp \
    -b file:$biolegend_fasta \
    --info-file=$info_file \
    --discard-untrimmed \
    --report='minimal' \
    --action='none' \
    --json=$json_file \
    -o "$output_dir/{name}.fastq.gz" $merged_barcodes_fastq > $summary_file